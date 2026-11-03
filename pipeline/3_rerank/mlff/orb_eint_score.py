"""Stage E3 -- Orb interaction energy for CYP3A4 poses: E_int = E(complex) - E(receptor) - E(substrate).

RECEPTOR = protein pocket + HEM.  SUBSTRATE = the drug only.  That sentence is the whole job: a
heme left on the substrate side turns E_int into "drug+heme vs apo protein", a large and entirely
plausible number that crashes nothing.  The prep stage already asserted the split; this re-asserts
it from the masks rather than trusting the file.

WHY ORB AND NOT EGRET-1.  Egret-1's embedding table is ten elements (H C N O F P S Cl Br I) in
every EGRET_1*.model checkpoint and raises "26 is not in list" on iron -- it failed 74/80 on an
Fe-bearing set where Orb failed 0/80.  Orb orbmol_v2 additionally PASSES the four falsification
tests in orb_fe_falsification.py (crystal heme inside the 20 eV/A gate at 11.09; Fe-N minimum at
the crystal 2.064 A; Fe out-of-plane minimum at 0.0 A; spin spread 1.36 eV).  That detects gross
failure only -- it does not establish that iron is described WELL.

TWO SPIN CONVENTIONS, BOTH REPORTED, NEITHER CHOSEN HERE.  Substrate-bound ferric P450 heme is
high spin (mult 6), but an azole/pyridine N coordinating the iron drives it LOW (mult 2); 62% of
this manifest coordinates.  So for a coordinating pose:

    E_int_crossover   = E(complex, mult 2) - E(receptor, mult 6) - E(substrate, 1)
    E_int_conserving  = E(complex, mult 6) - E(receptor, mult 6) - E(substrate, 1)

The first includes the spin-state change the ligand causes -- physically complete for Type II
binding, but it loads a large model-dependent spin gap into the ranking, and Orb's accuracy on
spin gaps is untested.  The second is a cleaner interaction energy that omits the very effect
making those ligands potent.  Non-coordinating poses have the two equal by construction.

THE ASE CACHE.  atoms.get_potential_energy() routes through check_state(), which compares
positions, numbers, cell and pbc -- NOT atoms.info, where Orb reads charge and spin.  At fixed
geometry a changed spin returns the STALE energy; this was measured, and it silently reported one
multiplicity three times in the falsification run.  Every single point here goes through
calc.reset() + calc.calculate(..., all_changes), which is what egret_energy.py:27 does.
analysis/egret_ranking/orb_energy.py:46 uses get_potential_energy and carries the bug.

FORCE GATE IS A FLAG, NEVER A DROP (force_gate_pass).  max|F_component| over the COMPLEX against
20 eV/A.  A strained heme region can fail the gate for reasons unrelated to the pose being ranked,
and an unmodified crystal heme already sits at 11.09, so the headroom is thin.  Dropping on it
would silently delete Type II poses.

FAILURES ARE ROWS, NOT ABSENCES.  "an empty energy table and a finished one are indistinguishable
downstream once the pooling stage imputes NaNs".

IO IS STDLIB csv, NOT pandas.  orb-env has no pandas (the existing analysis/egret_ranking/
orb_energy.py uses csv for the same reason).  Installing it into a working MLIP env to save a few
lines is not worth the risk of perturbing the torch/cuda stack this model depends on.

Env: orb-env.  GPU (CPU works, ~50x slower).
"""
from __future__ import annotations

import argparse
import os
import time
from pathlib import Path

os.environ.setdefault("TORCHDYNAMO_DISABLE", "1")   # cluster g++ is c++2a -> CppCompileError

import csv

import numpy as np

ROOT = Path(os.environ.get("COFOLD_ROOT", "/project/bme-20255165-cofolding-physics"))
FORCE_GATE = 20.0
SPIN_SUBSTRATE, SPIN_RECEPTOR_HIGH = 1, 6


def _parity_ok(z, charge, mult) -> bool:
    """multiplicity 2S+1 must match the electron-count parity."""
    nelec = int(np.sum(z)) - int(charge)
    return (nelec % 2) == ((int(mult) - 1) % 2)


def build(sym, xyz, charge, spin):
    import ase
    a = ase.Atoms(symbols=list(sym), positions=xyz)
    a.set_cell(np.ptp(xyz, axis=0) + 50.0)   # Orb needs det(cell) != 0 even non-periodic
    a.set_pbc(False)
    a.center()
    a.info["charge"] = int(charge)
    a.info["spin"] = int(spin)
    return a


def single_point(calc, sym, xyz, charge, spin, want_forces=False):
    """Energy and max|F_component|; with want_forces, also the (N,3) force array.

    `want_forces` is a keyword with a False default so every existing caller keeps its arity and
    its behaviour: the forces are ALREADY requested on every single point below, so returning them
    costs nothing and omitting them costs nothing either.
    """
    from ase.calculators.calculator import all_changes
    a = build(sym, xyz, charge, spin)
    calc.reset()                                     # see THE ASE CACHE above
    calc.calculate(a, ["energy", "forces"], all_changes)
    r = calc.results
    e, fmax = float(r["energy"]), float(np.abs(r["forces"]).max())
    if want_forces:
        return e, fmax, np.asarray(r["forces"], dtype=np.float64)
    return e, fmax


IFACE_CUT = 4.5          # A; a receptor atom this close to any ligand atom is "interface"


def interface_forces(dump_dir, tag, Z, F_cx, F_sub, F_rec, lig, rec_mask):
    """Per-atom INTERACTION forces across the protein-ligand interface.

        F_int = F(complex) - F(receptor alone) - F(substrate alone)

    with each subsystem's forces scattered back into complex index space. F(complex) on its own is
    the wrong object for this question: its intramolecular part is two to three orders of
    magnitude larger than the interaction term, so its argmax reports whichever bond inside the
    ligand is most strained, not which contact ACROSS the interface is. The same subtraction is
    what model/mlff/backends/mace.py:188-193 does, for the same reason.

    Note what `max_force` in the CSV is, by contrast: max|F_component| over every atom AND every
    Cartesian axis of the whole 12 A crop. It is not a per-atom magnitude and it is not restricted
    to the interface, so it cannot answer "which atom is most strained" -- that is why this exists.

    Returns scalar summary columns and writes the full arrays to <dump_dir>/<tag>.npz.
    """
    coords = Z["coords"].astype(np.float64)
    sym = Z["symbols"].astype(str)
    res = Z["res_name"].astype(str)
    F_ref = np.zeros_like(F_cx)
    F_ref[lig] = F_sub
    F_ref[rec_mask] = F_rec
    F_int = F_cx - F_ref

    lig_idx, rec_idx = np.flatnonzero(lig), np.flatnonzero(rec_mask)
    # (n_lig, n_rec); n_lig is ~50 and n_rec ~3000, so this is a few hundred thousand floats.
    D = np.linalg.norm(coords[lig_idx][:, None, :] - coords[rec_idx][None, :, :], axis=-1)
    d_rec = D.min(axis=0)                       # each receptor atom's distance to the ligand
    iface_rec = rec_idx[d_rec <= IFACE_CUT]
    # res_name carries no residue number, so adjacent residues of the same type would be
    # indistinguishable; a run-length ordinal over the crop disambiguates them.
    res_ord = np.cumsum(np.r_[True, res[1:] != res[:-1]]) - 1

    mag_int, mag_cx = np.linalg.norm(F_int, axis=1), np.linalg.norm(F_cx, axis=1)

    def label(i):
        return f"{sym[i]}:{res[i]}#{int(res_ord[i])}"

    out = {}
    # ligand side: the ligand atom pulled hardest by the receptor, and who is pulling it
    il = lig_idx[int(np.argmax(mag_int[lig_idx]))]
    jr = rec_idx[int(np.argmin(D[int(np.flatnonzero(lig_idx == il)[0])]))]
    out.update(fint_lig_max=float(mag_int[il]), fint_lig_atom=label(il),
               fint_lig_partner=label(jr),
               fint_lig_pair_d_A=float(np.linalg.norm(coords[il] - coords[jr])))
    # receptor side: the interface receptor atom pulled hardest by the ligand, and by which atom
    if len(iface_rec):
        jr2 = iface_rec[int(np.argmax(mag_int[iface_rec]))]
        il2 = lig_idx[int(np.argmin(D[:, int(np.flatnonzero(rec_idx == jr2)[0])]))]
        out.update(fint_rec_max=float(mag_int[jr2]), fint_rec_atom=label(jr2),
                   fint_rec_partner=label(il2),
                   fint_rec_pair_d_A=float(np.linalg.norm(coords[jr2] - coords[il2])))
    iface = lig.copy()
    iface[iface_rec] = True
    out.update(fint_iface_max=float(mag_int[iface].max()), n_iface=int(iface.sum()),
               fmag_cx_max=float(mag_cx.max()), fmag_cx_lig_max=float(mag_cx[lig_idx].max()),
               # translational invariance: both must sum to ~0 over all atoms. A scatter bug
               # breaks this long before it produces an implausible-looking maximum.
               fsum_cx_abs=float(np.abs(F_cx.sum(axis=0)).max()),
               fsum_int_abs=float(np.abs(F_int.sum(axis=0)).max()))

    d_full = np.zeros(len(sym))          # 0 on the ligand itself, distance-to-ligand elsewhere
    d_full[rec_idx] = d_rec
    dump_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(dump_dir / f"{tag}.npz", F_complex=F_cx.astype(np.float32),
                        F_int=F_int.astype(np.float32), coords=coords.astype(np.float32),
                        symbols=sym, res_name=res, res_ord=res_ord.astype(np.int32),
                        is_ligand=lig, is_interface=iface,
                        d_to_ligand=d_full.astype(np.float32))
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", type=Path,
                    default=ROOT / "data/cyp3a4_relax/eint_manifest.csv")
    ap.add_argument("--prep-dir", type=Path,
                    default=Path("/scratch-shared/20255165/cyp3a4_relax/eint_prep"))
    ap.add_argument("--out-csv", type=Path, required=True)
    ap.add_argument("--model", default="orbmol_v2")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0)
    # Opt-in, default off: this module is run by at least four campaigns' sbatch stages
    # (cyp3a4_bioemu2 60/60b, cyp3a4_eintval 13/13b), and unconditional extra columns or extra
    # I/O would change their outputs for a question they never asked.
    ap.add_argument("--dump-forces", type=Path, default=None,
                    help="write per-atom F_complex/F_int npz here and add interface-force columns")
    a = ap.parse_args(argv)

    import torch
    from orb_models.forcefield import pretrained
    from orb_models.forcefield.inference.calculator import ORBCalculator

    with open(a.manifest, newline="") as fh:
        M = list(csv.DictReader(fh))
    M = M[a.shard::a.nshards]
    if a.limit:
        M = M[:a.limit]
    device = "cuda" if torch.cuda.is_available() else "cpu"
    orbff, adapter = getattr(pretrained, a.model)(device=device, precision="float32-high")
    calc = ORBCalculator(orbff, atoms_adapter=adapter, device=device)
    print(f"[orb] device={device} model={a.model}  shard {a.shard}/{a.nshards}  "
          f"{len(M)} poses", flush=True)

    rows, t0 = [], time.time()
    truthy = {"true", "1", "1.0", "True"}
    for n, r in enumerate(M, 1):
        coord = str(r["fe_coordinating"]).strip() in truthy
        rec = {"tag": r["tag"], "arm": r["arm"], "query_id": r["query_id"],
               "seed": int(r["seed"]), "sample": int(r["sample"]), "base_id": r["base_id"],
               "min_ipae": float(r.get("min_ipae") or "nan"), "fe_coordinating": coord, "status": ""}
        try:
            npz = a.prep_dir / f"{r['tag']}.npz"
            Z = np.load(npz, allow_pickle=False)
            sym, xyz, z = Z["symbols"].astype(str), Z["coords"].astype(np.float64), Z["z"]
            lig = Z["is_ligand"].astype(bool)
            # re-assert the split from the masks; the one error that crashes nothing
            if bool((Z["is_hem"].astype(bool) & lig).any()) or bool((Z["is_fe"].astype(bool) & lig).any()):
                raise ValueError("HEM/Fe flagged is_ligand -- receptor/substrate split inverted")
            if not lig.any() or lig.all():
                raise ValueError(f"degenerate substrate mask: {int(lig.sum())} of {len(lig)}")
            rec_mask = ~lig
            q_sub = int(Z["ligand_charge"])
            q_rec = int(Z["pocket_charge"]) + int(Z["hem_charge"])
            q_cx = int(Z["total_charge"])
            if q_rec + q_sub != q_cx:
                raise ValueError(f"charges do not sum: receptor {q_rec} + substrate {q_sub} "
                                 f"!= complex {q_cx}")
            spin_cx = int(r["spin_complex"])
            rec.update(n_atoms=int(len(sym)), n_substrate=int(lig.sum()),
                       charge_complex=q_cx, charge_receptor=q_rec, charge_substrate=q_sub,
                       spin_complex=spin_cx, pocket_cut_A=float(Z["pocket_cut"]),
                       parity_ok=bool(_parity_ok(z, q_cx, spin_cx)))

            want = a.dump_forces is not None
            e_sub, _, *fs = single_point(calc, sym[lig], xyz[lig], q_sub,
                                         SPIN_SUBSTRATE, want_forces=want)
            e_rec, _, *fr = single_point(calc, sym[rec_mask], xyz[rec_mask], q_rec,
                                         SPIN_RECEPTOR_HIGH, want_forces=want)
            e_cx, f_cx, *fc = single_point(calc, sym, xyz, q_cx, spin_cx, want_forces=want)
            rec.update(E_substrate=e_sub, E_receptor=e_rec, E_complex=e_cx,
                       max_force=f_cx, force_gate_eV_A=FORCE_GATE,
                       force_gate_pass=int(f_cx <= FORCE_GATE),
                       E_int_crossover=e_cx - e_rec - e_sub)
            if want:
                rec.update(interface_forces(a.dump_forces, r["tag"], Z,
                                            fc[0], fs[0], fr[0], lig, rec_mask))
            if spin_cx != SPIN_RECEPTOR_HIGH:
                e_cx6, f_cx6 = single_point(calc, sym, xyz, q_cx, SPIN_RECEPTOR_HIGH)
                rec.update(E_complex_high=e_cx6, max_force_high=f_cx6,
                           E_int_conserving=e_cx6 - e_rec - e_sub)
            else:                      # non-coordinating: the two conventions coincide
                rec.update(E_complex_high=e_cx, max_force_high=f_cx,
                           E_int_conserving=e_cx - e_rec - e_sub)
            rec["status"] = "ok"
        except Exception as e:  # noqa: BLE001 -- a row, never an absence
            rec["status"] = f"{type(e).__name__}: {e}"[:200]
            print(f"  FAIL {r['tag']}: {rec['status']}", flush=True)
        rows.append(rec)
        if n % 10 == 0:
            ok = sum(1 for x in rows if x["status"] == "ok")
            mem = (torch.cuda.max_memory_allocated() / 2**30) if device == "cuda" else 0.0
            print(f"  {n}/{len(M)}  ok={ok}  {time.time()-t0:.0f}s  peakGPU={mem:.2f} GiB",
                  flush=True)

    cols, seen = [], set()
    for rec in rows:                       # union of keys, first-seen order; failures lack most
        for k in rec:
            if k not in seen:
                seen.add(k); cols.append(k)
    a.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for rec in rows:
            w.writerow(rec)

    ok = [x for x in rows if x["status"] == "ok"]
    print(f"\nshard {a.shard}: {len(ok)}/{len(rows)} scored  ({time.time()-t0:.0f}s)")
    bad = {}
    for x in rows:
        if x["status"] != "ok":
            bad[x["status"][:80]] = bad.get(x["status"][:80], 0) + 1
    for k, v in sorted(bad.items(), key=lambda kv: -kv[1])[:5]:
        print(f"  {v:4d}  {k}")
    if ok:
        cx = np.array([x["E_int_crossover"] for x in ok])
        cs = np.array([x["E_int_conserving"] for x in ok])
        mf = np.array([x["max_force"] for x in ok])
        cd = np.array([x["fe_coordinating"] for x in ok], dtype=bool)
        print(f"E_int_crossover  : median {np.median(cx):.3f} eV  min {cx.min():.3f}  "
              f"max {cx.max():.3f}")
        print(f"E_int_conserving : median {np.median(cs):.3f} eV")
        if cd.any():
            print(f"|crossover - conserving| on the {int(cd.sum())} coordinating poses: "
                  f"median {np.median(np.abs(cx - cs)[cd]):.3f} eV  "
                  f"max {np.abs(cx - cs)[cd].max():.3f} eV")
        print(f"force gate pass  : {sum(x['force_gate_pass'] for x in ok)}/{len(ok)}  "
              f"(max_force median {np.median(mf):.2f} eV/A)")
        npar = sum(1 for x in ok if not x.get("parity_ok", True))
        if npar:
            print(f"SPIN PARITY MISMATCH on {npar} poses -- that multiplicity is impossible for "
                  "the electron count; the charge accounting is wrong")
    print(f"wrote {a.out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
