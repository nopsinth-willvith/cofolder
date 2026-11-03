# Vendored from $COFOLD_ROOT/analysis/cyp3a4_bioemu2/mmgbsa_modes.py -- unchanged.
"""Single-snapshot MM-GBSA for BioEmu2 mode medoids.

The campaign's ``mmgbsa_clusters.py`` averages over frames of an MD trajectory in
the system's own tleap topology. BioEmu2 modes have no such trajectory -- each
mode contributes one medoid pose, freshly parameterised by ``build_system.py`` --
so this scores that single snapshot instead.

That is the weaker of the two framings and is reported as such, but it is the one
``mmgbsa_clusters.py``'s own docstring already endorses: "more defensible as a
RELATIVE ranking within one ligand's own modes, which is exactly the use here...
Comparing across ligands is NOT supported by this calculation."

Same engine and settings as the GaMD stage: ante-MMPBSA splits on ``:LIG`` with
``--radii=mbondi2``, then MMPBSA.py with igb=2 and 0.15 M salt.
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
from pathlib import Path

MMGBSA_IN = "&general\n  startframe=1, endframe=1, interval=1, verbose=2,\n/\n&gb\n  igb=2, saltcon=0.15,\n/\n"

# MM-GBSA on an UNMINIMISED generative pose scores clash severity, not binding:
# the pilot returned +16.6 and +27.0 kcal/mol on two of three modes, against the
# GaMD campaign's -42.8..-4.5 range from MD-equilibrated frames. build_system
# MMFF-repairs the ligand's INTERNAL geometry only and never relaxes the
# receptor-ligand contacts, so the force field has to see the complex first.
# Backbone restrained so the pose is relieved, not re-docked.
MIN_IN = """minimise before scoring
&cntrl
  imin=1, maxcyc=5000, ncyc=2500, ntmin=1,
  ntb=1, ntp=0, cut=9.0, ntc=1, ntf=1,
  ntpr=500, ntwx=0, ntwr=0,
  ntr=1, restraint_wt=2.0, restraintmask='@CA,C,N,O',
/
"""


def run(cmd, cwd, what):
    p = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"{what} rc={p.returncode}: {(p.stderr or p.stdout)[-300:]}")
    return p


def parse_delta(dat: Path) -> dict:
    """Pull the DELTA (complex - receptor - ligand) block, not the Complex absolutes.

    Anchoring on 'Differences' rather than on 'GENERALIZED BORN' is the fix the
    campaign already made once: the earlier parser reported Complex-block
    absolutes as vdw/eel while delta_total stayed correct.
    """
    text = dat.read_text()
    i = text.find("Differences")
    if i < 0:
        raise RuntimeError("no Differences block in MMPBSA output")
    block, out = text[i:], {}
    for key, label in (("vdw", "VDWAALS"), ("eel", "EEL"), ("egb", "EGB"),
                       ("esurf", "ESURF"), ("epb", "EPB"), ("enpolar", "ENPOLAR"),
                       ("edisper", "EDISPER"), ("delta_total", "DELTA TOTAL")):
        m = re.search(rf"^{re.escape(label)}\s+(-?\d+\.\d+)", block, re.M)
        if m:
            out[key] = float(m.group(1))
    if "delta_total" not in out:
        m = re.search(r"DELTA\s+TOTAL\s+(-?\d+\.\d+)", block)
        if m:
            out["delta_total"] = float(m.group(1))
    return out


def lig_shift(prm, crd0, crd1, wd) -> float | str:
    """Ligand heavy-atom RMSD between the input and minimised coordinates."""
    try:
        (wd / "rms.in").write_text(
            f"parm {prm}\ntrajin {crd1}\nreference {crd0}\n"
            f"rms LIG :LIG&!@H= reference out {wd}/rms.dat\ngo\nquit\n")
        run(["cpptraj", "-i", str(wd / "rms.in")], wd, "cpptraj_rms")
        vals = [float(l.split()[1]) for l in (wd / "rms.dat").read_text().splitlines()
                if l.strip() and not l.startswith("#")]
        return round(vals[0], 3) if vals else ""
    except Exception:  # noqa: BLE001
        return ""


def fe_lig_min(prm, crd) -> float | str:
    """Shortest heme Fe -- ligand heavy-atom distance in one frame.

    Logged before and after minimisation because a --no-fe-bond build carries the Fe-N contact
    by LJ + electrostatics alone, so the minimiser is free to pull a Type-II donor off the iron.
    """
    try:
        import numpy as np
        import parmed
        s = parmed.load_file(str(prm), xyz=str(crd))
        fe = [a.idx for a in s.atoms if a.residue.name == "HEM" and a.name.upper() == "FE"]
        lig = [a.idx for a in s.atoms if a.residue.name == "LIG" and a.element > 1]
        xyz = s.coordinates
        return round(float(np.linalg.norm(xyz[lig] - xyz[fe[0]], axis=1).min()), 3)
    except Exception:  # noqa: BLE001
        return ""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--poselist", type=Path, required=True)
    ap.add_argument("--systems-root", type=Path, required=True)
    ap.add_argument("--work-dir", type=Path, required=True)
    ap.add_argument("--out-csv", type=Path, required=True)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--igb", type=int, default=2)
    ap.add_argument("--salt", type=float, default=0.15)
    args = ap.parse_args()

    poses = list(csv.DictReader(args.poselist.open(newline="")))
    poses = [p for i, p in enumerate(poses) if i % args.nshards == args.shard]
    args.work_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for n, p in enumerate(poses, 1):
        pid = p["pose_id"]
        sysd = args.systems_root / pid / pid
        # carry the poselist's own columns through (inputs excepted), whatever the campaign
        rec = {k: v for k, v in p.items() if k not in ("pdb", "smiles")}
        rec["status"] = ""
        wd = args.work_dir / pid
        try:
            prm, crd = sysd / "system.prmtop", sysd / "system.inpcrd"
            if not prm.is_file():
                raise RuntimeError("no system.prmtop (build failed)")
            wd.mkdir(parents=True, exist_ok=True)
            (wd / "min.in").write_text(MIN_IN)
            run(["sander", "-O", "-i", "min.in", "-p", str(prm), "-c", str(crd),
                 "-ref", str(crd), "-r", "min.rst", "-o", "min.out"], wd, "sander_min")
            crd = wd / "min.rst"
            # one-frame trajectory of the SOLVATED, MINIMISED system; MMPBSA strips it
            (wd / "mktraj.in").write_text(
                f"parm {prm}\ntrajin {crd}\ntrajout {wd}/frame.nc netcdf\ngo\nquit\n")
            run(["cpptraj", "-i", str(wd / "mktraj.in")], wd, "cpptraj")
            run(["ante-MMPBSA.py", "-p", str(prm), "-c", "complex.prmtop",
                 "-r", "receptor.prmtop", "-l", "ligand.prmtop",
                 "-s", ":WAT,Na+,Cl-", "-n", ":LIG", "--radii=mbondi2"], wd, "ante-MMPBSA")
            (wd / "mmgbsa.in").write_text(
                MMGBSA_IN.replace("igb=2", f"igb={args.igb}").replace("0.15", str(args.salt)))
            run(["MMPBSA.py", "-O", "-i", "mmgbsa.in", "-o", "FINAL_RESULTS_MMPBSA.dat",
                 "-sp", str(prm), "-cp", "complex.prmtop", "-rp", "receptor.prmtop",
                 "-lp", "ligand.prmtop", "-y", "frame.nc"], wd, "MMPBSA")
            rec.update(parse_delta(wd / "FINAL_RESULTS_MMPBSA.dat"))
            rec["min_lig_rmsd_A"] = lig_shift(prm, sysd / "system.inpcrd", wd / "min.rst", wd)
            rec["fe_lig_in_A"] = fe_lig_min(prm, sysd / "system.inpcrd")
            rec["fe_lig_out_A"] = fe_lig_min(prm, wd / "min.rst")
            rec["status"] = "ok"
        except Exception as exc:  # noqa: BLE001
            rec["status"] = f"FAIL {type(exc).__name__}: {str(exc)[:160]}"
        rows.append(rec)
        print(f"[{n}/{len(poses)}] {pid}: {rec.get('delta_total', rec['status'])}", flush=True)

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for r in rows:
        for k in r:
            if k not in fields:
                fields.append(k)
    with args.out_csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields); w.writeheader(); w.writerows(rows)
    ok = sum(1 for r in rows if r["status"] == "ok")
    print(f"\n{ok}/{len(rows)} modes scored -> {args.out_csv}")


if __name__ == "__main__":
    main()
