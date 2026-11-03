# Vendored from $COFOLD_ROOT/analysis/cyp3a4_relax/eint_prep_shard.py -- unchanged apart from the COFOLD_ROOT default and the sibling import of prep_cofold_system.
"""Stage E2 -- prep one shard of the E_int manifest into per-pose npz files.

Drives `analysis/cyp_challenge/prep_cofold_system.prep_one` in-process so the imports and the
490 MB CCD scan are paid once per shard rather than once per pose (the HEM template is cached to
json after the first call, but the rdkit/openmm import cost is not).

PER-POSE FAILURE IS RECORDED, NEVER FATAL. The PXR run lost ~15% of its pool to RDKit bond-order
perception (`ValueError: No matching found` was 2642 of its failures) and the authors flagged the
resulting coverage gap as a confound -- candidates that failed protonation were never scored, so
a model could win by walkover. Here every pose writes a row either way, with the exception text
when it fails, so the coverage of the energy table is visible rather than inferred.

The spin written into each npz comes from the MANIFEST, not from prep_cofold_system's default.
prep_one has no coordination detection and would stamp hem_spin_mult=6 on every pose; 62.3% of
this manifest coordinates the heme iron through N and belongs at multiplicity 2, a 1.24 eV
difference measured on this very system.

Env: openmm83 (the only env with gemmi + rdkit + openmm together). CPU.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("COFOLD_ROOT", "/project/bme-20255165-cofolding-physics"))
sys.path.insert(0, str(Path(__file__).resolve().parent))  # prep_cofold_system.py, vendored alongside


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", type=Path,
                    default=ROOT / "data/cyp3a4_relax/eint_manifest.csv")
    ap.add_argument("--out-dir", type=Path,
                    default=Path("/scratch-shared/20255165/cyp3a4_relax/eint_prep"))
    ap.add_argument("--status-csv", type=Path, required=True)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--overwrite", action="store_true")
    a = ap.parse_args(argv)

    from prep_cofold_system import prep_one   # noqa: E402  (needs the sys.path insert above)

    # Was hardcoded "propionate". It is env-driven because the heme state has to be able to
    # TRACK the pH: PREP_PH only reaches PDBFixer (the protein), so setting PREP_PH=4.8 while
    # this stayed pinned left the heme in its pH-7.4 form -- both propionates ionised -- inside
    # a system protonated for pH 4.8, with nothing raising. Propionic acid pKa is 4.87, so at
    # pH 4.8 the two states bracket the truth rather than one being correct.
    # Default unchanged, so every existing campaign reproduces byte for byte.
    HEM_PROTONATION = os.environ.get("HEM_PROTONATION", "propionate")
    if HEM_PROTONATION not in ("propionate", "acid"):
        sys.exit(f"FATAL: HEM_PROTONATION={HEM_PROTONATION!r}, expected propionate|acid")
    print(f"[prep] PREP_PH={os.environ.get('PREP_PH', '7.4')} "
          f"HEM_PROTONATION={HEM_PROTONATION}", flush=True)

    M = pd.read_csv(a.manifest, low_memory=False)
    M = M.iloc[a.shard::a.nshards]
    a.out_dir.mkdir(parents=True, exist_ok=True)
    print(f"shard {a.shard}/{a.nshards}: {len(M)} poses -> {a.out_dir}", flush=True)

    rows, t0 = [], time.time()
    for n, r in enumerate(M.itertuples(), 1):
        out = a.out_dir / f"{r.tag}.npz"
        rec = {"tag": r.tag, "arm": r.arm, "query_id": r.query_id, "seed": int(r.seed),
               "sample": int(getattr(r, "sample")), "base_id": r.base_id,
               "npz": str(out), "spin_complex": int(r.spin_complex),
               "fe_coordinating": bool(r.fe_coordinating), "status": "", "n_atoms": -1}
        if out.is_file() and not a.overwrite:
            rec["status"] = "skip_exists"
            rows.append(rec)
            continue
        try:
            prep_one(str(r.cif), str(r.smiles), str(out), lig_res=str(r.lig_res),
                     hem_h="ccd", hem_protonation=HEM_PROTONATION,
                     hem_spin=int(r.spin_complex),     # per-pose, from the manifest
                     cys_thiolate="auto", fe_sg_max=4.0, label_fe=True)
            Z = np.load(out, allow_pickle=False)
            rec["n_atoms"] = int(Z["n_atoms"])
            # the heme MUST be on the receptor side; a heme swept onto the ligand turns E_int
            # into "drug+heme vs apo protein", which is large, plausible and silently wrong
            if bool((Z["is_hem"] & Z["is_ligand"]).any()):
                raise ValueError("HEM atoms flagged is_ligand -- receptor/substrate split inverted")
            if bool((Z["is_fe"] & Z["is_ligand"]).any()):
                raise ValueError("Fe flagged is_ligand -- receptor/substrate split inverted")
            if int(Z["hem_spin_mult"]) != int(r.spin_complex):
                raise ValueError(f"npz spin {int(Z['hem_spin_mult'])} != manifest "
                                 f"{int(r.spin_complex)}")
            rec["status"] = "ok"
        except Exception as e:  # noqa: BLE001 -- recorded per pose; coverage must be visible
            rec["status"] = f"{type(e).__name__}: {e}"[:200]
            print(f"  FAIL {r.tag}: {rec['status']}", flush=True)
            traceback.print_exc(limit=2)
        rows.append(rec)
        if n % 25 == 0:
            ok = sum(1 for x in rows if x["status"] in ("ok", "skip_exists"))
            print(f"  {n}/{len(M)}  ok={ok}  {time.time()-t0:.0f}s", flush=True)

    D = pd.DataFrame(rows)
    a.status_csv.parent.mkdir(parents=True, exist_ok=True)
    D.to_csv(a.status_csv, index=False)
    good = D.status.isin(["ok", "skip_exists"]).sum()
    print(f"\nshard {a.shard}: {good}/{len(D)} prepped  ({time.time()-t0:.0f}s)")
    print(D.status.value_counts().head(8).to_string())
    if good:
        print(f"atoms: median {D.loc[D.n_atoms > 0, 'n_atoms'].median():.0f}  "
              f"max {D.n_atoms.max()}")
    print(f"wrote {a.status_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
