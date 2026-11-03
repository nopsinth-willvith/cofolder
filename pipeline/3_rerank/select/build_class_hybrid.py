"""Assemble a submission by ligand CLASS: type-II members from one zip, non-type-II from another.

The two classes want different things from a pose. Type-II ligands coordinate the heme iron, so a
pose source that ENFORCES that coordination (FRED with a Chelator constraint on Fe) is the natural
supplier. Non-type-II ligands must stay off the iron and are better served by a source that sampled
the pocket physically. Taking each class from the arm built for it is the whole idea; nothing is
re-docked or re-relaxed here.

Members are byte copies of members of zips that were each PB-audited from their own bytes; a
chosen member whose source audit row is not pb_pass True is refused (when that audit exists). The
finished zip is then re-audited from scratch by the caller -- that audit is the claim.

The class split is the MD build's own: data/cyp3a4_relax/nontypeii_base_ids.txt lists the nine
non-type-II ligands; everything else in the zip is type-II.
"""
from __future__ import annotations

import argparse
import os
import csv
import hashlib
import sys
import zipfile
from pathlib import Path

ROOT = Path(os.environ.get("COFOLD_ROOT", "/project/bme-20255165-cofolding-physics"))
SUB = ROOT / "data/cyp3a4_relax/submission"
AUD = ROOT / "analysis/cyp3a4_gamd/results"


def members(name: str, sub_dir: Path) -> dict[str, bytes]:
    with zipfile.ZipFile(sub_dir / f"{name}.zip") as zf:
        return {Path(n).stem: zf.read(n) for n in zf.namelist() if n.endswith(".pdb")}


def audit(name: str, sub_dir: Path, audit_dir: Path) -> dict[str, str] | None:
    """Find the source zip's PoseBusters audit, under either campaign's naming.

    The MD/GaMD stages write analysis/cyp3a4_gamd/results/zipaudit_<name>.csv; the Mol-1 and OB0
    stages write <name>.zipaudit.csv NEXT TO the zip (slurm/3_rerank/select.sbatch). Both are tried so a
    caller does not have to know which convention its parents used -- and returning None (meaning
    "no audit found") silently downgrades the pb_pass refusal to the finished zip's own re-audit,
    so guessing one pattern would weaken the gate without saying so.
    """
    for p in (audit_dir / f"zipaudit_{name}.csv", sub_dir / f"{name}.zipaudit.csv"):
        if p.is_file():
            with p.open(newline="") as fh:
                return {r["base_id"]: r.get("pb_pass", "") for r in csv.DictReader(fh)}
    return None


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--typeii-from", required=True, help="submission zip basename")
    ap.add_argument("--nontypeii-from", required=True, help="submission zip basename")
    ap.add_argument("--name", required=True)
    ap.add_argument("--nontypeii-list", type=Path,
                    default=ROOT / "data/cyp3a4_relax/nontypeii_base_ids.txt")
    # Both source and destination directories: this script is not MD-specific -- the same member
    # copy serves any pair of parent zips, e.g. Mol-1's combined-pool conf and MM-GBSA submissions
    # in data/cyp3a4_mol1/submission.
    ap.add_argument("--sub-dir", type=Path, default=SUB,
                    help="directory holding the parent zips; the output zip is written here too")
    ap.add_argument("--audit-dir", type=Path, default=AUD,
                    help="directory holding zipaudit_<name>.csv (the <name>.zipaudit.csv form "
                         "next to the zip is also tried)")
    args = ap.parse_args(argv)
    SUB_DIR, AUD_DIR = args.sub_dir, args.audit_dir

    nont = {x.strip() for x in args.nontypeii_list.read_text().split() if x.strip()}
    src = {"type-II": args.typeii_from, "non-type-II": args.nontypeii_from}
    mem = {k: members(v, SUB_DIR) for k, v in src.items()}
    aud = {k: audit(v, SUB_DIR, AUD_DIR) for k, v in src.items()}
    ids = sorted(set(mem["type-II"]) | set(mem["non-type-II"]))
    if len(ids) != 20:
        sys.exit(f"expected 20 ligands across the sources, found {len(ids)}")

    rows, out = [], {}
    for pid in ids:
        cls = "non-type-II" if pid in nont else "type-II"
        b = mem[cls].get(pid)
        if b is None:
            sys.exit(f"{pid} ({cls}) is missing from {src[cls]}")
        a = aud[cls]
        if a is not None and str(a.get(pid, "")).lower() != "true":
            sys.exit(f"{pid} ({cls}): source {src[cls]} audit says pb_pass={a.get(pid)!r}")
        out[pid] = b
        rows.append({"pose_id": pid, "class": cls, "source": src[cls],
                     "source_audited": a is not None,
                     "sha256_16": hashlib.sha256(b).hexdigest()[:16]})

    zp = SUB_DIR / f"{args.name}.zip"
    with zipfile.ZipFile(zp, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for pid in ids:
            zf.writestr(f"{pid}.pdb", out[pid])
    prov = SUB_DIR / f"{args.name}.provenance.csv"
    with prov.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    n2 = sum(r["class"] == "type-II" for r in rows)
    print(f"{args.name}: {n2} type-II from {args.typeii_from}, "
          f"{len(rows) - n2} non-type-II from {args.nontypeii_from}")
    print(f"  {zp}  sha256 {hashlib.sha256(zp.read_bytes()).hexdigest()[:16]}")
    for k, v in aud.items():
        if v is None:
            print(f"  note: {src[k]} has no zipaudit yet -- its members are gated only by the "
                  "re-audit of the finished zip")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
