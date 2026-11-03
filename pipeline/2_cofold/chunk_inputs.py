"""Pack an arm's remaining ligands into a few multi-job inputs instead of one input per ligand.

One ligand per array task was the wrong unit on a contended cluster. Each task is its own queue
position AND reloads the 2.39 GB checkpoint (~1 min of a ~15 min task), so 20 ligands cost 20 waits
and 20 loads. Mol-1's input JSON is a LIST of jobs (docs/infer_json_format.md: "Even if you are
modeling only one set of sequences, the top-level structure should still be a list"), and the
runner iterates it with the model loaded once -- so N ligands in one file cost ONE wait and ONE
load. With five arms queued behind other users, that is the difference that matters, not GPU-seconds.

Chunking rather than one-file-per-arm keeps three properties worth having: some parallelism when
slots do free up, a bounded blast radius (inference.py raises after finishing if ANY sample failed,
so a chunk is the unit that has to be re-run), and a wall time that stays inside the 8 h limit.

Only ligands that are NOT already complete are packed, so this is safe to re-run after a partial
arm: it reads what is on disk rather than assuming the arm started empty.

Env: any python3 (stdlib). CPU, instant.
"""

from __future__ import annotations

import argparse
import os
import json
import sys
from pathlib import Path

R = Path(os.environ.get("COFOLD_ROOT", "/vast.mnt/project/bme-20255165-cofolding-physics"))


def complete(gen_dir: Path, name: str, expect: int) -> bool:
    d = gen_dir / name
    return d.is_dir() and len(list(d.glob("seed_*/predictions/*_sample_*.cif"))) >= expect


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, allow_abbrev=False)
    ap.add_argument("--inputs-dir", type=Path, required=True,
                    help="a per-ligand input dir written by make_inputs.py")
    ap.add_argument("--pool", required=True, help="arm label; its gen dir is checked for progress")
    ap.add_argument("--chunk-size", type=int, default=5)
    ap.add_argument("--expect", type=int, default=50, help="poses per ligand when complete")
    ap.add_argument("--gen-root", type=Path, default=R / "data/cyp3a4_mol1/gen")
    a = ap.parse_args()

    per_ligand = sorted(p for p in a.inputs_dir.glob("*.json") if not p.name.startswith("chunk_"))
    if not per_ligand:
        sys.exit(f"FATAL: no per-ligand inputs in {a.inputs_dir}")
    gen = a.gen_root / a.pool

    todo, done = [], []
    for p in per_ligand:
        (done if complete(gen, p.stem, a.expect) else todo).append(p)
    print(f"{a.pool}: {len(done)} ligand(s) already complete, {len(todo)} to run")
    if not todo:
        print("  nothing to do")
        (a.inputs_dir / f"chunks_index_{a.pool}.txt").write_text("")
        return 0

    # Per-POOL, not per-inputs-dir: the MSA-depth arms all share one per-ligand inputs dir, so a
    # single "chunks/" would have each arm overwrite the previous arm's index.
    out_dir = a.inputs_dir / f"chunks_{a.pool}"
    out_dir.mkdir(exist_ok=True)
    for stale in out_dir.glob("chunk_*.json"):
        stale.unlink()                      # a stale chunk would silently refold finished ligands

    written = []
    for i in range(0, len(todo), a.chunk_size):
        group = todo[i:i + a.chunk_size]
        jobs = []
        for p in group:
            loaded = json.loads(p.read_text())
            if not isinstance(loaded, list) or len(loaded) != 1:
                sys.exit(f"FATAL: {p} is not a single-job list")
            jobs.append(loaded[0])
        c = out_dir / f"chunk_{i // a.chunk_size}.json"
        c.write_text(json.dumps(jobs, indent=2))
        written.append(c)
        print(f"  {c.name}: {len(jobs)} ligands ({', '.join(j['name'] for j in jobs)})")

    index = a.inputs_dir / f"chunks_index_{a.pool}.txt"
    # absolute: the generation sbatch cds into the Mol-1 repo before calling the runner
    index.write_text("".join(f"{c.resolve()}\n" for c in written))
    print(f"{len(written)} chunk(s) -> {index}")
    print(f"  {len(written)} array task(s); slurm/run.sh cofold reads the count from {index.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
