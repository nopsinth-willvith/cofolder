# Vendored from $COFOLD_ROOT/analysis/cyp_challenge/prep_cofold_system.py -- unchanged apart from the COFOLD_ROOT default and the sibling import of prep_system.
"""Stage 50 -- co-folded CYP CIF -> MLFF-ready protonated npz (protein + HEM + substrate).

This is the CYP fork of ``analysis/affinity_bench/prep_cofold_system.py``.  Everything
about hydrogens, the 12 A pocket cut and the ``protenix_atom_index`` join is kept
IDENTICAL to that reference so the descriptors stay comparable with the affinity_bench
campaign.  The one thing that is genuinely new is the heme, and the heme is where every
silent-correctness trap lives.

WHAT A CYP CO-FOLD CIF CONTAINS (verified against the CIF, never assumed -- see
``split_groups``).  ``build_cofold_inputs.py`` emits three entities:

    chain A  proteinChain  the isoform receptor          -> standard amino-acid res_names
    chain B  ligand SMILES the challenge compound        -> res_name "l01"   (renamed)
    chain C  ligand CCD_HEM the heme cofactor            -> res_name "HEM"   (NOT renamed)

Protenix renames only SMILES ligands (l01..l99) and explicitly skips that rename for
``CCD_`` entities, so substrate and heme are separable by res_name.  ``split_groups``
still CHECKS this instead of trusting it: it auto-detects the substrate residue when the
expected ``l01`` is absent and refuses to guess when more than one candidate exists.

--------------------------------------------------------------------------------------
TRAP 1 -- HEM IS A *LIGAND ENTITY* AT CO-FOLD TIME AND A *RECEPTOR ATOM* HERE.
--------------------------------------------------------------------------------------
``model/mlff/backends/base.py:73-80`` builds the three single points by slicing on the
``is_ligand`` mask::

    PROTEIN = ~is_ligand      LIGAND = is_ligand      E_int = E_cx - E_prot - E_lig

so a heme counted as ligand silently redefines E_int as "drug+heme vs APO protein"
instead of "drug vs HOLO protein".  Nothing crashes; every energy feature is simply
wrong, by ~the entire drug-heme interaction.  Therefore:

    is_ligand := (res_name == <substrate res>)          <- NEVER "the non-protein atoms"

and ``prep_one`` ASSERTS that the ligand group contains no Fe and no HEM residue, and
prints per-group atom counts on every system.

--------------------------------------------------------------------------------------
TRAP 2 -- SPIN.  ``base.py:63-64`` derives multiplicity from electron-count parity and can
only ever return 1 or 2.  A substrate-bound ferric P450 heme is HIGH SPIN, S = 5/2,
multiplicity 6.  We cannot fix base.py from here, so the npz carries explicit
``hem_spin_mult`` / ``spin_hint_{complex,protein,ligand}`` fields for the Orb path to use
instead of the parity rule.  (POLAR takes charge only and is unaffected.)
Caveat recorded on purpose: an azole/pyridine inhibitor that coordinates the iron drives
it back to LOW spin (mult 2).  ``--hem_spin`` overrides.

--------------------------------------------------------------------------------------
TRAP 3 -- CHARGE IS A DECISION, NOT A DEFAULT.  POLAR conditions on the total charge
(``extract_polar_latents.py:79`` feeds ``ligand_charge + pocket_charge`` into
``atoms.info["charge"]``), so a wrong total silently changes every descriptor.  The heme
contribution is derived, printed and stored, never defaulted:

    Fe(III)  +3
    porphyrin dianion  -2      (the CCD HEM block has no H on NA/NB/NC/ND)
    2 x propionate     -2      (pH 7.4; we DROP the CCD acid protons H2A/H2D)
    ------------------------------
    HEM net             -1     = HEM_CHARGE["propionate"], the default

    --hem_protonation acid  keeps H2A/H2D and gives +1 instead.  The protonation state and
    the charge are read from ONE table so they can never disagree.

    The proximal cysteine is a THIOLATE when it ligates the iron, worth another -1.
    PDBFixer protonates it as neutral Cys-SH at pH 7.4 and ff14SB then reports it neutral,
    so ``--cys_thiolate auto`` (default) finds the Cys SG within ``--fe_sg_max`` of the Fe,
    DELETES the HG that PDBFixer just put on it (it would otherwise sit inside the iron)
    and books the -1 explicitly.

    DOWNSTREAM CONTRACT -- THE THREE CHARGES ARE SEPARABLE AND MUST BE SUMMED:

        total = ligand_charge + pocket_charge + hem_charge

    ``pocket_charge`` EXCLUDES the heme and INCLUDES the thiolate correction, which is the
    contract ``extract_polar_cyp.py`` and ``orb_energies_cyp.py`` (stages 60/65) read:
    both do ``q_lig + q_pock + q_hem`` and both take ``hem_charge`` straight from this npz.
    Folding the heme into ``pocket_charge`` would DOUBLE-COUNT it there.

    Corollary, and it matters: DO NOT run the affinity_bench
    ``analysis/affinity_bench/extract_polar_batch.py`` on these npz.  That script predates
    the heme and computes ``ligand_charge + pocket_charge``, so it would silently drop the
    heme's -1 and hand POLAR the wrong conditioning charge.  Use the CYP stage-60 script.

--------------------------------------------------------------------------------------
TRAP 4 -- THE Fe ROW, NOT THE POOLED HEME.  The heme is bit-identical across every complex
of an isoform, so pooling all 43 heme atoms produces near-constant descriptor dimensions on
a problem that already has only ~3 label rows per dimension.  ``labelled`` therefore covers
the substrate atoms PLUS the single Fe atom (``--no_label_fe`` turns the Fe off).
``label_group`` and ``fe_protenix_index`` are stored so the pooling stage can keep the Fe
row OUT of the ligand mean/max and give it its own block.  IF YOU POOL ``labelled`` AS ONE
GROUP YOU WILL SILENTLY MIX THE Fe ROW INTO THE LIGAND MEAN -- split on ``label_group``.

--------------------------------------------------------------------------------------
HYDROGENS.  Three different sources, deliberately:
  * amino acids -- PDBFixer @ pH 7.4 + ff14SB H-only relax with heavy atoms frozen
    (``analysis/energy_labels/prep_system.py:_minimize_frozen_heavy``), i.e. byte-for-byte
    the affinity_bench recipe.
  * substrate   -- RDKit ``AddHs(addCoords=True)`` + H-only MMFF/UFF relax
    (``prep_system.py:_relax_hydrogens``), again the affinity_bench recipe.
  * HEM         -- NEITHER of those works: amber14-all.xml has no HEM template, so feeding
    the heme to PDBFixer makes ``createSystem`` raise, and the reference code CATCHES that
    exception and falls back with ``pocket_charge = 0`` -- a silent charge bug.  The heme is
    therefore excluded from the PDBFixer block and protonated from its own CCD template by
    rigid local-frame superposition (``place_hem_hydrogens``): for each template hydrogen,
    Kabsch-fit its parent heavy atom and that parent's heavy neighbours from the CCD ideal
    geometry onto the pose, then carry the hydrogen through the same transform.  Bond
    lengths and angles are exactly ideal by construction; only 2nd-shell torsions (methyl
    rotamers, propionate tails) are approximate.

POCKET_CUT is THE MEMORY LEVER, not the heme.  POLAR's long-range module is superlinear in
atom count and a 707-atom system measured 22.6 GB; the heme's +43 heavy atoms are noise
next to that.  Default 12 A = MACE-POLAR-1-M's receptive field (r_max 6.0 x 2 layers).
The cut is measured from the UNION of the substrate and the heme heavy atoms, because the
whole heme is kept unconditionally: measuring from the substrate alone leaves the heme
periphery (and hence the headline Fe row, trap 4) with a TRUNCATED environment, which is
exactly what the receptive-field justification for 12 A is supposed to rule out.  Every
system prints the extra atoms/residues that the heme reference adds, its total atom count,
and warns above WARN_ATOMS (default 700, the measured 22.6 GB point).

INPUTS THIS SCRIPT REQUIRES (it fails loudly and by name if any is absent -- it never
emits an empty or partial npz):
  data/cyp_challenge/cofold_inputs/<tag>/manifest.csv   stage 30; pair_id -> pose_glob,
        smiles, isoform, has_hem_requested, fe_cys_bond_requested
  /scratch-shared/20255165/cyp_challenge/predictions/<tag>/<pair_id>/seed_101/predictions/
        <pair_id>_sample_*.cif                          stage 40 -- NOT YET RUN as of writing
  /vast.mnt/home/20255165/common/components.cif         the CCD, for the HEM H template

Env: miniforge3/envs/openmm83 -- the only env on this cluster carrying gemmi (the CIF
parser) AND rdkit AND openmm/pdbfixer together.  The `cyp` env has rdkit but no gemmi and
no openmm, so it cannot parse a CIF or protonate anything.

Usage:
  # one system
  prep_cofold_system.py --single --cif <cif> --smiles <smi> --out <npz>
  # array shard (what slurm/cyp_challenge/50_prep_array.sbatch does)
  prep_cofold_system.py --manifest <manifest.csv> --out_dir <dir> --shard i --nshards n
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
import traceback
from pathlib import Path

import numpy as np

ROOT = Path(os.environ.get("COFOLD_ROOT", "/project/bme-20255165-cofolding-physics"))
DATA = ROOT / "data/cyp_challenge"
SCRATCH = Path("/scratch-shared/20255165/cyp_challenge")

# --- knobs (env-overridable so the sbatch can set them without editing code) ------------
POCKET_CUT = float(os.environ.get("POCKET_CUT", "12.0"))
MAX_ATOMS = int(os.environ.get("MAX_ATOMS", "8000"))     # hard refusal
WARN_ATOMS = int(os.environ.get("WARN_ATOMS", "700"))    # 707 atoms measured 22.6 GB in POLAR
PH = float(os.environ.get("PREP_PH", "7.4"))
CCD_COMPONENTS = Path(os.environ.get(
    "CCD_COMPONENTS", "/vast.mnt/home/20255165/common/components.cif"))

# POCKET_REF / POCKET_FORCE_RES -- additive knobs for the 5 A substrate-centred cluster arm
# (slurm/cyp3a4_fgfree_eint).  BOTH DEFAULT TO THE EXISTING BEHAVIOUR: POCKET_REF
# "substrate_hem" reproduces the hardcoded union below line for line, and an empty
# POCKET_FORCE_RES adds nothing, so a default-knob re-prep must reproduce an existing npz.
#   POCKET_REF="substrate_hem" (default)  cut measured from the UNION of the substrate and
#     the heme heavy atoms.  That union is what the receptive-field argument for 12 A
#     (r_max 6.0 x 2 layers) requires: every retained atom, the Fe included, then has a
#     complete environment.
#   POCKET_REF="substrate"  cut measured from the substrate heavy atoms ALONE.  This is the
#     5 A cluster's reference.  At 5 A a heme-inclusive cut still drags in the whole distal
#     face and defeats the purpose, which is a system POLAR's superlinear long-range module
#     can hold (707 atoms measured 22.6 GB).  prep_one still keeps the heme WHOLE and
#     unconditionally, so in this mode the heme PERIPHERY has a truncated environment: a
#     deliberate, recorded trade for a tractable cluster, not an oversight.
#   POCKET_FORCE_RES  comma-separated auth_seq numbers kept regardless of distance, as WHOLE
#     residues through the SAME PDBFixer/ff14SB path as the distance-selected ones, so
#     pocket_charge stays correct.  "442,309,308" = Cys442 axial thiolate + the Thr309/Glu308
#     acid-alcohol pair.  At 5 A from the substrate the thiolate can fall out of the crop
#     entirely, which leaves the Fe coordination sphere incomplete while still producing a
#     perfectly plausible energy -- so a requested residue that is ABSENT is fatal, below.
POCKET_REF = os.environ.get("POCKET_REF", "substrate_hem")
POCKET_REF_CHOICES = ("substrate_hem", "substrate")
if POCKET_REF not in POCKET_REF_CHOICES:
    raise SystemExit(f"POCKET_REF={POCKET_REF!r} is not one of {POCKET_REF_CHOICES}")
POCKET_FORCE_RES = os.environ.get("POCKET_FORCE_RES", "")


def _parse_force_res(spec: str) -> list[int]:
    """"442,309,308" -> [442, 309, 308].  A non-integer token is fatal, never skipped: a
    typo that silently dropped Cys442 is exactly the failure this knob exists to prevent."""
    out = []
    for tok in str(spec).replace(";", ",").split(","):
        tok = tok.strip()
        if not tok:
            continue
        try:
            out.append(int(tok))
        except ValueError:
            raise SystemExit(f"POCKET_FORCE_RES={spec!r}: {tok!r} is not an integer "
                             f"auth_seq residue number") from None
    return out


POCKET_FORCE_RES_LIST = _parse_force_res(POCKET_FORCE_RES)

# PREP_SEED -- OPT-IN reproducibility for the pocket hydrogens.  DEFAULT UNSET = the historical
# unseeded behaviour, so the cyp_challenge campaign and the existing 12 A npz are unaffected.
#
# WHY IT EXISTS (measured, slurm/cyp3a4_fgfree_eint/01b-01d).  Re-prepping the SAME pose twice
# with the SAME code gives bitwise-identical heavy atoms, masks and charges but ~1441 DIFFERENT
# pocket hydrogen positions (mean ~0.6 A, max ~2.5 A), because
# openmm/app/modeller.py:1016-1023 -- inside Modeller.addHydrogens, which
# PDBFixer.addMissingHydrogens calls -- places every added hydrogen at a RANDOMLY perturbed
# position and relies on the following minimisation to fix it up ("The hydrogens were added at
# random positions", modeller.py:1038).  Python's `random` is seeded from OS entropy per process,
# so each prep run is a different realisation; PYTHONHASHSEED does NOT control it.  Seeding
# `random` made two preps bitwise identical on all 59 npz keys.
#
# The cost of leaving it unseeded, also measured: two independent 5 A realisations of the same
# pose move MACE-POLAR E_int by a median 0.027 eV (max 0.056), against a 0.481 eV across-pose
# spread.  Small, but it makes an energy table unreproducible, so any NEW campaign should set it.
PREP_SEED = os.environ.get("PREP_SEED", "")
PREP_SEED_VAL: int | None = None
if PREP_SEED.strip() != "":
    try:
        PREP_SEED_VAL = int(PREP_SEED)
    except ValueError:
        raise SystemExit(f"PREP_SEED={PREP_SEED!r} is not an integer") from None

# --- heme constants (CCD data_HEM, verified by scanning components.cif) -----------------
HEM_RES = "HEM"
HEM_FE_ATOM = "FE"
HEM_N_ATOMS = ("NA", "NB", "NC", "ND")
HEM_N_HEAVY = 43                      # C34 Fe N4 O4
HEM_ACID_H = ("H2A", "H2D")           # the two propionic-acid OH protons
# protonation state -> (H atom_ids to DROP, net formal charge). See TRAP 3.
HEM_PROTONATION = {
    "propionate": (HEM_ACID_H, -1),   # pH 7.4 default: both propionates ionised
    "acid": ((), +1),                 # neutral propionic acids
    # heavy-atom-only heme (--hem_h none): charge is None ON PURPOSE. -1 is derived for
    # C34H30FeN4O4(-); a bare C34FeN4O4 skeleton is 32 H lighter and no tabulated charge
    # applies, so prep_one REFUSES this mode unless --hem_charge is given explicitly.
    "none": (None, None),
}
HEM_SPIN_MULT = 6                     # substrate-bound ferric P450: S = 5/2 (TRAP 2)

CYS_SG_ATOM = "SG"
LIG_RES_RE = re.compile(r"^l\d{2}$")          # Protenix SMILES-ligand rename
SAMPLE_RE = re.compile(r"_sample_(\d+)\.cif$")

STD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU",
    "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "MSE", "SEC", "PYL",
}
TMP_CHAINS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"

GROUP_NAMES = ("protein", "hem", "substrate")
G_PROTEIN, G_HEM, G_SUBSTRATE = 0, 1, 2
CODE_VERSION = "cyp-prep-1"

_HEM_TEMPLATE: dict | None = None      # process-level cache; the CCD scan is ~490 MB


# --------------------------------------------------------------------------- CIF reading
def _as_int(v, default: int = -1) -> int:
    """mmCIF numeric field that may legitimately be the null token '.' or '?'."""
    try:
        return int(str(v).strip())
    except (TypeError, ValueError):
        return default


def read_atom_site(cif_path):
    """Read ``_atom_site`` in FILE ORDER -> dict of arrays + row index.

    Verbatim from ``analysis/affinity_bench/prep_cofold_system.py`` apart from the
    tolerant auth_seq parse.  It is copied rather than imported because that file has the
    SAME MODULE NAME as this one, so ``import prep_cofold_system`` would resolve to
    whichever directory happened to be first on sys.path.  (Proper fix: lift it into
    model/utils/ when a third copy appears.)

    Row order IS the Protenix featuriser's atom order -- ``save_structure_cif`` only
    overwrites ``.coord`` -- so ``row`` is the join key back into the trunk features and
    no 5-key (chain, seq, comp, atom, alt) join is ever needed.
    """
    import gemmi

    block = gemmi.cif.read(str(cif_path)).sole_block()
    cols = ["group_PDB", "label_atom_id", "label_comp_id", "label_asym_id",
            "auth_seq_id", "type_symbol", "Cartn_x", "Cartn_y", "Cartn_z"]
    tab = block.find("_atom_site.", cols)
    n = len(tab)
    if n == 0:
        raise ValueError(f"no _atom_site rows in {cif_path}")
    out = {
        "group": np.array([tab[i][0] for i in range(n)], dtype=object),
        "atom_name": np.array([tab[i][1].strip('"') for i in range(n)], dtype=object),
        "res_name": np.array([tab[i][2] for i in range(n)], dtype=object),
        "asym": np.array([tab[i][3] for i in range(n)], dtype=object),
        "auth_seq": np.array([_as_int(tab[i][4]) for i in range(n)], dtype=np.int32),
        "element": np.array([tab[i][5].capitalize() for i in range(n)], dtype=object),
        "coord": np.array([[float(tab[i][6]), float(tab[i][7]), float(tab[i][8])]
                           for i in range(n)], dtype=np.float64),
    }
    out["row"] = np.arange(n, dtype=np.int32)
    return out


def sorted_poses(paths):
    """CIF paths ordered by their numeric ``_sample_<k>`` index, not lexicographically.

    ``sorted()`` puts sample_10 before sample_2.  With N_sample=5 that never bites, which
    is exactly why it would go unnoticed if N_sample were ever raised.
    """
    def key(p):
        m = SAMPLE_RE.search(str(p))
        return (0, int(m.group(1))) if m else (1, str(p))
    return sorted(paths, key=key)


def kept_sample(cif_path) -> int | None:
    """The sample index the PRUNER kept, from ``<seeddir>/kept.json``, or None.

    ``analysis/affinity_bench/prune_samples.py:70-72`` writes ``kept_sample`` = the sample
    with the HIGHEST ``ranking_score``, i.e. Protenix's own verdict on which pose is best.
    Taking the lowest ``_sample_<k>`` index instead is only equivalent because PRUNE=1
    normally leaves exactly one CIF behind -- with PRUNE=0, or if pruning failed after
    inference succeeded, it silently picks an arbitrary sample.
    """
    kj = Path(cif_path).parent / "kept.json"
    if not kj.exists():
        return None
    try:
        return int(json.loads(kj.read_text())["kept_sample"])
    except (OSError, ValueError, KeyError, TypeError) as exc:
        print(f"  WARNING: {kj} is unreadable ({type(exc).__name__}: {exc}); "
              f"falling back to the lowest sample index")
        return None


def split_groups(at, lig_res: str | None = None, verbose: bool = True) -> dict:
    """protein / HEM / substrate masks, VERIFIED against the CIF rather than assumed.

    Returns ``{"protein", "hem", "substrate"}`` boolean masks plus the resolved
    ``lig_res``, ``has_hem`` and the Fe row.  Prints per-group atom counts: a silent
    mis-split here poisons every energy feature downstream (TRAP 1), so it is made loud.
    """
    res = at["res_name"]
    is_h = at["element"] == "H"          # Protenix predicts heavy atoms only; belt and braces
    protein = np.array([r in STD_AA for r in res]) & ~is_h
    hem = (res == HEM_RES) & ~is_h

    other_res = sorted({str(r) for r in res[~protein & ~hem & ~is_h]})
    if lig_res is None:
        cands = [r for r in other_res if LIG_RES_RE.match(r)]
        if len(cands) == 1:
            lig_res = cands[0]
        elif len(cands) == 0:
            raise ValueError(
                f"no substrate residue found: non-protein non-HEM res_names are "
                f"{other_res or '[]'} and none matches the Protenix SMILES-ligand rename "
                f"l01..l99. Was this pose built with --hem_mode none and no ligand, or is "
                f"the substrate a CCD entity?")
        else:
            raise ValueError(
                f"ambiguous substrate residue: {cands} all match l01..l99. This builder "
                f"emits exactly ONE SMILES ligand, so refusing to guess.")
    substrate = (res == lig_res) & ~is_h
    if not substrate.any():
        raise ValueError(f"no atoms with res_name {lig_res!r}; present: {other_res}")

    leftover = ~protein & ~hem & ~substrate & ~is_h
    if leftover.any():
        drop = sorted({str(r) for r in res[leftover]})
        print(f"  WARNING: {int(leftover.sum())} atoms in unclassified residues {drop} "
              f"are DROPPED (not protein, not HEM, not the substrate)")

    fe = hem & (at["element"] == "Fe")
    has_hem = bool(hem.any())
    if verbose:
        print(f"  groups: protein {int(protein.sum())} | HEM {int(hem.sum())} "
              f"| substrate[{lig_res}] {int(substrate.sum())} | dropped {int(leftover.sum())}")
    if has_hem and int(hem.sum()) != HEM_N_HEAVY:
        print(f"  WARNING: HEM has {int(hem.sum())} heavy atoms, expected {HEM_N_HEAVY}")
    if has_hem and not fe.any():
        raise ValueError("HEM residue present but it has no Fe atom -- the pose is corrupt")

    return {"protein": protein, "hem": hem, "substrate": substrate, "leftover": leftover,
            "lig_res": lig_res, "has_hem": has_hem, "fe": fe}


# ----------------------------------------------------------------------------- substrate
def map_ligand(at, lig_res: str, smiles: str):
    """Rebuild substrate chemistry from SMILES and attach the predicted coordinates.

    Protenix names ligand atoms ``f"{element}{k}"`` following the RDKit atom order of
    ``MolFromSmiles``, so heavy atoms map by that same counter -- no substructure match.
    Matching is done on UPPERCASE because the CIF writes two-letter elements as CL1/BR1
    while RDKit's GetSymbol() gives "Cl"/"Br"; without folding the case every halogen-
    bearing ligand fails to map.

    Returns ``(mol_heavy_with_conformer, cif_row_of_mol_idx)``.
    """
    from rdkit import Chem

    sel = at["res_name"] == lig_res
    names = [str(nm).upper() for nm in at["atom_name"][sel]]
    coords = at["coord"][sel]
    rows = at["row"][sel]

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles}")
    counter, name_of_idx = {}, {}
    for a in mol.GetAtoms():
        sym = a.GetSymbol()
        counter[sym] = counter.get(sym, 0) + 1
        name_of_idx[a.GetIdx()] = f"{sym}{counter[sym]}".upper()
    idx_of_name = {v: k for k, v in name_of_idx.items()}

    missing = [nm for nm in names if nm not in idx_of_name]
    if missing or mol.GetNumAtoms() != len(names):
        raise ValueError(
            f"ligand/SMILES mismatch: cif has {len(names)} atoms, SMILES has "
            f"{mol.GetNumAtoms()}; unmatched names {missing[:5]}")

    conf = Chem.Conformer(mol.GetNumAtoms())
    row_of_idx = np.full(mol.GetNumAtoms(), -1, dtype=np.int32)
    for nm, xyz, r in zip(names, coords, rows):
        conf.SetAtomPosition(idx_of_name[nm], xyz.tolist())
        row_of_idx[idx_of_name[nm]] = r
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    return mol, row_of_idx


def build_ligand(at, lig_res: str, smiles: str):
    """``map_ligand`` + RDKit AddHs + the affinity_bench H-only relax."""
    from rdkit import Chem

    mol, row_of_idx = map_ligand(at, lig_res, smiles)
    molH = Chem.AddHs(mol, addCoords=True)
    heavy_rows = np.full(molH.GetNumAtoms(), -1, dtype=np.int32)
    heavy_rows[: len(row_of_idx)] = row_of_idx     # AddHs appends H after the heavy atoms

    _relax_hydrogens = _energy_labels("_relax_hydrogens")
    try:
        _relax_hydrogens(molH)
    except Exception as exc:  # noqa: BLE001 -- AddHs geometry is an acceptable fallback
        print(f"  ligand H relax failed ({exc}); keeping AddHs coordinates")
    return molH, heavy_rows


def _energy_labels(name: str):
    """Pull one helper out of ``analysis/energy_labels/prep_system.py``.

    Reused verbatim (not reimplemented) so hydrogen placement is bit-identical to the
    affinity_bench prep and the MLFF descriptors stay comparable across campaigns.
    """
    sys.path.insert(0, str(Path(__file__).resolve().parent))  # prep_system.py, vendored alongside
    import prep_system                                   # noqa: PLC0415
    return getattr(prep_system, name)


# --------------------------------------------------------------------- HEM CCD template
def parse_ccd_template(path: Path, code: str = HEM_RES) -> dict:
    """Single-pass text scan of components.cif for one ``data_<code>`` block.

    Returns ``{"atom_ids", "elements", "xyz" (ideal, [n,3]), "bonds" [(i, j), ...]}``.
    Same deliberate text-scan approach as ``build_cofold_inputs.py:parse_ccd_block`` (which
    only needs the atom ids); this one additionally keeps the ideal coordinates and the
    bond graph, which is all ``place_hem_hydrogens`` needs.  ~490 MB scan, so the result is
    cached to JSON and to a process-level global.
    """
    want = f"data_{code}"
    in_block = False
    atom_hdr: list[str] = []
    bond_hdr: list[str] = []
    mode = None
    atom_ids: list[str] = []
    elements: list[str] = []
    xyz: list[list[float]] = []
    bonds: list[tuple[str, str]] = []
    if not path.exists():
        raise FileNotFoundError(
            f"CCD components file {path} is missing, so the HEM hydrogen template cannot "
            f"be built. Point CCD_COMPONENTS at a components.cif, or run with "
            f"--hem_h none (heavy-atom-only heme, documented degradation).")
    with path.open() as fh:
        for ln in fh:
            line = ln.rstrip("\n")
            if line.startswith("data_"):
                if in_block:
                    break
                in_block = line.strip() == want
                continue
            if not in_block:
                continue
            if line.startswith("_chem_comp_atom."):
                atom_hdr.append(line.strip())
                mode = "atom"
                continue
            if line.startswith("_chem_comp_bond."):
                bond_hdr.append(line.strip())
                mode = "bond"
                continue
            f = line.split()
            if not f or f[0] != code:
                continue
            if mode == "atom":
                i_id = atom_hdr.index("_chem_comp_atom.atom_id")
                i_sym = atom_hdr.index("_chem_comp_atom.type_symbol")
                try:
                    ix = atom_hdr.index("_chem_comp_atom.pdbx_model_Cartn_x_ideal")
                except ValueError:
                    ix = atom_hdr.index("_chem_comp_atom.model_Cartn_x")
                if len(f) <= max(i_id, i_sym, ix + 2):
                    continue
                try:
                    p = [float(f[ix]), float(f[ix + 1]), float(f[ix + 2])]
                except ValueError:                    # '?' ideal coords -> use the model set
                    im = atom_hdr.index("_chem_comp_atom.model_Cartn_x")
                    p = [float(f[im]), float(f[im + 1]), float(f[im + 2])]
                atom_ids.append(f[i_id])
                elements.append(f[i_sym].capitalize())
                xyz.append(p)
            elif mode == "bond":
                i1 = bond_hdr.index("_chem_comp_bond.atom_id_1")
                i2 = bond_hdr.index("_chem_comp_bond.atom_id_2")
                if len(f) > max(i1, i2):
                    bonds.append((f[i1], f[i2]))
    if not atom_ids:
        raise ValueError(f"no data_{code} block (or no atoms in it) in {path}")
    return {"code": code, "atom_ids": atom_ids, "elements": elements,
            "xyz": xyz, "bonds": [list(b) for b in bonds]}


def hem_template(cache: Path | None = None) -> dict:
    """Cached CCD HEM template.  Verifies the invariants the rest of the file relies on."""
    global _HEM_TEMPLATE
    if _HEM_TEMPLATE is not None:
        return _HEM_TEMPLATE
    cache = cache or (DATA / "receptors" / "hem_template.json")
    tmpl = None
    if cache.exists():
        try:
            tmpl = json.loads(cache.read_text())
        except json.JSONDecodeError:
            print(f"  WARNING: {cache} is not valid JSON; rebuilding from the CCD")
    if tmpl is None:
        t0 = time.time()
        tmpl = parse_ccd_template(CCD_COMPONENTS)
        print(f"  CCD scan for data_{HEM_RES}: {time.time() - t0:.1f} s")
        try:
            # Atomic: 20 array tasks race to build this on their first system, and a
            # half-written JSON read by a sibling would look like a corrupt template.
            cache.parent.mkdir(parents=True, exist_ok=True)
            tmp = cache.with_suffix(f".{os.getpid()}.tmp")
            tmp.write_text(json.dumps(tmpl))
            os.replace(tmp, cache)
            print(f"  cached HEM template -> {cache}")
        except OSError as exc:
            print(f"  WARNING: could not cache the HEM template ({exc})")
    heavy = [a for a, e in zip(tmpl["atom_ids"], tmpl["elements"]) if e != "H"]
    if len(heavy) != HEM_N_HEAVY:
        raise ValueError(f"HEM template has {len(heavy)} heavy atoms, expected {HEM_N_HEAVY}")
    if HEM_FE_ATOM not in tmpl["atom_ids"]:
        raise ValueError(f"HEM template has no atom named {HEM_FE_ATOM!r}")
    for a in HEM_ACID_H:
        if a not in tmpl["atom_ids"]:
            raise ValueError(
                f"HEM template has no atom {a!r}; HEM_ACID_H is what makes "
                f"HEM_CHARGE 'propionate' = -1 correct, so it must not be stale")
    _HEM_TEMPLATE = tmpl
    return tmpl


def _kabsch(P: np.ndarray, Q: np.ndarray):
    """Rigid transform with ``Q ~ P @ R.T + t``.  Returns ``(R, t, rmsd)`` or None if the
    anchor set is degenerate (collinear), which would leave the rotation undetermined."""
    pc, qc = P.mean(0), Q.mean(0)
    A, B = P - pc, Q - qc
    if np.linalg.svd(A, compute_uv=False)[1] < 1e-3:      # 2nd singular value ~ 0
        return None
    U, _, Vt = np.linalg.svd(A.T @ B)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    t = qc - R @ pc
    rmsd = float(np.sqrt((((A @ R.T) - B) ** 2).sum(axis=1).mean()))
    return R, t, rmsd


def place_hem_hydrogens(names, xyz, tmpl, protonation: str = "propionate"):
    """Hydrogens for the heme, by rigid local-frame transfer from the CCD ideal geometry.

    amber14-all.xml has no HEM template and RDKit cannot reliably sanitise an
    Fe-porphyrin, so neither of the two protonation paths used elsewhere applies.  For each
    template hydrogen h with parent heavy atom p we Kabsch-fit {p} + p's heavy neighbours
    (extending to the 2nd shell when that is fewer than 3 anchors) from the ideal geometry
    onto the pose, and carry h through the same rigid transform.  Bond lengths and angles
    are therefore EXACTLY ideal; only the 2nd-shell torsion (methyl rotamer, propionate
    tail) is approximate, which is the same accuracy class as RDKit's AddHs.

    ``protonation`` selects which template hydrogens are dropped -- and the same table
    fixes the formal charge, so the two can never disagree (TRAP 3).
    """
    drop, charge = HEM_PROTONATION[protonation]
    if drop is None:
        return np.zeros(0, dtype=object), np.zeros((0, 3)), {"placed": 0, "failed": 0,
                                                             "max_rmsd": 0.0, "charge": charge}
    dropset = set(drop)
    tid = {a: i for i, a in enumerate(tmpl["atom_ids"])}
    txyz = np.asarray(tmpl["xyz"], dtype=np.float64)
    elem = dict(zip(tmpl["atom_ids"], tmpl["elements"]))
    pos_of = {str(n): i for i, n in enumerate(names)}

    unknown = sorted(set(pos_of) - set(tid))
    if unknown:
        raise ValueError(
            f"pose HEM carries atom names absent from the CCD template: {unknown[:6]} -- "
            f"either the pose is not a CCD_HEM entity or the template is stale")

    nbrs: dict[str, list[str]] = {a: [] for a in tmpl["atom_ids"]}
    for a, b in tmpl["bonds"]:
        if a in nbrs and b in nbrs:
            nbrs[a].append(b)
            nbrs[b].append(a)
    heavy_nbrs = {a: [n for n in v if elem.get(n) != "H"] for a, v in nbrs.items()}

    out_xyz, failed, max_rmsd = [], [], 0.0
    for h in tmpl["atom_ids"]:
        if elem[h] != "H" or h in dropset:
            continue
        parents = [p for p in nbrs[h] if elem.get(p) != "H"]
        if len(parents) != 1 or parents[0] not in pos_of:
            failed.append(h)
            continue
        p = parents[0]
        anchors = [p] + [n for n in heavy_nbrs[p] if n in pos_of]
        if len(anchors) < 3:
            for n in heavy_nbrs[p]:
                for m in heavy_nbrs.get(n, ()):
                    if m in pos_of and m not in anchors:
                        anchors.append(m)
        if len(anchors) < 3:
            failed.append(h)
            continue
        fit = _kabsch(txyz[[tid[a] for a in anchors]],
                      np.asarray([xyz[pos_of[a]] for a in anchors], dtype=np.float64))
        if fit is None:
            failed.append(h)
            continue
        R, t, rmsd = fit
        max_rmsd = max(max_rmsd, rmsd)
        out_xyz.append(R @ txyz[tid[h]] + t)

    xyz_h = np.asarray(out_xyz, dtype=np.float64).reshape(-1, 3)
    info = {"placed": len(out_xyz), "failed": len(failed), "failed_ids": failed[:8],
            "max_rmsd": max_rmsd, "charge": charge}
    print(f"  HEM protonation [{protonation}]: {info['placed']} H placed, "
          f"{info['failed']} failed, max local-frame RMSD {max_rmsd:.3f} A, "
          f"formal charge {charge:+d}")
    if failed:
        print(f"    WARNING: unplaced HEM hydrogens {failed[:8]} -- the heme is left with "
              f"dangling valences at those atoms")
    return np.full(len(out_xyz), "H", dtype=object), xyz_h, info


# -------------------------------------------------------------------------------- pocket
def build_pocket(at, prot_mask, lig_heavy_xyz, fe_xyz=None, cys_thiolate: str = "auto",
                 fe_sg_max: float = 4.0, hem_heavy_xyz=None):
    """Whole amino-acid residues within POCKET_CUT of the POCKET_REF reference set.

    POCKET_REF (module-level, env-overridable) chooses that reference set: the DEFAULT
    ``substrate_hem`` is the substrate-or-heme union described below, and ``substrate``
    measures from the substrate heavy atoms alone for the small-cluster arm.  POCKET_FORCE_RES
    additionally keeps named auth_seq residues regardless of distance -- as whole residues
    through this same PDBFixer/ff14SB path, so ``charge_ff`` covers them -- and RAISES, naming
    the residue, if one is absent from the protein chain.

    The default reference set is the UNION of the substrate and the heme, because ``prep_one``
    keeps the whole heme unconditionally: POCKET_CUT = 12 A is justified as MACE-POLAR's
    receptive field (r_max 6.0 x 2 layers), i.e. every retained atom has a COMPLETE
    environment -- but that only holds for the atoms the cut is MEASURED FROM.  Measuring
    from the substrate alone would hand POLAR a heme periphery (and hence the headline Fe
    row, TRAP 4) with a truncated environment.  The extra atoms/residues are printed as an
    explicit delta because POCKET_CUT is the stage-60 memory lever.

    HEM IS DELIBERATELY NOT IN HERE (as ATOMS -- only as a distance reference).
    amber14-all.xml has no HEM template, so ``createSystem`` would raise, the reference
    code would catch it and fall back with ``pocket_charge = 0`` -- and a wrong total
    charge silently changes every POLAR descriptor.  The heme is protonated separately by
    ``place_hem_hydrogens``.

    Returns ``(symbols, res_names, coords_A, charge_ff, charge_ok, cys)`` where
    ``charge_ok`` is False when the ff14SB relax raised and ``charge_ff`` is therefore a
    fallback 0 rather than a measured charge, and ``cys`` describes the axial cysteine and
    the thiolate bookkeeping.  ``res_names`` exists because stages 60/65 accept
    ``res_name == "HEM"`` as a fallback heme mask, and a heme mask that cannot be
    recovered from the npz is trap 1 waiting to happen.
    """
    import io

    from openmm.unit import angstrom
    from pdbfixer import PDBFixer

    if not prot_mask.any():
        raise ValueError("no protein atoms found")
    n_lig_ref = len(lig_heavy_xyz)
    ref_xyz = np.asarray(lig_heavy_xyz, dtype=np.float64)
    if hem_heavy_xyz is not None and len(hem_heavy_xyz):
        ref_xyz = np.vstack([ref_xyz, np.asarray(hem_heavy_xyz, dtype=np.float64)])
    d = np.linalg.norm(at["coord"][prot_mask][:, None, :] - ref_xyz[None, :, :], axis=-1)
    p_asym, p_seq = at["asym"][prot_mask], at["auth_seq"][prot_mask]

    def _res_set(near_mask):
        return set(zip(p_asym[near_mask], p_seq[near_mask]))

    def _atom_sel(res_set):
        return np.array([(a, s) in res_set
                         for a, s in zip(at["asym"], at["auth_seq"])], dtype=bool) & prot_mask

    # POCKET_REF chooses WHICH COLUMNS of d the cut is measured from; the substrate columns
    # are the first n_lig_ref by construction above.  "substrate_hem" is `d.min(axis=1)`,
    # i.e. byte-for-byte the previous behaviour.
    d_ref = d if POCKET_REF == "substrate_hem" else d[:, :n_lig_ref]
    keep_dist = _res_set(d_ref.min(axis=1) <= POCKET_CUT)
    if not keep_dist:
        raise ValueError(
            f"no protein residue within POCKET_CUT={POCKET_CUT} A of the "
            f"{'substrate/heme' if POCKET_REF == 'substrate_hem' else 'substrate-only'} "
            f"reference set")
    keep_res = set(keep_dist)

    # --- POCKET_FORCE_RES: whole residues kept regardless of distance -------------------
    # They join keep_res BEFORE the PDB round-trip, so they go through exactly the same
    # whole-residue selection, PDBFixer protonation and ff14SB charge sum as everything
    # else and pocket_charge stays correct.  An ABSENT residue is fatal rather than
    # skipped: a 5 A crop that quietly lost Cys442 leaves the Fe coordination sphere
    # incomplete and still returns a plausible energy, which is unfalsifiable downstream.
    forced, forced_new = [], []
    if POCKET_FORCE_RES_LIST:
        p_res = at["res_name"][prot_mask]
        for num in POCKET_FORCE_RES_LIST:
            hit = np.asarray(p_seq) == num
            if not hit.any():
                raise ValueError(
                    f"POCKET_FORCE_RES={POCKET_FORCE_RES}: auth_seq residue {num} is ABSENT "
                    f"from the protein chain of this pose, so it cannot be force-included. "
                    f"It is on that list because its absence is not survivable -- a missing "
                    f"Cys442 leaves the Fe coordination sphere incomplete while still "
                    f"producing a plausible energy. Refusing to emit a silently truncated "
                    f"pocket.")
            pairs = set(zip(np.asarray(p_asym)[hit], np.asarray(p_seq)[hit]))
            label = f"{str(np.asarray(p_res)[hit][0])}{num}"
            if pairs - keep_res:
                forced_new.append(label)
            keep_res |= pairs
            forced.append(label)

    sel = _atom_sel(keep_res)
    n_forced_atoms = (int(sel.sum()) - int(_atom_sel(keep_dist).sum())) if forced else 0
    print(f"  pocket: {int(sel.sum())} heavy atoms in {len(keep_res)} residues "
          f"(cut {POCKET_CUT} A, POCKET_REF={POCKET_REF}, "
          f"POCKET_FORCE_RES={POCKET_FORCE_RES or 'none'})")
    if forced:
        print(f"    POCKET_FORCE_RES kept {len(forced)} residue(s) {forced} regardless of "
              f"distance; {len(forced_new)} of them were NOT already inside the "
              f"{POCKET_CUT} A cut {forced_new} (+{n_forced_atoms} heavy atoms)")
    if len(ref_xyz) > n_lig_ref:
        # The substrate-only vs substrate+HEM delta, printed in BOTH modes (POCKET_CUT is
        # the stage-60 memory lever, so the size of the heme's contribution has to stay
        # visible even when the cut is not measured from it).  Both figures are
        # distance-only; POCKET_FORCE_RES is reported separately above.
        keep_l = _res_set(d[:, :n_lig_ref].min(axis=1) <= POCKET_CUT)
        keep_u = _res_set(d.min(axis=1) <= POCKET_CUT)
        sel_l, sel_u = _atom_sel(keep_l), _atom_sel(keep_u)
        print(f"    reference = substrate + HEM heavy atoms: {int(sel_l.sum())} -> "
              f"{int(sel_u.sum())} atoms (+{int(sel_u.sum()) - int(sel_l.sum())}), "
              f"{len(keep_l)} -> {len(keep_u)} residues "
              f"(+{len(keep_u) - len(keep_l)}) vs a substrate-only cut "
              f"[POCKET_REF={POCKET_REF} is the one IN FORCE]")

    # One temp single-char PDB chain per label_asym_id. The reference hardcodes 'A', which
    # would silently MERGE residues if a pocket ever spanned two chains (it does not today
    # -- one proteinChain entity -- but a merged residue is invisible until the FF explodes).
    asyms = sorted({str(a) for a in at["asym"][sel]})
    if len(asyms) > len(TMP_CHAINS):
        raise ValueError(f"pocket spans {len(asyms)} chains; too many for the PDB round-trip")
    chain_of = {a: TMP_CHAINS[i] for i, a in enumerate(asyms)}

    lines, src_rows = [], []
    for i, r in enumerate(np.where(sel)[0]):
        nm, rn = at["atom_name"][r], at["res_name"][r]
        x, y, z = at["coord"][r]
        lines.append(
            f"ATOM  {i+1:5d} {nm:^4s} {rn:>3s} {chain_of[str(at['asym'][r])]}"
            f"{at['auth_seq'][r]:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {at['element'][r]:>2s}")
        src_rows.append(int(r))
    lines.append("END")

    fixer = PDBFixer(pdbfile=io.StringIO("\n".join(lines)))
    fixer.findMissingResidues()
    fixer.missingResidues = {}       # never build unresolved loops
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(PH)

    syms = np.array([a.element.symbol if a.element is not None else "X"
                     for a in fixer.topology.atoms()], dtype=object)
    atom_meta = [(a.residue.name, int(a.residue.id), a.name)
                 for a in fixer.topology.atoms()]
    res_names = np.array([m[0] for m in atom_meta], dtype=object)
    _minimize_frozen_heavy = _energy_labels("_minimize_frozen_heavy")
    charge_ok = True
    try:
        xyz, charge_ff = _minimize_frozen_heavy(fixer.topology, fixer.positions)
        xyz = np.asarray(xyz, dtype=np.float64)
    except Exception as exc:  # noqa: BLE001
        # LOUD, not silent: the reference swallows this and reports charge 0, which is a
        # charge bug rather than a geometry degradation.  charge_ok=False is what makes the
        # fallback 0 DISTINGUISHABLE from a genuinely neutral pocket downstream -- without
        # it a ~60-residue CYP pocket that should carry, say, -4 looks identical to one
        # that really is 0, and stage 60 feeds q_lig+q_pock+q_hem straight into POLAR.
        print(f"  WARNING: pocket H relax FAILED ({type(exc).__name__}: {exc}); keeping "
              f"PDBFixer geometry AND falling back to charge 0 -- the POLAR conditioning "
              f"charge for this system is NOT trustworthy (pocket_charge_ok=False)")
        xyz = np.array(fixer.positions.value_in_unit(angstrom), dtype=np.float64)
        charge_ff = 0
        charge_ok = False
    if xyz.shape != (len(syms), 3):
        raise ValueError(f"pocket coord shape {xyz.shape} != ({len(syms)}, 3)")

    # --- the proximal cysteine thiolate (TRAP 3) ---------------------------------------
    cys = {"fe_sg_dist": float("nan"), "resid": -1, "thiolate": False,
           "charge": 0, "n_hg_dropped": 0}
    if fe_xyz is not None and cys_thiolate != "off":
        sg = [i for i, (rn, _, an) in enumerate(atom_meta)
              if rn == "CYS" and an == CYS_SG_ATOM]
        if sg:
            dists = np.linalg.norm(xyz[sg] - np.asarray(fe_xyz)[None, :], axis=1)
            j = int(np.argmin(dists))
            cys["fe_sg_dist"] = float(dists[j])
            cys["resid"] = atom_meta[sg[j]][1]
            if cys_thiolate == "on" or dists[j] <= fe_sg_max:
                # Drop the HG that PDBFixer put on the coordinating SG: it does not exist
                # (the sulfur's lone pair is on the iron) and it sits inside the Fe.
                keep = np.ones(len(syms), bool)
                for i, (rn, rid, an) in enumerate(atom_meta):
                    if (syms[i] == "H" and rn == "CYS" and rid == cys["resid"]
                            and np.linalg.norm(xyz[i] - xyz[sg[j]]) < 1.6):
                        keep[i] = False
                n_drop = int((~keep).sum())
                # The -1 is the charge RELATIVE to what ff14SB just reported, and ff14SB
                # reported a neutral Cys-SH. If no HG was actually there to drop, the
                # sulfur was not counted as protonated in the first place and booking -1
                # would invent a charge -- so the correction is tied to the deletion.
                cys.update(thiolate=n_drop > 0, charge=-1 if n_drop else 0,
                           n_hg_dropped=n_drop)
                syms, xyz = syms[keep], xyz[keep]
                res_names = res_names[keep]
                if n_drop:
                    print(f"  axial CYS{cys['resid']} SG-Fe {dists[j]:.2f} A -> THIOLATE: "
                          f"dropped {n_drop} H, charge {cys['charge']:+d}")
                else:
                    print(f"  WARNING: axial CYS{cys['resid']} SG-Fe {dists[j]:.2f} A but "
                          f"PDBFixer placed no HG on it; NOT booking the thiolate -1, "
                          f"because ff14SB did not count a proton there either")
            else:
                print(f"  nearest CYS{cys['resid']} SG-Fe {dists[j]:.2f} A > "
                      f"{fe_sg_max} A -> NOT treated as a thiolate")
        else:
            print("  no CYS SG in the pocket crop -- no thiolate correction")

    print(f"  pocket protonated: {len(syms)} atoms, ff14SB charge {charge_ff:+d}"
          f"{'' if charge_ok else ' (UNTRUSTED FALLBACK)'}")
    # psel records WHICH pocket definition produced this crop, so a scored CSV can never be
    # confused about whether it came from the 12 A union or the 5 A substrate cluster.
    psel = {"ref": POCKET_REF, "force_res": POCKET_FORCE_RES,
            "forced": list(forced), "forced_new": list(forced_new),
            "n_forced_atoms": int(n_forced_atoms),
            "n_res": len(keep_res), "n_res_dist": len(keep_dist)}
    return syms, res_names, xyz, int(charge_ff), bool(charge_ok), cys, psel


# ----------------------------------------------------------------------------- prep_one
def prep_one(cif, smiles, out_npz, *, pair_id="", isoform="", lig_res=None,
             hem_h="ccd", hem_protonation="propionate", hem_charge=None,
             hem_spin=HEM_SPIN_MULT, cys_thiolate="auto", fe_sg_max=4.0,
             label_fe=True, pose_rank=-1, n_poses=-1,
             has_hem_requested=True, fe_cys_bond_requested=True) -> dict:
    """One co-folded CIF -> one prep npz.  Importable so a shard driver pays imports once."""
    # Seeded PER POSE rather than once per shard, so a pose's hydrogens depend only on PREP_SEED
    # and not on how many poses preceded it in the shard -- otherwise re-prepping one pose alone
    # would not reproduce the same pose inside an array.
    if PREP_SEED_VAL is not None:
        import random  # noqa: PLC0415 -- only needed on the opt-in path
        random.seed(PREP_SEED_VAL)
    print(f"[in ] {cif}")
    at = read_atom_site(cif)
    print(f"  _atom_site rows: {len(at['row'])}")
    g = split_groups(at, lig_res)
    lig_res = g["lig_res"]

    # --- substrate -----------------------------------------------------------------------
    molH, lig_rows = build_ligand(at, lig_res, smiles)
    lconf = molH.GetConformer()
    lig_sym = np.array([a.GetSymbol() for a in molH.GetAtoms()], dtype=object)
    lig_xyz = np.array([list(lconf.GetAtomPosition(i)) for i in range(molH.GetNumAtoms())])
    lig_charge = int(sum(a.GetFormalCharge() for a in molH.GetAtoms()))
    lig_heavy_xyz = lig_xyz[[a.GetIdx() for a in molH.GetAtoms() if a.GetAtomicNum() > 1]]
    print(f"  substrate: {molH.GetNumAtoms()} atoms with H (formal charge {lig_charge:+d})")

    # --- heme ----------------------------------------------------------------------------
    hem_sym = np.zeros(0, dtype=object)
    hem_xyz = np.zeros((0, 3))
    hem_rows = np.zeros(0, dtype=np.int32)
    hem_q, hem_info = 0, {"placed": 0, "failed": 0, "max_rmsd": 0.0}
    hem_heavy_xyz = None
    fe_xyz = None
    fe_row = -1
    if g["has_hem"]:
        hm = g["hem"]
        hem_names = [str(n) for n in at["atom_name"][hm]]
        hem_heavy_xyz = at["coord"][hm]
        hem_rows = at["row"][hm].astype(np.int32)
        fe_local = int(np.argmax(g["fe"][hm]))
        fe_xyz = hem_heavy_xyz[fe_local]
        fe_row = int(at["row"][g["fe"]][0])
        dmin = float(np.linalg.norm(
            hem_heavy_xyz[:, None, :] - lig_heavy_xyz[None, :, :], axis=-1).min())
        print(f"  HEM: {len(hem_names)} heavy atoms, closest heme-substrate contact "
              f"{dmin:.2f} A (Fe row {fe_row})")
        if dmin > POCKET_CUT:
            print(f"    WARNING: heme is further than POCKET_CUT={POCKET_CUT} A from the "
                  f"substrate; it is kept anyway (a CYP without its heme is not a CYP)")
        mode = hem_protonation if hem_h == "ccd" else "none"
        h_sym, h_xyz, hem_info = place_hem_hydrogens(
            hem_names, hem_heavy_xyz, hem_template(), mode)
        hem_sym = np.concatenate([np.array(
            [at["element"][r] for r in at["row"][hm]], dtype=object), h_sym])
        hem_xyz = np.concatenate([hem_heavy_xyz, h_xyz])
        hem_rows = np.concatenate([hem_rows, np.full(len(h_sym), -1, np.int32)])
        if hem_info["failed"]:
            # TRAP 3: hem_charge is derived for exactly C34H30FeN4O4(-). A heme missing
            # some of its hydrogens has a different electron count, so the charge and the
            # geometry disagree -- the very thing the single HEM_PROTONATION table exists
            # to prevent. The shard driver logs this pair and carries on.
            raise ValueError(
                f"HEM protonation incomplete: {hem_info['failed']} template hydrogen(s) "
                f"could not be placed {hem_info.get('failed_ids', [])}, so the geometry is "
                f"a partially protonated heme while the {mode}-implied charge "
                f"{hem_info['charge']} is derived for the fully protonated one. Refusing "
                f"to emit a system whose charge and geometry disagree.")
        if hem_info["charge"] is None:
            # --hem_h none strips ALL 32 CCD hydrogens; no tabulated charge applies to the
            # bare C34FeN4O4 skeleton, so it must be stated rather than inherited.
            if hem_charge is None:
                raise ValueError(
                    "--hem_h none leaves a heavy-atom-only heme (32 H removed), for which "
                    "the HEM_PROTONATION charges (derived for C34H30FeN4O4) are wrong. "
                    "Pass an explicit --hem_charge for this species.")
            hem_q = int(hem_charge)
            print(f"  HEM charge {hem_q:+d} taken from --hem_charge (heavy-atom-only heme "
                  f"has no tabulated charge)")
        else:
            hem_q = int(hem_charge) if hem_charge is not None else int(hem_info["charge"])
            if hem_charge is not None and hem_charge != hem_info["charge"]:
                print(f"  WARNING: --hem_charge {hem_charge:+d} OVERRIDES the "
                      f"{hem_protonation}-implied {hem_info['charge']:+d}")
    else:
        msg = ("no HEM residue in the pose. CONTEXT.md risk 1 is exactly this: a CYP "
               "co-fold without its heme has no active-site physics to measure.")
        if has_hem_requested:
            raise ValueError(msg + " The manifest says has_hem_requested=True, so this "
                                   "pose is a co-folding failure, not an ablation.")
        print(f"  {msg} (has_hem_requested=False -> ablation arm, continuing)")

    # --- pocket --------------------------------------------------------------------------
    pock_sym, pock_res, pock_xyz, charge_ff, charge_ok, cys, psel = build_pocket(
        at, g["protein"], lig_heavy_xyz, fe_xyz=fe_xyz,
        cys_thiolate=cys_thiolate, fe_sg_max=fe_sg_max, hem_heavy_xyz=hem_heavy_xyz)

    # --- assemble ------------------------------------------------------------------------
    symbols = np.concatenate([pock_sym, hem_sym, lig_sym])
    coords = np.concatenate([pock_xyz, hem_xyz, lig_xyz])
    n_p, n_h, n_l = len(pock_sym), len(hem_sym), len(lig_sym)
    group = np.concatenate([np.full(n_p, G_PROTEIN, np.int8),
                            np.full(n_h, G_HEM, np.int8),
                            np.full(n_l, G_SUBSTRATE, np.int8)])
    # TRAP 1: is_ligand is the SUBSTRATE alone. The heme rides with the receptor.
    is_ligand = (group == G_SUBSTRATE)
    is_hem = (group == G_HEM)
    res_name = np.concatenate([pock_res,
                               np.full(n_h, HEM_RES, dtype=object),
                               np.full(n_l, lig_res, dtype=object)])
    prot_index = np.concatenate([np.full(n_p, -1, np.int32), hem_rows, lig_rows])
    is_fe = np.zeros(len(symbols), bool)
    fe_index = -1
    if fe_row >= 0:
        hit = np.where((group == G_HEM) & (prot_index == fe_row))[0]
        if len(hit):
            fe_index = int(hit[0])
            is_fe[fe_index] = True

    # TRAP 1 assertions -- a violation poisons every energy feature, so it is fatal. These
    # are `raise`, not `assert`: `python -O` strips asserts and the whole point is that this
    # check can never be optimised away.
    for bad, msg in (
        (bool(is_fe[is_ligand].any()), "ligand group contains the heme Fe atom"),
        (bool((group[is_ligand] == G_HEM).any()), "ligand group contains HEM atoms"),
        ("Fe" in set(symbols[is_ligand]), "ligand group contains an Fe atom"),
        (int(is_ligand.sum()) != n_l or n_l == 0, "empty or over-wide ligand group"),
    ):
        if bad:
            raise ValueError(
                f"is_ligand mask is WRONG for {pair_id or cif}: {msg}. "
                f"model/mlff/backends/base.py:73-80 slices the subsystems on this mask, so "
                f"E_int would silently become 'drug+heme vs apo protein'. Refusing to emit.")

    # TRAP 4: the Fe gets its OWN descriptor row; the 42 other heme atoms do not, because
    # the heme is identical across every complex of an isoform and would just add
    # near-constant dimensions.
    labelled = (prot_index >= 0) & (is_ligand | (is_fe if label_fe else False))

    # TRAP 3: pocket_charge = amino acids + thiolate, heme EXCLUDED. Stages 60/65 sum
    # ligand + pocket + hem themselves, so folding the heme in here would double-count it.
    pocket_charge = charge_ff + int(cys["charge"])
    total_charge = lig_charge + pocket_charge + hem_q

    n_atoms = len(symbols)
    print(f"  charge: substrate {lig_charge:+d} | ff14SB pocket {charge_ff:+d} | "
          f"Cys thiolate {cys['charge']:+d} -> pocket_charge {pocket_charge:+d} "
          f"(heme EXCLUDED) | HEM {hem_q:+d}  =>  TOTAL {total_charge:+d}")
    if n_atoms > MAX_ATOMS:
        raise ValueError(f"system too large: {n_atoms} > MAX_ATOMS={MAX_ATOMS}")
    if n_atoms > WARN_ATOMS:
        print(f"  WARNING: {n_atoms} atoms > WARN_ATOMS={WARN_ATOMS}. POLAR's long-range "
              f"module is superlinear in atom count (707 atoms measured 22.6 GB); lower "
              f"POCKET_CUT (currently {POCKET_CUT} A) if the extract stage OOMs.")

    # Cheap clash audit -- PDBFixer never saw the heme, so a protein H can be placed inside
    # it, and a bad contact makes the MLFF energies garbage rather than merely noisy.
    min_gap = float("inf")
    if n_h and (n_p or n_l):
        other = np.concatenate([pock_xyz, lig_xyz]) if n_p and n_l else (
            pock_xyz if n_p else lig_xyz)
        min_gap = float(np.linalg.norm(
            hem_xyz[:, None, :] - other[None, :, :], axis=-1).min())
        if min_gap < 0.6:
            raise ValueError(
                f"heme/environment clash: closest contact {min_gap:.2f} A. The pose or the "
                f"protonation is broken; refusing to emit a system whose MLFF energy would "
                f"be meaningless.")
        if min_gap < 0.9:
            print(f"  WARNING: closest heme-environment contact {min_gap:.2f} A")

    out_npz = str(out_npz)
    Path(out_npz).parent.mkdir(parents=True, exist_ok=True)
    # Atomic: the shard driver resumes on plain EXISTENCE, so a task killed at the walltime
    # mid-write would otherwise leave a truncated npz that every resubmission skips and the
    # pair would be lost for good. Same dance as the HEM template cache above. The temp name
    # ends in .npz so numpy does not append a second suffix, and it does NOT end in
    # `_prep.npz` so stages 60/65 (`glob("*_prep.npz")`) never pick a partial file up.
    tmp_npz = f"{out_npz}.{os.getpid()}.tmp.npz"
    np.savez_compressed(
        tmp_npz,
        # --- geometry / topology ---------------------------------------------------------
        symbols=symbols.astype("U2"),
        z=np.array([_Z.get(str(s), 0) for s in symbols], dtype=np.int16),
        coords=coords.astype(np.float32),
        is_ligand=is_ligand, is_hem=is_hem, group=group,
        res_name=res_name.astype("U8"),
        group_names=np.array(GROUP_NAMES, dtype="U9"),
        labelled=labelled, label_group=group.copy(),
        protenix_atom_index=prot_index.astype(np.int32),
        is_fe=is_fe,
        fe_index=np.int32(fe_index), fe_protenix_index=np.int32(fe_row),
        # --- charge / spin ---------------------------------------------------------------
        ligand_charge=np.int32(lig_charge),
        pocket_charge=np.int32(pocket_charge),        # = ff14SB + thiolate, NO heme
        pocket_charge_ff=np.int32(charge_ff),
        # False => the ff14SB relax raised and pocket_charge_ff is a FALLBACK 0, not a
        # measured charge. Without this flag an untrusted pocket is indistinguishable from
        # a genuinely neutral one (FINDING 1); stage 60 refuses on it.
        pocket_charge_ok=np.bool_(charge_ok),
        hem_charge=np.int32(hem_q),
        cys_thiolate_charge=np.int32(cys["charge"]),
        total_charge=np.int32(total_charge),
        hem_spin_mult=np.int32(hem_spin if g["has_hem"] else 1),
        spin_hint_complex=np.int32(hem_spin if g["has_hem"] else 1),
        spin_hint_protein=np.int32(hem_spin if g["has_hem"] else 1),
        spin_hint_ligand=np.int32(1),
        # --- provenance ------------------------------------------------------------------
        pair_id=np.array(str(pair_id)), isoform=np.array(str(isoform)),
        smiles=np.array(str(smiles)), source_cif=np.array(str(cif)),
        lig_res=np.array(str(lig_res)),
        pose_rank=np.int32(pose_rank), n_poses=np.int32(n_poses),
        pocket_cut=np.float32(POCKET_CUT), ph=np.float32(PH),
        # --- pocket-definition knobs (ADDITIVE; the defaults reproduce the 12 A union) ---
        pocket_ref=np.array(str(psel["ref"])),
        pocket_force_res=np.array(str(psel["force_res"])),
        pocket_force_res_included=np.array(",".join(psel["forced"])),
        n_pocket_force_res=np.int32(len(psel["forced"])),
        n_pocket_force_res_outside_cut=np.int32(len(psel["forced_new"])),
        n_pocket_force_atoms=np.int32(psel["n_forced_atoms"]),
        n_pocket_res=np.int32(psel["n_res"]),
        n_pocket_res_within_cut=np.int32(psel["n_res_dist"]),
        # -1 = unseeded (the historical default): the pocket hydrogens of this npz are ONE
        # stochastic realisation and re-prepping will not reproduce them.
        prep_seed=np.int32(-1 if PREP_SEED_VAL is None else PREP_SEED_VAL),
        has_hem=bool(g["has_hem"]),
        has_hem_requested=bool(has_hem_requested),
        fe_cys_bond_requested=bool(fe_cys_bond_requested),
        hem_h_mode=np.array(str(hem_h)),
        hem_protonation=np.array(str(hem_protonation if hem_h == "ccd" else "none")),
        n_hem_h_placed=np.int32(hem_info["placed"]),
        n_hem_h_failed=np.int32(hem_info["failed"]),
        hem_h_max_rmsd=np.float32(hem_info["max_rmsd"]),
        fe_sg_dist=np.float32(cys["fe_sg_dist"]),
        cys_thiolate_resid=np.int32(cys["resid"]),
        min_hem_contact=np.float32(min_gap),
        n_protein=np.int32(n_p), n_hem=np.int32(n_h), n_ligand=np.int32(n_l),
        n_pocket=np.int32(n_p),      # affinity_bench-era alias for n_protein
        n_atoms=np.int32(n_atoms),
        code_version=np.array(CODE_VERSION),
    )
    os.replace(tmp_npz, out_npz)
    print(f"[out] {out_npz}  total {n_atoms} atoms ({n_p} protein + {n_h} HEM + "
          f"{n_l} substrate), {int(labelled.sum())} labelled "
          f"({int((labelled & is_ligand).sum())} substrate + "
          f"{int((labelled & is_fe).sum())} Fe)")
    return {"pair_id": pair_id, "n_atoms": n_atoms, "n_protein": n_p, "n_hem": n_h,
            "n_ligand": n_l, "n_labelled": int(labelled.sum()),
            "has_hem": bool(g["has_hem"]), "total_charge": total_charge,
            "pocket_charge_ok": bool(charge_ok),
            "fe_sg_dist": cys["fe_sg_dist"], "min_hem_contact": min_gap,
            "n_hem_h_failed": hem_info["failed"], "source_cif": str(cif)}


_Z = {"H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9,
      "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17,
      "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25,
      "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Se": 34, "Br": 35, "I": 53}


# ------------------------------------------------------------------------- shard driver
def load_work(args):
    """manifest.csv (preferred) or pairs.csv + pred_root -> a DataFrame of work.

    NOTHING is inferred silently: if the manifest is absent the message names the exact
    path and the stage that writes it.  An empty work list is an error, never an empty
    output directory.
    """
    import pandas as pd

    if args.manifest is not None:
        if not args.manifest.exists():
            raise SystemExit(
                f"FATAL: manifest {args.manifest} does not exist.\n"
                f"       It is written by stage 30 (analysis/cyp_challenge/"
                f"build_cofold_inputs.py, via slurm/cyp_challenge/30_build_cofold_inputs"
                f".sbatch) as data/cyp_challenge/cofold_inputs/<TAG>/manifest.csv.\n"
                f"       Alternatively pass --pairs_csv and --pred_root.")
        df = pd.read_csv(args.manifest)
        need = {"pair_id", "smiles"}
        if not need.issubset(df.columns):
            raise SystemExit(f"FATAL: {args.manifest} lacks columns {need - set(df.columns)}")
        if "pose_glob" not in df.columns:
            if args.pred_root is None:
                raise SystemExit(
                    f"FATAL: {args.manifest} has no pose_glob column and --pred_root was "
                    f"not given, so the co-folded CIFs cannot be located.")
            df["pose_glob"] = [str(args.pred_root / p / f"seed_{args.seed}" / "predictions"
                                   / f"{p}_sample_*.cif") for p in df["pair_id"]]
    else:
        if args.pairs_csv is None or args.pred_root is None:
            raise SystemExit("FATAL: pass --manifest, or both --pairs_csv and --pred_root")
        if not args.pairs_csv.exists():
            raise SystemExit(f"FATAL: --pairs_csv {args.pairs_csv} does not exist")
        df = pd.read_csv(args.pairs_csv)
        df["pose_glob"] = [str(args.pred_root / p / f"seed_{args.seed}" / "predictions"
                               / f"{p}_sample_*.cif") for p in df["pair_id"]]
    for c, default in (("isoform", ""), ("has_hem_requested", True),
                       ("fe_cys_bond_requested", True)):
        if c not in df.columns:
            df[c] = default
    return df.sort_values("pair_id").reset_index(drop=True)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--single", action="store_true", help="prep one CIF and exit")
    p.add_argument("--cif", type=Path, help="--single: the pose CIF")
    p.add_argument("--out", type=Path, help="--single: the output npz")
    p.add_argument("--smiles", help="--single: the substrate SMILES")

    p.add_argument("--manifest", type=Path, default=None,
                   help="data/cyp_challenge/cofold_inputs/<TAG>/manifest.csv (stage 30)")
    p.add_argument("--pairs_csv", type=Path, default=None, help="fallback work list")
    p.add_argument("--pred_root", type=Path, default=None,
                   help="/scratch-shared/.../cyp_challenge/predictions/<TAG>")
    p.add_argument("--out_dir", type=Path, default=None)
    p.add_argument("--seed", type=int, default=101)
    p.add_argument("--shard", type=int, default=0)
    p.add_argument("--nshards", type=int, default=1)
    p.add_argument("--limit", type=int, default=0)
    p.add_argument("--overwrite", action="store_true")

    p.add_argument("--lig_res", default=None,
                   help="force the substrate res_name instead of auto-detecting l01..l99")
    p.add_argument("--pose", type=int, default=-1,
                   help="sample index to prep; -1 = the lowest available index (matches "
                        "affinity_bench, and after PRUNE=1 only one sample survives)")
    p.add_argument("--hem_h", choices=["ccd", "none"], default="ccd",
                   help="'none' leaves the heme as a bare heavy-atom skeleton and REQUIRES "
                        "an explicit --hem_charge (the tabulated charges are derived for "
                        "the fully protonated C34H30FeN4O4)")
    p.add_argument("--hem_protonation", choices=["propionate", "acid"], default="propionate",
                   help="propionate (pH 7.4, net -1) or neutral diacid (net +1)")
    p.add_argument("--hem_charge", type=int, default=None,
                   help="override the protonation-implied HEM formal charge (prints a "
                        "warning; nothing about the charge is ever defaulted silently)")
    p.add_argument("--hem_spin", type=int, default=HEM_SPIN_MULT,
                   help="heme spin multiplicity hint for the Orb path; 6 = substrate-bound "
                        "high-spin ferric, 2 = low-spin (water- or azole-coordinated)")
    p.add_argument("--cys_thiolate", choices=["auto", "on", "off"], default="auto")
    p.add_argument("--fe_sg_max", type=float, default=4.0)
    p.add_argument("--no_label_fe", action="store_true",
                   help="exclude the Fe row from `labelled` (see TRAP 4)")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    common = dict(lig_res=args.lig_res, hem_h=args.hem_h,
                  hem_protonation=args.hem_protonation, hem_charge=args.hem_charge,
                  hem_spin=args.hem_spin, cys_thiolate=args.cys_thiolate,
                  fe_sg_max=args.fe_sg_max, label_fe=not args.no_label_fe)
    print(f"POCKET_CUT={POCKET_CUT} A  POCKET_REF={POCKET_REF}  "
          f"POCKET_FORCE_RES={POCKET_FORCE_RES or 'none'}  MAX_ATOMS={MAX_ATOMS}  "
          f"WARN_ATOMS={WARN_ATOMS}  pH={PH}  "
          f"PREP_SEED={'unset (unseeded hydrogens)' if PREP_SEED_VAL is None else PREP_SEED_VAL}"
          f"  code={CODE_VERSION}")

    if args.single:
        if not (args.cif and args.smiles and args.out):
            raise SystemExit("FATAL: --single needs --cif, --smiles and --out")
        if not args.cif.exists():
            raise SystemExit(f"FATAL: --cif {args.cif} does not exist")
        prep_one(str(args.cif), args.smiles, str(args.out), **common)
        return 0

    import pandas as pd

    if args.out_dir is None:
        raise SystemExit("FATAL: --out_dir is required outside --single")
    work = load_work(args)
    mine = work.iloc[args.shard::args.nshards]
    if args.limit:
        mine = mine.head(args.limit)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    print(f"shard {args.shard}/{args.nshards}: {len(mine)} of {len(work)} pairs "
          f"-> {args.out_dir}")
    if len(mine) == 0:
        raise SystemExit(f"FATAL: shard {args.shard} has no work; the work list has "
                         f"{len(work)} rows. Check --nshards.")

    import glob as _glob

    ok = skipped = nocif = failed = 0
    rows, fails = [], []
    for row in mine.itertuples(index=False):
        pid = str(row.pair_id)
        out_npz = args.out_dir / f"{pid}_prep.npz"
        if out_npz.exists() and not args.overwrite:
            skipped += 1
            continue
        cifs = sorted_poses(_glob.glob(str(row.pose_glob)))
        if not cifs:
            nocif += 1                      # not co-folded yet; a later pass picks it up
            continue
        pick = 0
        if args.pose >= 0:
            hit = [i for i, c in enumerate(cifs)
                   if (m := SAMPLE_RE.search(c)) and int(m.group(1)) == args.pose]
            if not hit:
                fails.append({"pair_id": pid, "error": f"sample {args.pose} not among "
                                                       f"{[Path(c).name for c in cifs]}"})
                failed += 1
                continue
            pick = hit[0]
        else:
            # Prefer the PRUNER's verdict (highest ranking_score) over the lowest sample
            # index; they only coincide because PRUNE=1 usually leaves one CIF behind.
            keep = kept_sample(cifs[0])
            hit = [] if keep is None else [
                i for i, c in enumerate(cifs)
                if (m := SAMPLE_RE.search(c)) and int(m.group(1)) == keep]
            if hit:
                pick = hit[0]
            elif keep is not None:
                print(f"  WARNING: {pid}: kept.json names sample {keep} but only "
                      f"{[Path(c).name for c in cifs]} are present; taking "
                      f"{Path(cifs[0]).name}")
            elif len(cifs) > 1:
                print(f"  WARNING: {pid}: {len(cifs)} pose CIFs and no kept.json (no "
                      f"pruner verdict available); taking the lowest index "
                      f"{Path(cifs[0]).name} out of {[Path(c).name for c in cifs]}")
        m = SAMPLE_RE.search(cifs[pick])
        try:
            rows.append(prep_one(
                cifs[pick], str(row.smiles), str(out_npz), pair_id=pid,
                isoform=str(getattr(row, "isoform", "")),
                pose_rank=int(m.group(1)) if m else -1, n_poses=len(cifs),
                has_hem_requested=bool(getattr(row, "has_hem_requested", True)),
                fe_cys_bond_requested=bool(getattr(row, "fe_cys_bond_requested", True)),
                **common))
            ok += 1
        except Exception as exc:  # noqa: BLE001 -- one bad pocket must not cost the shard
            failed += 1
            fails.append({"pair_id": pid, "cif": cifs[pick],
                          "error": f"{type(exc).__name__}: {exc}"})
            print(f"  FAIL {pid}: {type(exc).__name__}: {exc}")
            traceback.print_exc(limit=2)

    tag = f"{args.shard:04d}"
    if rows:
        pd.DataFrame(rows).to_csv(args.out_dir / f"prep_summary_{tag}.csv", index=False)
    if fails:
        pd.DataFrame(fails).to_csv(args.out_dir / f"prep_failures_{tag}.csv", index=False)
    print(f"\nshard {args.shard}: ok={ok} skipped={skipped} no_cif={nocif} failed={failed}")
    if rows:
        n = np.array([r["n_atoms"] for r in rows])
        print(f"  atoms: min {n.min()} median {int(np.median(n))} max {n.max()} "
              f"({int((n > WARN_ATOMS).sum())} above WARN_ATOMS={WARN_ATOMS})")
        nohem = [r["pair_id"] for r in rows if not r["has_hem"]]
        if nohem:
            print(f"  *** {len(nohem)}/{len(rows)} prepped systems have NO HEM ***")
        nocharge = [r["pair_id"] for r in rows if not r["pocket_charge_ok"]]
        if nocharge:
            print(f"  *** {len(nocharge)}/{len(rows)} prepped systems have an UNTRUSTED "
                  f"pocket charge (pocket_charge_ok=False) ***")
    # A shard that prepped nothing must NOT report success: with stage 40 not yet run every
    # pair falls into `nocif`, and an unconditional `return 0` would show 20 green array
    # tasks over an empty directory to anyone chaining --dependency=afterok.
    if ok == 0 and skipped == 0:
        print(f"FATAL: shard {args.shard} produced NO npz "
              f"({nocif} pairs with no pose CIF, {failed} failed).")
        if nocif:
            example = str(mine.iloc[0].pose_glob) if len(mine) else "<no rows>"
            print(f"       The pose CIFs are written by stage 40 "
                  f"(slurm/cyp_challenge/40_cofold_array.sbatch, which runs Protenix-v2) "
                  f"to\n         {example}\n"
                  f"       Run stage 40 for this tag first; stage 50 has nothing to prep.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
