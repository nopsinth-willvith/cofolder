# Vendored from $COFOLD_ROOT/analysis/cyp3a4_gamd/build_system.py -- unchanged apart from the COFOLD_ROOT default and the sibling import of substructure_priors.
"""Stage G1 -- build an AMBER system (prmtop/inpcrd) for ONE CYP3A4 pose: protein + ferric heme
+ blinded ligand, solvated and neutralised, ready for equilibration and GaMD.

FORCE FIELD, and where it departs from the paper. Boulaamane et al. build their heme by hand with
MCPB.py: B3LYP-D3 geometry + Hessian in Gaussian09, RESP charges, Seminario bonded terms. There is
no QM engine on this cluster (no g09/g16/orca/psi4/xtb), so that build cannot be reproduced. We use
the published equivalent that AmberTools already ships:

    $AMBERHOME/dat/contrib/Shahrokh_heme/IC6   Shahrokh, Orendt, Yost & Cheatham, J Comput Chem
                                               (2011) -- QM-derived AMBER-compatible heme
                                               parameters for the P450 catalytic cycle

IC6 is the ferric penta-coordinate high-spin state: the resting enzyme a substrate docks into, and
the same chemistry the paper parameterised. It is a BONDED model -- fe-nc/fe-nd (Fe-porphyrin) and
fe-SH (Fe-thiolate) bond and angle terms, B3LYP/LACVP geometry, RED-III RESP charges -- shipped as
HEM.mol2 plus CYP.mol2, a drop-in replacement residue for the axial cysteinate.

Everything else follows the paper: ff14SB protein, GAFF + AM1-BCC ligand (antechamber/sqm, no QM
code needed), TIP3P in a cubic box with a 10 A solute-wall buffer, Na+/Cl- neutralisation. GAFF1
not GAFF2 is deliberate: the Shahrokh frcmod is written in GAFF1 types and loading leaprc.gaff2
beside it would silently overwrite shared type parameters.

THREE THINGS THIS SCRIPT FIXES IN THE INPUT, because the poses are model output, not structures:

1. N-TERMINAL TRANSMEMBRANE ANCHOR. The poses are full-length 503-residue CYP3A4, so they include
   the TM helix (res ~1-22) that every CYP3A4 crystal construct truncates. Solvated in bulk water
   that helix has nothing to insert into; it collapses onto the surface and dominates any global
   RMSD/RMSF/PCA. --first-residue (default 25) excises it, matching 5TE8/5VC0/9BV5.

2. DISTORTED LIGAND GEOMETRY. The predictor emits C=O at 1.16 A and C#N at 1.07 A -- not physical
   heavy-atom bonds. Taking those into AM1-BCC risks bad charges or an sqm failure. The ligand is
   therefore MMFF-minimised under harmonic POSITION RESTRAINTS (not fixed points) before
   antechamber: bond lengths relax to physical values while every heavy atom stays on its pose.
   The displacement this costs is measured and recorded in build.json.

3. TYPE-II Fe COORDINATION. Despite the 'nontypeii' label, most of these poses put a ligand
   nitrogen 2.0-2.5 A from Fe. The Shahrokh model has no Fe-N(ligand) term, so plain MD sees only
   LJ + Coulomb there and the contact relaxes out to ~3.6 A within picoseconds -- erasing the very
   geometry the pose encodes. Where a STRONG donor (azole or pyridine N, classified by this
   project's own analysis/cyp3a4_heme/substructure_priors.py) sits within --fe-donor-cutoff of Fe,
   an explicit harmonic Fe-N bond is added to the prmtop with ParmEd AFTER tleap.
   Why ParmEd and not tleap: tleap would demand bond AND angle parameters for every fe-N-X triple,
   none of which exist for an exogenous donor, so they would have to be invented. A bare bond term
   invents exactly one number pair (k, r0) and lets sterics set the orientation. It also removes
   the Fe/N pair from the nonbonded list, which is what makes the 2.2 A contact stable at all.
   This ENFORCES the predicted coordination -- MD can no longer falsify it. That is a deliberate
   choice, recorded per pose in build.json as fe_ligand_bond.

Outputs under <out-dir>/<pose-id>/:  system.prmtop  system.inpcrd  complex.pdb  tleap.in/.out
                                     LIG.mol2 LIG.frcmod  build.json
"""
from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

# --------------------------------------------------------------------------- protonation
# PREP_PH -- opt-in pH-aware residue naming, same env name the MLFF prep already uses
# (analysis/cyp_challenge/prep_cofold_system.py:158), so one variable sets both halves of the
# campaign. DEFAULT 7.4 reproduces tleap's historical behaviour EXACTLY: nothing is renamed and
# leaprc.protein.ff14SB:126 maps HIS -> HIE as it always has, so every existing system rebuilds
# byte-identical.
#
# There is no pH anywhere else in this pipeline -- tleap assigns protonation from the residue
# NAME and has no pH concept -- so renaming at write time is the only lever that exists.
PROTONATION_PH = float(os.environ.get("PREP_PH", "7.4"))
HIP_PH_MAX = 6.0          # below this, His is predominantly protonated


def ph_resname(resname: str) -> str:
    """His is the ONLY residue this renames, and that is a deliberate, narrow choice.

    At pH 4.8: His (pKa ~6.5) is ~98% protonated, so HIP is right. But Asp (pKa 3.9) is still
    ~89% ionised and Glu (pKa 4.25) ~78% -- a blanket ASP->ASH / GLU->GLH rename would be WRONG
    for most of them, and this pipeline has no pKa calculator to decide per residue, so they stay
    ionised. Lys, Arg and Tyr are unchanged at both pH.

    Called only for ordinary residues: CYP (Fe-proximal thiolate) and CYX (disulfide) are fixed
    by their bonding, not by pH, and keep precedence at the call site.
    """
    if PROTONATION_PH <= HIP_PH_MAX and resname == "HIS":
        return "HIP"
    return resname


ROOT = Path(os.environ.get("COFOLD_ROOT", "/project/bme-20255165-cofolding-physics"))
AMBERHOME = Path(os.environ.get("AMBERHOME", ROOT / "miniforge3/envs/openmm83"))
HEME_DIR = AMBERHOME / "dat/contrib/Shahrokh_heme/IC6"
sys.path.insert(0, str(Path(__file__).resolve().parent))  # substructure_priors.py (donor classification), vendored alongside

SS_CUTOFF = 2.5        # A, SG-SG: a disulfide is ~2.05 A; the next-nearest pair is far beyond
FE_SG_CUTOFF = 3.2     # A, Fe-SG of the proximal cysteinate (observed 2.23-2.28 A)
PEPTIDE_CUTOFF = 2.0   # A, C(i)-N(i+1): a real peptide bond is 1.33 A
FE_N_K = 100.0         # kcal/mol/A^2, same stiffness the IC6 set uses for fe-nc/fe-nd
FE_N_R0 = 2.15         # A, crystallographic Fe-N(pyridine/azole) in type-II P450 complexes.
#                        A chemical constant, deliberately NOT the per-pose predicted distance:
#                        that would bake prediction noise into the restraint.
RESTRAINT_K = 10.0     # kcal/mol/A^2, MMFF position restraint while ligand bonds are repaired


# --------------------------------------------------------------------------- PDB handling

class Atom:
    __slots__ = ("record", "name", "altloc", "resname", "chain", "resseq", "icode", "xyz",
                 "element")

    def __init__(self, line: str):
        self.record = line[:6].strip()
        self.name = line[12:16].strip()
        self.altloc = line[16]
        self.resname = line[17:20].strip()
        self.chain = line[21]
        self.resseq = int(line[22:26])
        self.icode = line[26]
        self.xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        self.element = (line[76:78].strip() or self.name[0]).capitalize()

    def key(self):
        return (self.chain, self.resseq, self.icode)

    def to_line(self, serial: int, resname: str | None = None) -> str:
        rn = resname or self.resname
        name = self.name if len(self.name) >= 4 else f" {self.name:<3s}"
        return (f"{'ATOM  ' if self.record == 'ATOM' else 'HETATM'}{serial:5d} "
                f"{name:<4s} {rn:>3s} {self.chain}{self.resseq:4d}{self.icode}   "
                f"{self.xyz[0]:8.3f}{self.xyz[1]:8.3f}{self.xyz[2]:8.3f}"
                f"{1.0:6.2f}{0.0:6.2f}          {self.element.upper():>2s}\n")


def dist(a, b) -> float:
    return math.dist(a, b)


def read_pose(path: Path, first_residue: int):
    """Split a pose PDB into protein / heme / ligand, dropping the N-terminal anchor."""
    protein, heme, ligand, dropped_alt, truncated = [], [], [], 0, 0
    for line in path.read_text().splitlines():
        if line[:6] not in ("ATOM  ", "HETATM"):
            continue
        atom = Atom(line)
        if atom.element == "H":              # poses are heavy-atom only; H come from tleap/RDKit
            continue
        if atom.altloc not in (" ", "A"):
            dropped_alt += 1
            continue
        if atom.record == "ATOM":
            if atom.resseq < first_residue:
                truncated += 1
                continue
            protein.append(atom)
        elif atom.resname == "HEM":
            heme.append(atom)
        else:
            ligand.append(atom)
    if not heme:
        raise SystemExit(f"{path.name}: no HEM found")
    if not ligand:
        raise SystemExit(f"{path.name}: no ligand HETATM found")
    return protein, heme, ligand, dropped_alt, truncated


def find_proximal_cys(protein: list[Atom], heme: list[Atom]):
    fe = next((a for a in heme if a.name == "FE"), None)
    if fe is None:
        raise SystemExit("heme has no FE atom")
    best, best_d = None, 1e9
    for a in protein:
        if a.name == "SG" and a.resname in ("CYS", "CYX"):
            d = dist(fe.xyz, a.xyz)
            if d < best_d:
                best, best_d = a, d
    if best is None or best_d > FE_SG_CUTOFF:
        raise SystemExit(f"no cysteine SG within {FE_SG_CUTOFF} A of FE (closest {best_d:.2f} A)")
    return best.key(), best_d, fe.xyz


def find_disulfides(protein: list[Atom], skip_key) -> list[tuple]:
    sgs = [a for a in protein if a.name == "SG" and a.key() != skip_key]
    pairs = []
    for i, a in enumerate(sgs):
        for b in sgs[i + 1:]:
            d = dist(a.xyz, b.xyz)
            if d < SS_CUTOFF:
                pairs.append((a.key(), b.key(), d))
    return pairs


def find_chain_breaks(protein: list[Atom]) -> list[tuple]:
    """Residues after which a TER is needed -- tleap bonds residues in file order, so an unmarked
    gap is welded shut into a stretched peptide bond."""
    order, seen = [], set()
    for a in protein:
        if a.key() not in seen:
            seen.add(a.key())
            order.append(a.key())
    coords = {(a.key(), a.name): a.xyz for a in protein if a.name in ("C", "N")}
    breaks = []
    for cur, nxt in zip(order, order[1:]):
        if cur[0] != nxt[0]:
            breaks.append((cur, nxt, float("nan")))
            continue
        c, n = coords.get((cur, "C")), coords.get((nxt, "N"))
        if c is not None and n is not None and dist(c, n) > PEPTIDE_CUTOFF:
            breaks.append((cur, nxt, dist(c, n)))
    return breaks


# --------------------------------------------------------------------------- pocket truncation

def select_pocket(protein: list[Atom], ligand: list[Atom], free_r: float, buffer_r: float,
                  always_keep: set, gap_fill: int = 2):
    """Residues to keep, split into a mobile zone and a rigid buffer.

    free zone    any heavy atom within free_r of any ligand heavy atom -- these move
    buffer zone  within buffer_r -- kept so the pocket has something to hold its shape against,
                 then made rigid (see freeze_buffer)
    beyond       deleted, and the cut ends capped

    Short gaps are filled in: a single skipped residue between two kept ones would otherwise become
    two extra chain breaks and four extra caps, for no structural gain.
    """
    lig_xyz = [a.xyz for a in ligand]
    by_res: dict = {}
    for a in protein:
        by_res.setdefault(a.key(), []).append(a)
    min_d = {k: min(dist(a.xyz, l) for a in atoms for l in lig_xyz)
             for k, atoms in by_res.items()}

    free = {k for k, d in min_d.items() if d <= free_r} | set(always_keep)
    buffer = {k for k, d in min_d.items() if free_r < d <= buffer_r} - free
    keep = free | buffer

    kept_by_chain: dict = {}
    for chain, resseq, icode in keep:
        kept_by_chain.setdefault(chain, set()).add(resseq)
    for chain, numbers in kept_by_chain.items():
        ordered = sorted(numbers)
        for lo, hi in zip(ordered, ordered[1:]):
            if 1 < hi - lo <= gap_fill + 1:
                for n in range(lo + 1, hi):
                    key = (chain, n, " ")
                    if key in by_res:
                        buffer.add(key)
                        keep.add(key)
    return free, buffer, keep, min_d


def build_fragments(protein: list[Atom], keep: set) -> list[list]:
    """Contiguous runs of kept residues, in file order -- one capped fragment each."""
    order, seen = [], set()
    for a in protein:
        if a.key() not in seen:
            seen.add(a.key())
            order.append(a.key())
    fragments, current = [], []
    for key in order:
        if key not in keep:
            continue
        if current and key[0] == current[-1][0] and key[1] == current[-1][1] + 1:
            current.append(key)
        else:
            if current:
                fragments.append(current)
            current = [key]
    if current:
        fragments.append(current)
    return fragments


def make_cap(flank_atoms: list[Atom], kind: str) -> list[Atom]:
    """Turn the residue flanking a cut into an ACE or NME cap, reusing its real coordinates.

    ACE is CH3-CO-, so the preceding residue's CA/C/O become CH3/C/O.
    NME is -NH-CH3, so the following residue's N/CA become N/C -- note the methyl carbon of NME is
    named "C" in AmberTools' aminoct12.lib, NOT "CH3" as in ACE. Getting this wrong is silent until
    tleap reports "Atom .R<NME n>.A<CH3 7> does not have a type": it adds its own methyl from the
    template and leaves the mis-named one as an extra, untyped atom.
    tleap builds the missing hydrogens. Capping matters: without it every cut end becomes a charged
    NH3+/COO- terminus, and a truncated pocket would acquire a dozen spurious formal charges right
    next to the ligand.
    """
    wanted = {"ACE": {"CA": "CH3", "C": "C", "O": "O"},
              "NME": {"N": "N", "CA": "C"}}[kind]
    out = []
    for atom in flank_atoms:
        if atom.name not in wanted:
            continue
        cap = Atom(atom.to_line(1))
        cap.name = wanted[atom.name]
        cap.resname = kind
        cap.element = atom.element      # CA/C -> C, N -> N, O -> O
        out.append(cap)
    missing = set(wanted.values()) - {a.name for a in out}
    if missing:
        raise SystemExit(f"cannot build {kind} cap, flanking residue lacks {sorted(missing)}")
    return sorted(out, key=lambda a: list(wanted.values()).index(a.name))


def freeze_buffer(work: Path, pdb_residue_order: list, rigid_indices: set,
                  big_mass: float = 1.0e6) -> dict:
    """Make the buffer shell rigid by giving its heavy atoms a huge mass.

    WHY MASS AND NOT A RESTRAINT. gamdRunner builds the System itself, inside gamdSimulation.py,
    straight from the prmtop -- there is no hook to add a CustomExternalForce, and the XML has no
    restraint section. Mass is the one handle that travels in the prmtop.

    WHY A SEPARATE FILE. This is written to system_rigid.prmtop and system.prmtop is left with
    normal masses, because the two stages need different things. Stage G2 equilibrates under a
    barostat, and a barostat rescales MOLECULE CENTROIDS -- a coordinate transformation that mass
    does not resist. Truncation turns the protein into ~20 separate capped fragments, so running
    NPT against frozen masses would let the box scaling walk the fragments around each other while
    nothing pulled them back. G2 therefore uses normal masses with explicit position restraints on
    this same residue set, and only GaMD (NVT, no barostat) uses the rigid topology.

    WHY NOT ZERO MASS. The fork hardcodes createSystem(hydrogenMass=1.5*amu), which moves 0.49 amu
    off each heavy atom onto its hydrogens. From a zero mass that gives a NEGATIVE mass and the run
    dies. From 1e6 it is a rounding error.

    At 1e6 amu Langevin still thermostats these atoms correctly (equipartition is mass-independent,
    so the temperature stays honest) but positional diffusion goes as 1/m, i.e. ~1e5 times slower
    than carbon: rigid on any timescale here. Their hydrogens are still reset to 1.5 amu by the
    hardcoded HMR and can rotate about the frozen heavy scaffold, which is harmless in the buffer.
    """
    import parmed

    parm = parmed.load_file(str(work / "system.prmtop"), str(work / "system.inpcrd"))
    # tleap keeps input residue order and appends solvent, so the first N prmtop residues map 1:1
    # onto the complex PDB. That is asserted, not assumed.
    aliases = {"HIS": {"HID", "HIE", "HIP"}, "CYS": {"CYX", "CYM"}, "ASP": {"ASH"},
               "GLU": {"GLH"}, "LYS": {"LYN"}}
    if len(parm.residues) < len(pdb_residue_order):
        raise SystemExit(f"prmtop has {len(parm.residues)} residues, fewer than the "
                         f"{len(pdb_residue_order)} written")
    for i, (name, _) in enumerate(pdb_residue_order):
        got = parm.residues[i].name
        if got != name and got not in aliases.get(name, set()):
            raise SystemExit(f"prmtop residue {i} is {got}, expected {name} -- residue order is "
                             f"not what freezing assumes")
    n_frozen = 0
    for idx in sorted(rigid_indices):
        for atom in parm.residues[idx].atoms:
            if atom.atomic_number > 1:
                atom.mass = big_mass
                n_frozen += 1
    parm.save(str(work / "system_rigid.prmtop"), overwrite=True)
    return {"n_rigid_residues": len(rigid_indices), "n_rigid_heavy_atoms": n_frozen,
            "rigid_mass_amu": big_mass,
            "rigid_residue_indices": sorted(i + 1 for i in rigid_indices)}


# --------------------------------------------------------------------------- ligand

def build_ligand(smiles: str, ligand: list[Atom], workdir: Path, name: str = "LIG"):
    """SMILES + pose heavy-atom coordinates -> repaired geometry, GAFF/AM1-BCC mol2 + frcmod."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Geometry import Point3D

    ref = Chem.MolFromSmiles(smiles)
    if ref is None:
        raise SystemExit(f"RDKit could not parse SMILES: {smiles}")
    frags = Chem.GetMolFrags(ref, asMols=True, sanitizeFrags=True)
    if len(frags) > 1:                        # salts: keep the parent
        ref = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    if ref.GetNumAtoms() != len(ligand):
        raise SystemExit(f"ligand atom-count mismatch: SMILES {ref.GetNumAtoms()} heavy atoms vs "
                         f"PDB {len(ligand)}")

    # The pose writer emits heavy atoms in SMILES order, but that is CHECKED, not assumed: a silent
    # mis-mapping would parameterise the right molecule onto the wrong coordinates.
    if [a.GetSymbol() for a in ref.GetAtoms()] == [a.element for a in ligand]:
        order, mapping = list(range(len(ligand))), "positional"
    else:
        order, mapping = _match_by_graph(ref, ligand), "graph-match"

    conf = Chem.Conformer(ref.GetNumAtoms())
    for ref_idx, pdb_idx in enumerate(order):
        x, y, z = ligand[pdb_idx].xyz
        conf.SetAtomPosition(ref_idx, Point3D(x, y, z))
    ref.RemoveAllConformers()
    ref.AddConformer(conf, assignId=True)
    Chem.AssignStereochemistryFrom3D(ref)
    mol_h = Chem.AddHs(ref, addCoords=True)
    max_disp, rms_disp = _repair_geometry(mol_h)

    net_charge = Chem.GetFormalCharge(mol_h)
    sdf = workdir / "ligand_repaired.sdf"
    with Chem.SDWriter(str(sdf)) as w:
        w.write(mol_h)

    mol2 = workdir / f"{name}.mol2"
    run(["antechamber", "-i", str(sdf), "-fi", "sdf", "-o", str(mol2), "-fo", "mol2",
         "-c", "bcc", "-nc", str(net_charge), "-at", "gaff", "-rn", name, "-pf", "n", "-dr", "n"],
        cwd=workdir, what="antechamber")
    if not mol2.is_file():
        raise SystemExit(f"antechamber produced no mol2 -- see {workdir}/sqm.out")

    frcmod = workdir / f"{name}.frcmod"
    run(["parmchk2", "-i", str(mol2), "-f", "mol2", "-o", str(frcmod), "-s", "gaff"],
        cwd=workdir, what="parmchk2")
    lig_pdb = workdir / f"{name}.pdb"
    run(["antechamber", "-i", str(mol2), "-fi", "mol2", "-o", str(lig_pdb), "-fo", "pdb",
         "-pf", "y", "-dr", "n"], cwd=workdir, what="antechamber_pdb")

    donors = _strong_donors(mol_h, order, ligand)
    return {
        "mol2": mol2, "frcmod": frcmod, "pdb": lig_pdb, "net_charge": net_charge,
        "mapping": mapping, "n_atoms": mol_h.GetNumAtoms(), "donors": donors,
        "repair_max_disp": round(max_disp, 3), "repair_rms_disp": round(rms_disp, 3),
    }


def _match_by_graph(ref, ligand):
    from rdkit import Chem

    block = "".join(a.to_line(i + 1) for i, a in enumerate(ligand)) + "END\n"
    probe = Chem.MolFromPDBBlock(block, removeHs=False, sanitize=False)
    if probe is None:
        raise SystemExit("could not read ligand PDB block for graph matching")
    match = probe.GetSubstructMatch(Chem.MolFromSmarts(Chem.MolToSmarts(ref)))
    if len(match) != ref.GetNumAtoms():
        raise SystemExit("ligand graph match failed -- SMILES does not fit the pose connectivity")
    return list(match)


def _repair_geometry(mol_h):
    """Relax unphysical predicted bond lengths without re-docking: MMFF under harmonic position
    restraints on every heavy atom. Returns (max, rms) heavy-atom displacement in A."""
    from rdkit.Chem import AllChem

    before = mol_h.GetConformer().GetPositions().copy()
    props = AllChem.MMFFGetMoleculeProperties(mol_h)
    if props is None:                          # no MMFF types: leave the geometry alone
        return 0.0, 0.0
    ff = AllChem.MMFFGetMoleculeForceField(mol_h, props)
    for atom in mol_h.GetAtoms():
        if atom.GetAtomicNum() > 1:
            ff.MMFFAddPositionConstraint(atom.GetIdx(), 0.0, RESTRAINT_K)
    ff.Minimize(maxIts=2000)
    after = mol_h.GetConformer().GetPositions()
    heavy = [a.GetIdx() for a in mol_h.GetAtoms() if a.GetAtomicNum() > 1]
    d = [math.dist(before[i], after[i]) for i in heavy]
    return max(d), math.sqrt(sum(x * x for x in d) / len(d))


def _strong_donors(mol_h, order, ligand) -> list[dict]:
    """Strong Type-II donors (azole / pyridine N), mapped back to pose atom names.

    Classification is this project's own, from analysis/cyp3a4_heme/substructure_priors.py, so the
    MD setup and the earlier heme analyses agree on what counts as a donor.
    """
    from rdkit import Chem
    from substructure_priors import classify_heteroatoms

    heavy_only = Chem.RemoveHs(mol_h)
    rows = classify_heteroatoms(heavy_only)
    donors = []
    for row in rows:
        if row["donor"] != "strong":
            continue
        pdb_idx = order[row["idx"]]
        donors.append({"rdkit_idx": row["idx"], "pose_atom": ligand[pdb_idx].name,
                       "class": row["class"], "element": row["element"]})
    return donors


# --------------------------------------------------------------------------- shell helpers

def run(cmd, cwd: Path, what: str, check: bool = True):
    proc = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True, timeout=7200)
    (cwd / f"{what}.log").write_text(proc.stdout + "\n--- stderr ---\n" + proc.stderr)
    if check and proc.returncode != 0:
        sys.stderr.write(proc.stdout[-3000:] + "\n" + proc.stderr[-3000:] + "\n")
        raise SystemExit(f"{what} failed (rc={proc.returncode}) in {cwd}")
    return proc


TLEAP_TEMPLATE = """\
source leaprc.protein.ff14SB
source leaprc.water.tip3p
source leaprc.gaff

# --- ferric penta-coordinate high-spin heme (Shahrokh et al. 2011, AmberTools IC6 set) ---
loadamberparams {heme_dir}/IC6.frcmod
HEM = loadmol2 {heme_dir}/HEM.mol2
CYP = loadmol2 {heme_dir}/CYP.mol2
# CYP replaces the axial cysteine, so it needs backbone connectivity or it stays a floating unit
set CYP.1 restype protein
set CYP.1 connect0 CYP.1.N
set CYP.1 connect1 CYP.1.C
set CYP head CYP.1.N
set CYP tail CYP.1.C

# --- ligand (GAFF / AM1-BCC) ---
loadamberparams {lig_frcmod}
LIG = loadmol2 {lig_mol2}

mol = loadpdb {complex_pdb}
{bonds}
check mol

solvatebox mol TIP3PBOX 10.0
addions mol Na+ 0
addions mol Cl- 0

saveamberparm mol {prmtop} {inpcrd}
quit
"""


# --------------------------------------------------------------------------- driver

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pose", type=Path, required=True, help="pose PDB (protein + HEM + LIG)")
    ap.add_argument("--smiles", required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--pose-id", default=None, help="defaults to the pose file stem")
    ap.add_argument("--first-residue", type=int, default=25,
                    help="drop protein residues below this (25 excises the TM anchor)")
    ap.add_argument("--fe-donor-cutoff", type=float, default=2.6,
                    help="A; a strong donor nearer than this to Fe gets an explicit Fe-N bond")
    ap.add_argument("--no-fe-bond", action="store_true",
                    help="never add the Fe-ligand bond (the uniform/unrestrained control arm)")
    ap.add_argument("--pocket-free-radius", type=float, default=0.0,
                    help="A around the ligand that stays mobile. 0 = keep the whole protein")
    ap.add_argument("--pocket-buffer-radius", type=float, default=16.0,
                    help="A around the ligand kept as a RIGID buffer; beyond this is deleted")
    args = ap.parse_args(argv)
    pocket_mode = args.pocket_free_radius > 0

    pose_id = args.pose_id or args.pose.stem
    work = args.out_dir / pose_id
    work.mkdir(parents=True, exist_ok=True)
    for b in ("antechamber", "parmchk2", "tleap"):
        if shutil.which(b) is None:
            raise SystemExit(f"{b} not on PATH (AMBERHOME={AMBERHOME})")
    if not (HEME_DIR / "IC6.frcmod").is_file():
        raise SystemExit(f"Shahrokh IC6 heme parameters not found at {HEME_DIR}")

    protein, heme, ligand, dropped_alt, truncated = read_pose(args.pose, args.first_residue)
    prox_key, fe_sg, fe_xyz = find_proximal_cys(protein, heme)
    ss = find_disulfides(protein, prox_key)
    breaks = find_chain_breaks(protein)
    ss_keys = {k for pair in ss for k in pair[:2]}
    first_res = min(a.resseq for a in protein)
    last_res = max(a.resseq for a in protein)

    print(f"[{pose_id}] protein res {first_res}-{last_res} ({len(protein)} heavy atoms; "
          f"dropped {truncated} atoms below residue {args.first_residue})")
    print(f"[{pose_id}] proximal cysteine {prox_key[0]}{prox_key[1]} -> CYP, Fe-SG {fe_sg:.2f} A")
    print(f"[{pose_id}] disulfides {[(a[1], b[1]) for a, b, _ in ss] or 'none'}; "
          f"chain breaks {[(a[1], b[1]) for a, b, _ in breaks] or 'none'}")

    lig = build_ligand(args.smiles, ligand, work)
    print(f"[{pose_id}] ligand {lig['n_atoms']} atoms with H, formal charge {lig['net_charge']:+d}, "
          f"mapping {lig['mapping']}, geometry repair moved heavy atoms "
          f"max {lig['repair_max_disp']} A / rms {lig['repair_rms_disp']} A")

    # --- which strong donor, if any, is coordinating Fe ---
    by_name = {a.name: a for a in ligand}
    coordinating = None
    for d in lig["donors"]:
        atom = by_name.get(d["pose_atom"])
        if atom is None:
            continue
        r = dist(fe_xyz, atom.xyz)
        if r <= args.fe_donor_cutoff and (coordinating is None or r < coordinating["distance"]):
            coordinating = {**d, "distance": round(r, 3)}
    if coordinating and not args.no_fe_bond:
        print(f"[{pose_id}] type-II coordination: {coordinating['class']} "
              f"{coordinating['pose_atom']} at {coordinating['distance']} A -> adding Fe-N bond")
    elif coordinating:
        print(f"[{pose_id}] type-II coordination present ({coordinating['distance']} A) but "
              f"--no-fe-bond given: left unbonded")
    else:
        print(f"[{pose_id}] no strong donor within {args.fe_donor_cutoff} A of Fe")

    # --- select the pocket, if truncating ---
    by_res: dict = {}
    for atom in protein:
        by_res.setdefault(atom.key(), []).append(atom)
    if pocket_mode:
        free, buffer, keep, min_d = select_pocket(
            protein, ligand, args.pocket_free_radius, args.pocket_buffer_radius,
            always_keep={prox_key})
        ss = [(a, b, d) for a, b, d in ss if a in keep and b in keep]
        ss_keys = {k for pair in ss for k in pair[:2]}
        print(f"[{pose_id}] pocket: {len(free)} mobile + {len(buffer)} rigid residues kept of "
              f"{len(by_res)} (free {args.pocket_free_radius} A, buffer "
              f"{args.pocket_buffer_radius} A); Cys{prox_key[1]} is "
              f"{'mobile' if prox_key in free else 'rigid'} at {min_d[prox_key]:.1f} A")
    else:
        free, buffer, keep = set(by_res), set(), set(by_res)

    # --- assemble the complex PDB: protein (CYP/CYX renamed, capped) + heme + repaired ligand ---
    # Residues are renumbered 1..N here, deliberately. tleap does NOT honour the numbers in the
    # file: it numbers residues sequentially starting from the FIRST residue's number and ignores
    # gaps (confirmed with `desc` -- see slurm/cyp3a4_gamd/diag_tleap_index.sbatch, where the field
    # tleap prints as "PDB sequence number" is really the positional index). With a contiguous
    # chain the two coincide, which is why the full-system build worked; the moment truncation
    # punches gaps in the numbering they diverge and `bond mol.442.SG` selects the wrong residue.
    # Renumbering 1..N makes file number, positional index and tleap's number the same thing.
    # The original labels are kept in build.json's residue_map, because an Amber prmtop does not
    # store them and "PHE213" is how anyone reads a contact list.
    complex_pdb = work / "complex.pdb"
    serial, lines = 1, []
    residue_order: list = []            # (resname, original label) in file order
    residue_map: list = []
    rigid_indices: set = set()
    break_after = {a for a, _, _ in breaks}
    fragments = build_fragments(protein, keep)
    def emit(atoms, resname, orig_label, rigid):
        """Write one residue under its new sequential number, recording where it came from."""
        nonlocal serial
        index = len(residue_order) + 1
        if rigid:
            rigid_indices.add(index - 1)
        residue_order.append((resname, orig_label))
        residue_map.append({"index": index, "resname": resname, "original": orig_label,
                            "zone": "rigid" if rigid else "mobile"})
        for atom in atoms:
            if resname == "CYP" and atom.name == "HG":   # a thiolate carries no SH proton
                continue
            atom.resseq = index
            lines.append(atom.to_line(serial, resname))
            serial += 1
        return index

    def emit_cofactors():
        """HEM and LIG go right after the fragment carrying the proximal cysteine.

        They are one molecule with it -- Fe is bonded to that residue's SG, and to the ligand N
        where there is type-II coordination. Writing them at the end of the file instead would
        leave that molecule's atoms split around the other ~20 fragments, and ParmEd then warns
        "Molecule atoms are not contiguous! Attempting to reorder to fix" and renumbers the
        topology out from under the coordinates and the residue map.
        """
        hem = emit(heme, "HEM", "HEM", False)
        lines.append("TER\n")
        lig_atoms = [Atom(l) for l in lig["pdb"].read_text().splitlines()
                     if l[:6] in ("ATOM  ", "HETATM")]
        for a in lig_atoms:
            a.chain, a.record = "C", "HETATM"
        li = emit(lig_atoms, "LIG", "LIG", False)
        lines.append("TER\n")
        return hem, li

    prox_fragment = next((i for i, f in enumerate(fragments) if prox_key in f), None)
    if prox_fragment is None:
        raise SystemExit("the proximal cysteine is not in any kept fragment")
    hem_resseq = lig_resseq = None
    for frag_i, fragment in enumerate(fragments):
        chain = fragment[0][0]
        first, last = fragment[0][1], fragment[-1][1]
        prev_key = (chain, first - 1, " ")
        if pocket_mode and prev_key in by_res:          # cut end -> ACE cap
            emit(make_cap(by_res[prev_key], "ACE"), "ACE",
                 f"{by_res[prev_key][0].resname}{first - 1}", True)
        prev_written = None
        for key in fragment:
            if prev_written is not None and prev_written in break_after:
                lines.append("TER\n")
            resname = by_res[key][0].resname
            if key == prox_key:
                resname = "CYP"
            elif key in ss_keys:
                resname = "CYX"
            else:
                resname = ph_resname(resname)
            emit(by_res[key], resname, f"{by_res[key][0].resname}{key[1]}", key in buffer)
            prev_written = key
        next_key = (chain, last + 1, " ")
        if pocket_mode and next_key in by_res:          # cut end -> NME cap
            emit(make_cap(by_res[next_key], "NME"), "NME",
                 f"{by_res[next_key][0].resname}{last + 1}", True)
        lines.append("TER\n")
        if frag_i == prox_fragment:
            hem_resseq, lig_resseq = emit_cofactors()
    # HEM and LIG both arrive as residue 1 of their own chain. tleap ignores chain IDs, sees
    # residue 1 twice in a row under two names, declares it a "split residue" and renumbers the
    # whole unit -- after which mol.<N>.<atom> no longer refers to what the PDB says it does.
    # Numbering them on from the protein avoids that entirely.
    if hem_resseq is None:                      # whole-protein build: no fragmentation
        hem_resseq, lig_resseq = emit_cofactors()
    lines.append("END\n")
    complex_pdb.write_text("".join(lines))

    # --- tleap ---
    hem_num = _residue_number(complex_pdb, "HEM")
    cyp_num = _residue_number(complex_pdb, "CYP")
    bonds = [f"bond mol.{cyp_num}.SG mol.{hem_num}.FE"]
    for a, b, _ in ss:
        bonds.append(f"bond mol.{_residue_number(complex_pdb, 'CYX', a[1])}.SG "
                     f"mol.{_residue_number(complex_pdb, 'CYX', b[1])}.SG")
    tleap_in = work / "tleap.in"
    tleap_in.write_text(TLEAP_TEMPLATE.format(
        heme_dir=HEME_DIR, lig_frcmod=lig["frcmod"], lig_mol2=lig["mol2"],
        complex_pdb=complex_pdb, bonds="\n".join(bonds),
        prmtop=work / "system.prmtop", inpcrd=work / "system.inpcrd"))
    proc = run(["tleap", "-f", str(tleap_in)], cwd=work, what="tleap", check=False)
    (work / "tleap.out").write_text(proc.stdout)
    if not (work / "system.prmtop").is_file():
        sys.stderr.write(proc.stdout[-5000:])
        raise SystemExit(f"[{pose_id}] tleap wrote no prmtop -- see {work/'tleap.out'}")

    fe_bond = None
    if coordinating and not args.no_fe_bond:
        fe_bond = add_fe_ligand_bond(work, coordinating["pose_atom"])

    rigid = None
    if pocket_mode:
        rigid = freeze_buffer(work, residue_order, rigid_indices)
        rigid.update({"free_radius_A": args.pocket_free_radius,
                      "buffer_radius_A": args.pocket_buffer_radius,
                      "n_mobile_residues": len(free), "n_fragments": len(fragments),
                      "n_caps": sum(1 for name, _ in residue_order if name in ("ACE", "NME"))})
        print(f"[{pose_id}] {len(fragments)} capped fragments, {rigid['n_caps']} caps, "
              f"{rigid['n_rigid_heavy_atoms']} heavy atoms made rigid")

    summary = validate(work, pose_id, fe_sg, lig, ss, breaks,
                       (first_res, last_res), coordinating, fe_bond, rigid)
    summary["residue_map"] = residue_map
    (work / "build.json").write_text(json.dumps(summary, indent=2))
    print(f"[{pose_id}] built: {summary['n_atoms']} atoms, {summary['n_waters']} waters, "
          f"net charge {summary['net_charge']:+.4f}, box {summary['box']}")
    return 0


def _residue_number(pdb: Path, resname: str, resseq: int | None = None) -> int:
    """The number tleap's `mol.<N>.<atom>` selector expects.

    Verified against `desc` (slurm/cyp3a4_gamd/diag_tleap_index.sbatch): tleap's RESIDUE SEQUENCE
    NUMBER is the number written in the PDB, and the field it confusingly labels "PDB sequence
    number" is the positional index. Addressing residues by position therefore selects the wrong
    one -- with this construct `mol.418` is an ARG, not the CYP at 442. Residue numbers are unique
    across the complex by construction: protein 25-503, then HEM and LIG numbered on from there.
    """
    for line in pdb.read_text().splitlines():
        if line[:6] not in ("ATOM  ", "HETATM"):
            continue
        if line[17:20].strip() == resname and (resseq is None or int(line[22:26]) == resseq):
            return int(line[22:26])
    raise SystemExit(f"residue {resname} {resseq or ''} not found in {pdb}")


def add_fe_ligand_bond(work: Path, donor_atom_name: str) -> dict:
    """Add the harmonic Fe-N(ligand) term to the finished prmtop with ParmEd.

    tleap cannot do this: it would require bond AND angle parameters for every fe-N-X triple, none
    of which exist for an exogenous donor. A bare bond invents one (k, r0) pair, and -- because
    remake_parm rebuilds the exclusion list -- it also removes the Fe/N nonbonded pair, which is
    what makes a 2.2 A contact stable instead of a ~14 kcal/mol LJ wall.
    """
    import parmed

    prmtop = work / "system.prmtop"
    parm = parmed.load_file(str(prmtop), str(work / "system.inpcrd"))
    fe = [a for a in parm.atoms if a.name == "FE"]
    lig_res = [r for r in parm.residues if r.name == "LIG"]
    if len(fe) != 1 or len(lig_res) != 1:
        raise SystemExit(f"expected one FE and one LIG, found {len(fe)} / {len(lig_res)}")
    donor = next((a for a in lig_res[0].atoms if a.name == donor_atom_name), None)
    if donor is None:
        names = sorted(a.name for a in lig_res[0].atoms)
        raise SystemExit(f"donor atom {donor_atom_name} not in LIG; have {names}")

    observed = math.dist((fe[0].xx, fe[0].xy, fe[0].xz), (donor.xx, donor.xy, donor.xz))
    btype = parmed.BondType(FE_N_K, FE_N_R0, list=parm.bond_types)
    parm.bond_types.append(btype)
    parm.bonds.append(parmed.Bond(fe[0], donor, type=btype))
    parm.remake_parm()
    # Both files are rewritten from the SAME object. remake_parm can reorder atoms, and a prmtop
    # saved against a stale inpcrd would silently mis-assign every coordinate in the system.
    parm.save(str(prmtop), overwrite=True)
    parm.save(str(work / "system.inpcrd"), format="rst7", overwrite=True)
    return {"donor_atom": donor_atom_name, "k_kcal_per_mol_A2": FE_N_K, "r0_A": FE_N_R0,
            "observed_distance_A": round(observed, 3)}


def validate(work: Path, pose_id: str, fe_sg: float, lig: dict, ss, breaks, res_range,
             coordinating, fe_bond, rigid=None) -> dict:
    """Re-open the prmtop and prove the heme is really bonded -- not just that tleap exited 0."""
    import parmed

    parm = parmed.load_file(str(work / "system.prmtop"), str(work / "system.inpcrd"))
    fe = [a for a in parm.atoms if a.name == "FE"]
    if len(fe) != 1:
        raise SystemExit(f"[{pose_id}] expected exactly one FE, found {len(fe)}")
    partners = [(b.atom1 if b.atom2 is fe[0] else b.atom2) for b in parm.bonds
                if fe[0] in (b.atom1, b.atom2)]
    names = sorted(f"{a.residue.name}:{a.name}" for a in partners)
    n_porph = sum(1 for a in partners if a.residue.name == "HEM" and a.name.startswith("N"))
    n_thiolate = sum(1 for a in partners if a.name == "SG")
    if n_porph != 4 or n_thiolate != 1:
        raise SystemExit(f"[{pose_id}] Fe coordination wrong: bonded to {names} "
                         f"(want 4 porphyrin N + 1 thiolate S)")
    if fe_bond and not any(a.residue.name == "LIG" for a in partners):
        raise SystemExit(f"[{pose_id}] Fe-ligand bond was requested but is absent from the prmtop")
    charge = sum(a.charge for a in parm.atoms)
    # Tolerance is 0.01 e, not machine epsilon: a prmtop stores charges at limited precision, so
    # summing ~40,000 of them leaves a few 1e-3 e of rounding. A genuinely missing counter-ion
    # would show up as +-1 e, which this still catches by two orders of magnitude.
    if abs(charge) > 1e-2:
        raise SystemExit(f"[{pose_id}] system is not neutral: {charge:+.4f} e")
    box = [round(v, 2) for v in (parm.box[:3] if parm.box is not None else [])]
    return {
        "pose_id": pose_id,
        "n_atoms": len(parm.atoms), "n_residues": len(parm.residues),
        "n_waters": sum(1 for r in parm.residues if r.name in ("WAT", "HOH")),
        "net_charge": round(charge, 6), "box": box,
        "protein_residue_range": list(res_range),
        "fe_bonded_to": names, "fe_sg_input_distance": round(fe_sg, 3),
        "fe_ligand_bond": fe_bond,
        "pocket_truncation": rigid,
        "coordinating_donor": coordinating,
        "ligand_formal_charge": lig["net_charge"],
        "ligand_coordinate_mapping": lig["mapping"],
        "ligand_geometry_repair": {"max_disp_A": lig["repair_max_disp"],
                                   "rms_disp_A": lig["repair_rms_disp"]},
        "ligand_strong_donors": lig["donors"],
        "n_disulfides": len(ss), "n_chain_breaks": len(breaks),
        "heme_model": "Shahrokh IC6 ferric penta-coordinate high-spin",
        "protein_ff": "ff14SB", "ligand_ff": "GAFF/AM1-BCC", "water": "TIP3P, 10 A buffer",
        # Recorded so a finished system is not indistinguishable from a pH-7.4 one. Without this
        # the only way to tell them apart is dumping RESIDUE_LABEL from the prmtop.
        "protonation_ph": PROTONATION_PH,
        "protonation_renames": ({"HIS": "HIP"} if PROTONATION_PH <= HIP_PH_MAX else {}),
    }


if __name__ == "__main__":
    raise SystemExit(main())
