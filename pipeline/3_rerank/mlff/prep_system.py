# Vendored from $COFOLD_ROOT/analysis/energy_labels/prep_system.py -- unchanged.
"""Stage 1 of the per-atom interaction-energy label pipeline: turn a Protenix bioassembly
into the three protonated subsystems Orb needs, carrying Protenix join keys through.

For one (pdb_id, assembly_id, ligand label_asym_id, sym_id) system:
  1. load mmcif_bioassembly/<pdb_id>.pkl.gz -> biotite AtomArray (heavy atoms only, no H)
  2. ligand   = the indexed ligand chain (label_asym_id + sym_id_int)
  3. pocket   = WHOLE protein residues with any heavy atom within POCKET_CUT of a
                ligand heavy atom
  4. protonate pocket (PDBFixer @ pH 7.4 + ff14SB H-only relaxation, heavy frozen) and
     ligand (RDKit from the Protenix BondList + MMFF H-only relaxation, heavy frozen)
  5. emit ONE npz holding both subsystems at IDENTICAL geometry, plus join keys.

Join keys are the crux: the training loop must map labels onto Protenix atom ordering, so
every labellable atom carries (label_asym_id, sym_id_int, label_seq_id, auth_seq_id,
res_name, atom_name). Positional indices are deliberately NOT used as keys.

PDBFixer cannot round-trip multi-character label_asym_id (PDB chain is 1 char), so the
pocket is written with temporary single-char chain ids and sequential residue numbers, and
mapped back by (chain order, residue order, atom name) -- which PDBFixer preserves because
we forbid it from modelling missing residues.

Hydrogens are added by us and do not exist in Protenix, so they are never labelled
directly; instead each H is tagged with its parent heavy atom so stage 2 can fold its
energy into that heavy atom (united-atom aggregation).

Usage: prep_system.py <pdb_id> <assembly_id> <lig_asym> <lig_sym> <out_npz>
Run in the openmm83 env (biotite + rdkit + openmm + pdbfixer).
"""
import gzip
import io
import os
import pickle
import sys
import tempfile

import numpy as np
from openmm import (Context, LangevinIntegrator, LocalEnergyMinimizer,
                    NonbondedForce, Platform, unit)
from openmm.app import ForceField, PDBFile
from pdbfixer import PDBFixer
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

RDLogger.DisableLog("rdApp.*")

BIO_DIR = os.environ.get(
    "BIO_DIR", "/scratch-shared/20255165/protenix_data/mmcif_bioassembly")
POCKET_CUT = float(os.environ.get("POCKET_CUT", "10.0"))
MAX_ATOMS = int(os.environ.get("MAX_ATOMS", "6000"))
FF = ForceField("amber14-all.xml")

STD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU",
    "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "MSE", "SEC", "PYL",
}
# biotite BondType -> RDKit bond order. Biotite already stores the KEKULE structure of
# aromatic rings (AROMATIC_SINGLE=5 / AROMATIC_DOUBLE=6 / AROMATIC_TRIPLE=7), so those map
# to plain single/double/triple and RDKit perceives aromaticity during sanitisation.
# Handing RDKit BondType.AROMATIC instead makes Kekulize fail ("Can't kekulize mol").
BOND_ORDER = {
    0: Chem.BondType.SINGLE, 1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE,
    3: Chem.BondType.TRIPLE, 4: Chem.BondType.QUADRUPLE, 5: Chem.BondType.SINGLE,
    6: Chem.BondType.DOUBLE, 7: Chem.BondType.TRIPLE, 8: Chem.BondType.AROMATIC,
}
# single-char chain ids for the PDB round-trip
TMP_CHAINS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"


class _Stub:
    """Placeholder for Protenix classes we do not need (e.g. TokenArray)."""

    def __init__(self, *a, **k):
        pass

    def __setstate__(self, state):
        self.__dict__.update(state if isinstance(state, dict) else {"_state": state})


class _Unpickler(pickle.Unpickler):
    """Loads the bioassembly without requiring protenix on the path (we only want
    the biotite AtomArray, which biotite itself can reconstruct)."""

    def find_class(self, module, name):
        try:
            return super().find_class(module, name)
        except Exception:  # noqa: BLE001
            t = type(f"{module}.{name}", (_Stub,), {})
            t.__module__ = module
            return t


def load_assembly(pdb_id):
    with gzip.open(os.path.join(BIO_DIR, f"{pdb_id}.pkl.gz"), "rb") as fh:
        return _Unpickler(io.BytesIO(fh.read())).load()


def _as_int(v, default=-1):
    """Parse an mmCIF numeric field that may be the null token '.' or '?'.

    label_seq_id is null for every non-polymer atom, so ligand atoms legitimately have no
    label_seq_id and must be joined on auth_seq_id instead -- which is why both are stored.
    """
    try:
        return int(str(v).strip())
    except (TypeError, ValueError):
        return default


def _spin(atomic_numbers_sum, charge):
    return 1 if (atomic_numbers_sum - int(charge)) % 2 == 0 else 2


def _minimize_frozen_heavy(topology, positions):
    """ff14SB local minimisation with every heavy atom frozen (mass 0): relaxes only the
    hydrogens we just added, so the crystal heavy-atom geometry is preserved exactly.
    Returns (coords in A, net integer charge from the ff partial charges)."""
    system = FF.createSystem(topology, rigidWater=False)
    nbf = next(f for f in system.getForces() if isinstance(f, NonbondedForce))
    charge = int(round(sum(nbf.getParticleParameters(i)[0]._value
                           for i in range(system.getNumParticles()))))
    for i, atom in enumerate(topology.atoms()):
        if atom.element is not None and atom.element.symbol != "H":
            system.setParticleMass(i, 0.0)
    integ = LangevinIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond,
                               1.0 * unit.femtosecond)
    ctx = Context(system, integ, Platform.getPlatformByName("CPU"))
    ctx.setPositions(positions)
    LocalEnergyMinimizer.minimize(ctx, maxIterations=500)
    state = ctx.getState(getPositions=True)
    return state.getPositions(asNumpy=True).value_in_unit(unit.angstrom), charge


def select_system(d, lig_asym, lig_sym):
    """Return (ligand_idx, pocket_residue_groups) as index arrays into the AtomArray."""
    aa = d["atom_array"]
    an = aa._annot
    coord = aa._coord
    # Protenix marks unresolved/zero-occupancy atoms; their coordinates are meaningless.
    ok = an["is_resolved"].astype(bool) & (an["occupancy"] > 0.0)

    lig = ok & (an["label_asym_id"] == lig_asym) & (an["sym_id_int"].astype(int) == lig_sym)
    lig_idx = np.flatnonzero(lig)
    if lig_idx.size == 0:
        raise ValueError(f"ligand chain {lig_asym}{lig_sym} has no resolved atoms")
    if lig_idx.size < 2:
        raise ValueError("single-atom ligand (ion) -- not a co-folding ligand")

    prot = ok & an["is_protein"].astype(bool) & np.isin(an["res_name"], list(STD_AA))
    prot_idx = np.flatnonzero(prot)
    if prot_idx.size == 0:
        raise ValueError("no resolved standard-AA protein atoms")

    d2 = ((coord[prot_idx][:, None, :] - coord[lig_idx][None, :, :]) ** 2).sum(-1)
    near = prot_idx[d2.min(axis=1) <= POCKET_CUT ** 2]
    if near.size == 0:
        raise ValueError("no protein atoms within POCKET_CUT of the ligand")

    # expand to WHOLE residues (residue identity = chain + sym copy + seq id + ins code)
    def res_key(i):
        return (str(an["label_asym_id"][i]), int(an["sym_id_int"][i]),
                _as_int(an["label_seq_id"][i]), _as_int(an["auth_seq_id"][i]),
                str(an["ins_code"][i]))

    keep_keys = {res_key(i) for i in near}
    groups = {}
    for i in prot_idx:
        k = res_key(i)
        if k in keep_keys:
            groups.setdefault(k, []).append(i)

    # Per-residue min distance to the ligand, computed HERE on the ORIGINAL AtomArray.
    #
    # This exists so a downstream buffer crop (extract_interface_latents.py --label-cut)
    # can reproduce THIS residue selection at a SMALLER cutoff exactly, instead of
    # re-deriving it from the protonated output. Re-deriving is wrong: PDBFixer's
    # replaceNonstandardResidues() DELETES atoms it cannot map (MSE -> MET drops SE and
    # models a fresh SD, which carries labelled=False), so an MSE residue whose only
    # close contact is via SE is kept by a native 10 A run and dropped by the buffer
    # crop -- silently breaking the identical-row-set premise the multi-source
    # comparison rests on.
    d2_min = d2.min(axis=1)
    res_min_d = {}
    for j, i in enumerate(prot_idx):
        k = res_key(i)
        if k in keep_keys:
            v = float(np.sqrt(d2_min[j]))
            if v < res_min_d.get(k, np.inf):
                res_min_d[k] = v
    return lig_idx, groups, res_min_d


def build_pocket(d, groups, res_min_d):
    """Write pocket heavy atoms to a temp PDB with temporary single-char chain ids,
    protonate with PDBFixer, relax H, and map every output atom back to its Protenix
    join key. Returns dict with symbols, coords, charge and per-atom join keys."""
    an = d["atom_array"]._annot
    coord = d["atom_array"]._coord

    # order residues by (chain, seq) so the PDB has a sane, reproducible layout
    keys = sorted(groups, key=lambda k: (k[0], k[1], k[3], k[4]))
    chains = sorted({(k[0], k[1]) for k in keys})
    if len(chains) > len(TMP_CHAINS):
        raise ValueError("too many pocket chains for the PDB round-trip")
    tmp_of = {c: TMP_CHAINS[i] for i, c in enumerate(chains)}

    # residue_meta[(tmp_chain, tmp_resseq)] = (join key, {atom_name: source index})
    residue_meta, lines, serial = {}, [], 0
    per_chain_counter = {}
    for k in keys:
        tmp_c = tmp_of[(k[0], k[1])]
        per_chain_counter[tmp_c] = per_chain_counter.get(tmp_c, 0) + 1
        tmp_r = per_chain_counter[tmp_c]
        idxs = groups[k]
        residue_meta[(tmp_c, tmp_r)] = (k, {str(an["atom_name"][i]): int(i) for i in idxs})
        for i in idxs:
            serial += 1
            name = str(an["atom_name"][i])
            # PDB atom-name column convention: 4-char field, 1-char elements start at col 2
            nm = f" {name:<3s}" if len(name) < 4 else name
            x, y, z = coord[i]
            el = str(an["element"][i]).capitalize()
            lines.append(
                f"ATOM  {serial:5d} {nm}{str(an['res_name'][i]):>4s}{tmp_c:>2s}{tmp_r:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}")
    pdb_block = "\n".join(lines) + "\nEND\n"

    tmp = tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False,
                                      dir=os.environ.get("TMPDIR", "/tmp"))
    try:
        tmp.write(pdb_block)
        tmp.close()
        fixer = PDBFixer(filename=tmp.name)
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()          # MSE -> MET etc, required by ff14SB
        fixer.findMissingResidues()
        fixer.missingResidues = {}                  # never model unresolved loops
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()                     # terminal OXT / missing side-chain heavy
        fixer.addMissingHydrogens(7.4)
    finally:
        os.unlink(tmp.name)

    coords, charge = _minimize_frozen_heavy(fixer.topology, fixer.positions)

    # Map back by (chain order, residue order, atom name). PDBFixer preserves both orders
    # because we forbade missing-residue modelling.
    top = fixer.topology
    atoms = list(top.atoms())
    parent = _heavy_parent_map(top, atoms)

    n = len(atoms)
    syms = np.empty(n, dtype="<U2")
    key_asym = np.empty(n, dtype="<U8")
    key_sym = np.full(n, -1, dtype=np.int32)
    key_label_seq = np.full(n, -1, dtype=np.int32)
    key_auth_seq = np.full(n, -1, dtype=np.int32)
    key_ins = np.empty(n, dtype="<U4")
    key_res = np.empty(n, dtype="<U8")
    key_atom = np.empty(n, dtype="<U8")
    labelled = np.zeros(n, dtype=bool)
    parent_idx = np.arange(n, dtype=np.int32)
    # inf for anything with no source residue -> never inside a label crop
    res_min_lig = np.full(n, np.inf, dtype=np.float32)

    chain_ids = [c.id for c in top.chains()]
    res_seq_of = {}
    for ci, ch in enumerate(top.chains()):
        for ri, res in enumerate(ch.residues(), start=1):
            res_seq_of[res.index] = (chain_ids[ci], ri)

    for i, at in enumerate(atoms):
        syms[i] = at.element.symbol if at.element is not None else "X"
        ck = res_seq_of[at.residue.index]
        meta = residue_meta.get(ck)
        key_asym[i], key_sym[i] = "", -1
        key_ins[i], key_res[i], key_atom[i] = "", "", at.name
        if meta is None:
            continue
        k, name_map = meta
        key_asym[i], key_sym[i] = k[0], k[1]
        key_label_seq[i], key_auth_seq[i], key_ins[i] = k[2], k[3], k[4]
        key_res[i] = str(an["res_name"][next(iter(name_map.values()))])
        res_min_lig[i] = res_min_d.get(k, np.inf)
        if syms[i] != "H" and at.name in name_map:
            labelled[i] = True                       # a real Protenix heavy atom
        parent_idx[i] = parent[i]

    # H inherit their parent's join key so stage 2 can fold their energy in
    for i in range(n):
        if syms[i] == "H":
            p = parent_idx[i]
            key_asym[i], key_sym[i] = key_asym[p], key_sym[p]
            key_label_seq[i], key_auth_seq[i] = key_label_seq[p], key_auth_seq[p]
            key_ins[i], key_res[i], key_atom[i] = key_ins[p], key_res[p], key_atom[p]

    return dict(symbols=syms, coords=coords.astype(np.float32), charge=charge,
                key_asym=key_asym, key_sym=key_sym, key_label_seq=key_label_seq,
                key_auth_seq=key_auth_seq, key_ins=key_ins, key_res=key_res,
                key_atom=key_atom, labelled=labelled, parent=parent_idx,
                res_min_lig_dist=res_min_lig)


def _heavy_parent_map(topology, atoms):
    """For each atom index, the index of the heavy atom it belongs to (itself if heavy)."""
    parent = list(range(len(atoms)))
    index_of = {a: i for i, a in enumerate(atoms)}
    for a, b in topology.bonds():
        ia, ib = index_of[a], index_of[b]
        sa = a.element.symbol if a.element is not None else "X"
        sb = b.element.symbol if b.element is not None else "X"
        if sa == "H" and sb != "H":
            parent[ia] = ib
        elif sb == "H" and sa != "H":
            parent[ib] = ia
    return np.array(parent, dtype=np.int32)


def _relax_hydrogens(molH):
    """Relax only the added H against frozen heavy atoms. Tries MMFF then UFF.

    Both can fail on exotic/metal-containing CCD ligands -- missing parameters, or BFGS
    hitting a degenerate search direction ("bad direction in linearSearch") when almost
    every atom is fixed. RDKit's AddHs(addCoords=True) already places H at sensible bond
    lengths and angles, so falling back to that geometry is a mild degradation, whereas
    letting the exception through would drop the system entirely.
    """
    heavy = [a.GetIdx() for a in molH.GetAtoms() if a.GetAtomicNum() > 1]
    for make_ff in (
        lambda m: AllChem.MMFFGetMoleculeForceField(m, AllChem.MMFFGetMoleculeProperties(m)),
        AllChem.UFFGetMoleculeForceField,
    ):
        try:
            ff = make_ff(molH)
            if ff is None:
                continue
            for idx in heavy:
                ff.AddFixedPoint(idx)
            ff.Minimize(maxIts=500)
            return True
        except Exception:  # noqa: BLE001
            continue
    return False


def build_ligand(d, lig_idx):
    """RDKit mol straight from the Protenix BondList (bond orders + formal charges are
    already there, so no SMILES/CCD lookup is needed), AddHs, MMFF relax with heavy fixed."""
    aa = d["atom_array"]
    an, coord = aa._annot, aa._coord
    bonds = aa._bonds._bonds

    pos_of = {int(g): loc for loc, g in enumerate(lig_idx)}
    rw = Chem.RWMol()
    for g in lig_idx:
        el = str(an["element"][g]).capitalize()
        at = Chem.Atom(el)
        # Deliberately leave implicit-H perception ON: the bond orders + formal charges
        # below fully determine how many H each heavy atom is missing, and AddHs needs that
        # implicit count to place them. Setting NoImplicit here would make AddHs a no-op and
        # leave the ligand unprotonated (which then feeds MMFF radical-like valences).
        at.SetFormalCharge(int(round(float(an["ref_charge"][g]))))
        rw.AddAtom(at)
    for a, b, t in bonds:
        if int(a) in pos_of and int(b) in pos_of:
            rw.AddBond(pos_of[int(a)], pos_of[int(b)],
                       BOND_ORDER.get(int(t), Chem.BondType.SINGLE))
    mol = rw.GetMol()

    conf = Chem.Conformer(mol.GetNumAtoms())
    for g, loc in pos_of.items():
        x, y, z = coord[g]
        conf.SetAtomPosition(loc, Point3D(float(x), float(y), float(z)))
    mol.AddConformer(conf, assignId=True)
    Chem.SanitizeMol(mol)

    molH = Chem.AddHs(mol, addCoords=True)
    _relax_hydrogens(molH)

    cf = molH.GetConformer()
    n = molH.GetNumAtoms()
    syms = np.array([a.GetSymbol() for a in molH.GetAtoms()], dtype="<U2")
    coords = np.array([list(cf.GetAtomPosition(i)) for i in range(n)], dtype=np.float32)

    loc_to_g = {loc: g for g, loc in pos_of.items()}
    parent = np.arange(n, dtype=np.int32)
    for a in molH.GetAtoms():
        if a.GetAtomicNum() == 1:
            nbrs = [x.GetIdx() for x in a.GetNeighbors()]
            if nbrs:
                parent[a.GetIdx()] = nbrs[0]

    key_asym = np.empty(n, dtype="<U8")
    key_sym = np.full(n, -1, dtype=np.int32)
    key_label_seq = np.full(n, -1, dtype=np.int32)
    key_auth_seq = np.full(n, -1, dtype=np.int32)
    key_ins = np.empty(n, dtype="<U4")
    key_res = np.empty(n, dtype="<U8")
    key_atom = np.empty(n, dtype="<U8")
    labelled = np.zeros(n, dtype=bool)
    for i in range(n):
        g = loc_to_g.get(int(parent[i]))
        if g is None:
            continue
        key_asym[i] = str(an["label_asym_id"][g])
        key_sym[i] = int(an["sym_id_int"][g])
        key_label_seq[i] = _as_int(an["label_seq_id"][g])
        key_auth_seq[i] = _as_int(an["auth_seq_id"][g])
        key_ins[i] = str(an["ins_code"][g])
        key_res[i] = str(an["res_name"][g])
        key_atom[i] = str(an["atom_name"][g])
        labelled[i] = syms[i] != "H"

    return dict(symbols=syms, coords=coords, charge=Chem.GetFormalCharge(molH),
                key_asym=key_asym, key_sym=key_sym, key_label_seq=key_label_seq,
                key_auth_seq=key_auth_seq, key_ins=key_ins, key_res=key_res,
                key_atom=key_atom, labelled=labelled, parent=parent,
                # the ligand DEFINES the crop, so its distance to itself is 0
                res_min_lig_dist=np.zeros(n, dtype=np.float32))


def main():
    pdb_id, assembly_id, lig_asym, lig_sym, out_npz = sys.argv[1:6]
    lig_sym = int(lig_sym)

    d = load_assembly(pdb_id)
    lig_idx, groups, res_min_d = select_system(d, lig_asym, lig_sym)
    pocket = build_pocket(d, groups, res_min_d)
    ligand = build_ligand(d, lig_idx)

    n_tot = len(pocket["symbols"]) + len(ligand["symbols"])
    if n_tot > MAX_ATOMS:
        raise ValueError(f"system too large: {n_tot} atoms > MAX_ATOMS={MAX_ATOMS}")

    # concatenated complex: pocket block first, then ligand block, identical geometry
    cat = lambda k: np.concatenate([pocket[k], ligand[k]])  # noqa: E731
    is_lig = np.concatenate([np.zeros(len(pocket["symbols"]), bool),
                             np.ones(len(ligand["symbols"]), bool)])
    parent = np.concatenate([pocket["parent"],
                             ligand["parent"] + len(pocket["symbols"])]).astype(np.int32)

    np.savez_compressed(
        out_npz,
        pdb_id=pdb_id, assembly_id=str(assembly_id),
        lig_asym=lig_asym, lig_sym=lig_sym,
        symbols=cat("symbols"), coords=cat("coords"),
        is_ligand=is_lig, parent=parent, labelled=cat("labelled"),
        # per-atom min distance to the ligand, from the ORIGINAL AtomArray.
        # A downstream --label-cut is a pure THRESHOLD on this, which is exact
        # even for residues PDBFixer rewrites (MSE -> MET drops SE).
        res_min_lig_dist=cat("res_min_lig_dist"),
        key_asym=cat("key_asym"), key_sym=cat("key_sym"),
        key_label_seq=cat("key_label_seq"), key_auth_seq=cat("key_auth_seq"),
        key_ins=cat("key_ins"), key_res=cat("key_res"), key_atom=cat("key_atom"),
        pocket_charge=np.int32(pocket["charge"]), ligand_charge=np.int32(ligand["charge"]),
        n_pocket=np.int32(len(pocket["symbols"])), n_ligand=np.int32(len(ligand["symbols"])),
        pocket_cut=np.float32(POCKET_CUT),
    )
    print(f"{pdb_id}_{assembly_id}_{lig_asym}{lig_sym}: OK "
          f"pocket={len(pocket['symbols'])} (q={pocket['charge']}) "
          f"ligand={len(ligand['symbols'])} (q={ligand['charge']}) "
          f"labelled={int(cat('labelled').sum())}", flush=True)


if __name__ == "__main__":
    sys.exit(main())
