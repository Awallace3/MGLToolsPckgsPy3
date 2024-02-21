import MDAnalysis as mda
import warnings
import pandas as pd
import numpy as np
import os
import qcelemental as qcel
warnings.filterwarnings("ignore")


def mda_universe_to_xyz(u):
    u_xyz = u.atoms.positions
    u_elements = u.atoms.elements
    u_elements = [qcel.periodictable.to_Z(i) for i in u_elements]
    g = np.concatenate((np.reshape(u_elements, (-1, 1)), u_xyz), axis=1)
    return g


def read_pdbqt(lig_pdbqt, model=1, identifier=1):
    if not os.path.exists(lig_pdbqt):
        raise ValueError(f"File {lig_pdbqt} not found")
    # sed remove all lines outside of MODEL {model} to ENDMDL
    tmp_file = f"tmp{identifier}.pdbqt"
    os.system(f"sed -n '/^MODEL 1$/,/^ENDMDL/p' {lig_pdbqt} > {tmp_file}")
    u = mda.Universe(tmp_file)
    # os.system(f"rm {tmp_file}")
    return u.atoms.positions


def convert_bio_to_pdb(bf):
    num_count = bf.split(".bio")[-1]
    pdb_path = ".".join(bf.split(".")[:-1]) + "_" + num_count + ".pdb"
    with open(bf, "r") as f:
        lines = f.readlines()
    crystal_line_found = False
    new_lines = []
    for n, l in enumerate(lines):
        if l.startswith("CONECT"):
            continue
        if "CRYST1" in l:
            crystal_line_found = True
        if crystal_line_found and "CRYST1" in l:
            if lines[n + 1].startswith("MODEL"):
                continue
        new_lines.append(l)
    with open(pdb_path, "w") as f:
        f.writelines(new_lines)
    return pdb_path, num_count


def split_pdb_into_components(pdb_path, pdb_id=None, count=None, verbose=0, version='_pqr'):
    warnings.filterwarnings("ignore")
    if pdb_id is None:
        pdb_path_basename = ".".join(pdb_path.split(".")[:-1])
    else:
        pdb_path_basename = "/".join(pdb_path.split("/")[:-1]) + f"/{pdb_id.upper()}"

    if count is None:
        count = ""
    else:
        count = f"_{count}"
    print("Count: ", count)

    if verbose:
        print(f"PDB: {pdb_path}")
        print(f"Basename: {pdb_path_basename}")
    # pdb = mda.Universe(pdb_path, format='pdb')
    pdb = mda.Universe(pdb_path)
    protein = pdb.select_atoms("protein")
    if verbose:
        print(f"Protein: {protein.n_atoms}, ")
    with mda.Writer(f"{pdb_path_basename}_pro{count}.pdb", protein.n_atoms) as W:
        W.write(protein)
    waters = pdb.select_atoms("resname HOH")
    if verbose:
        print(f"Waters: {waters.n_atoms}")
    with mda.Writer(f"{pdb_path_basename}_wat{count}.pdb", waters.n_atoms) as W:
        W.write(waters)
    others = pdb.select_atoms("not resname HOH and not protein")
    chains = others.segments.segids
    if verbose:
        print(f"Others: {others.n_atoms}")
        print(f"Chains: {chains}")
    for c in chains:
        chain = others.select_atoms(f"segid {c}")
        with mda.Writer(
            f"{pdb_path_basename}_other_{c}{count}.pdb", chain.n_atoms
        ) as W:
            W.write(chain)
    return

def pdb_to_componets_mixed(pdb_path, pdb_others_path, pdb_base_out, ligand_resname=None, override=False, verbose=1):
    """
    Takes a pdb file and splits it into protein; protein and water; protein and other; and protein, water, and other
    """
    warnings.filterwarnings("ignore")
    if verbose:
        print(f"PDB: {pdb_path}")
        print(f"pdb_others: {pdb_others_path}")
        print(f"Basename: {pdb_base_out}")
    pro_out = f"{pdb_base_out}_pro.pdb"
    pro_wat_out = f"{pdb_base_out}_pro_wat.pdb"
    pro_oth_out = f"{pdb_base_out}_pro_oth.pdb"
    pro_wat_oth_out = f"{pdb_base_out}_pro_oth_wat.pdb"
    if ligand_resname is not None:
        lig_out = f"{pdb_base_out}_lig.pdb"

    pdb = mda.Universe(pdb_path)
    protein = pdb.select_atoms("protein")
    with mda.Writer(pro_out, protein.n_atoms) as W:
        W.write(protein)
    pro_wat = pdb.select_atoms("protein or resname HOH or resname WAT or resname TIP3")
    if verbose:
        print(f"Protein: {protein.n_atoms}; Waters: {pro_wat.n_atoms}")
    with mda.Writer(pro_wat_out, pro_wat.n_atoms) as W:
        W.write(pro_wat)
    if pdb_others_path is not None:
        others = mda.Universe(pdb_others_path)
        pro_oth = mda.Merge(protein, others.atoms)
        with mda.Writer(pro_oth_out, pro_oth.atoms.n_atoms) as W:
            W.write(pro_oth)
        pro_wat_oth = mda.Merge(pro_wat, others.atoms)
        with mda.Writer(pro_wat_oth_out, pro_wat_oth.atoms.n_atoms) as W:
            W.write(pro_wat_oth)
        pro_oth_atom_count = get_atom_count_pdb(pro_oth_out)
        pro_wat_oth_atom_count = get_atom_count_pdb(pro_wat_oth_out)
    else:
        pro_oth_out = None
        pro_wat_oth_out = None
        pro_oth_atom_count = None
        pro_wat_oth_atom_count = None
    if ligand_resname is not None:
        ligand = pdb.select_atoms(f"resname {ligand_resname}")
        with mda.Writer(lig_out, ligand.n_atoms) as W:
            W.write(ligand)
        ligand_atom_count = get_atom_count_pdb(lig_out)
    else:
        lig_out = None
        ligand_atom_count = None

    pro_atom_count = get_atom_count_pdb(pro_out)
    pro_wat_atom_count = get_atom_count_pdb(pro_wat_out)
    return pro_out, pro_wat_out, pro_oth_out, pro_wat_oth_out, pro_atom_count, pro_wat_atom_count, pro_oth_atom_count, pro_wat_oth_atom_count, lig_out, ligand_atom_count


def get_atom_count_pdb(pdb_path):
    pdb = mda.Universe(pdb_path)
    return int(pdb.atoms.n_atoms)


def split_pdb_into_components_identify_ligand(
    pdb_path: str,
    pdb_id: str,
    ligand_resnames: [str],
    count: str = None,
    verbose=0,
):
    """
    Incomplete.
    """
    warnings.filterwarnings("ignore")
    if pdb_id is None:
        pdb_path_basename = ".".join(pdb_path.split(".")[:-1])
    else:
        pdb_path_basename = "/".join(pdb_path.split("/")[:-1]) + f"/{pdb_id.upper()}"

    if count is None:
        count = ""
    else:
        count = f"_{count}"


    if verbose:
        print(f"PDB: {pdb_path}")
        print(f"Basename: {pdb_path_basename}")
        print(f"Ligand: {ligand_resnames}")
    # pdb = mda.Universe(pdb_path, format='pdb')
    pdb = mda.Universe(pdb_path)
    protein = pdb.select_atoms("protein")
    if verbose:
        print(f"Protein: {protein.n_atoms}, ")
    with mda.Writer(f"{pdb_path_basename}_pro{count}.pdb", protein.n_atoms) as W:
        W.write(protein)
    waters = pdb.select_atoms("resname HOH")
    if verbose:
        print(f"Waters: {waters.n_atoms}")
    if len(waters) > 0:
        with mda.Writer(f"{pdb_path_basename}_wat{count}.pdb", waters.n_atoms) as W:
            W.write(waters)
    others = pdb.select_atoms("not resname HOH and not protein")
    chains = others.segments.segids
    if verbose:
        print(f"Others: {others.n_atoms}")
        print(f"Chains: {chains}")
    if len(chains) > 0:
        for c in chains:
            for lig in ligand_resnames:
                if verbose:
                    print(f"Chain: {c}, Ligand: {lig}")
                chain = others.select_atoms(f"segid {c}")
                ligand = chain.select_atoms(f"resname {lig}")
                if len(ligand) > 0:
                    with mda.Writer(
                        f"{pdb_path_basename}_lig_{c}_{lig}{count}.pdb", ligand.n_atoms
                    ) as W:
                        W.write(ligand)
                non_ligand = chain.select_atoms(f"not resname {lig}")
                if len(non_ligand) > 0:
                    with mda.Writer(
                        f"{pdb_path_basename}_oth_{c}_{lig}{count}.pdb", non_ligand.n_atoms
                    ) as W:
                        W.write(non_ligand)
    else:
        for lig in ligand_resnames:
            ligand = others.select_atoms(f"resname {lig}")
            if len(ligand) > 0:
                with mda.Writer(
                    f"{pdb_path_basename}_lig_{lig}{count}.pdb", ligand.n_atoms
                ) as W:
                    W.write(ligand)
            non_ligand = others.select_atoms(f"not resname {lig}")
            if len(non_ligand) > 0:
                with mda.Writer(
                    f"{pdb_path_basename}_oth_{lig}{count}.pdb", non_ligand.n_atoms
                ) as W:
                    W.write(non_ligand)
    return

def remove_overlapping_atoms(
    df,
    verbose=1,
    n_jobs=-1,
):
    from joblib import Parallel, delayed
    def check_pdb(pdb, pdb_lig):
        residues_to_remove = []
        ligand_uni = mda.Universe(pdb_lig)
        ligand_center = ligand_uni.atoms.center_of_mass()
        pro_wat_uni = mda.Universe(pdb)
        for res in pro_wat_uni.residues:
            if res.resname != "HOH":
                continue
            res_center = res.atoms.center_of_mass()
            # check if all atoms are beyond the cutoff
            if np.linalg.norm(res_center - ligand_center) < cutoff:
                continue
            # check if water's atoms are too close to other residues
            for res2 in pro_wat_uni.residues:
                if res==res2:
                    continue
                for atom in res.atoms:
                    for atom2 in res2.atoms:
                        if atom == atom2:
                            continue
                        distance = np.linalg.norm(atom.position - atom2.position)
                        if distance < cutoff_overlapping:
                            print(f"    {atom.name} {atom.index} {res.resid} {atom2.name} {atom2.index} {res2.resid} {distance}")
                            residues_to_remove.append(res)
                            break
        if len(residues_to_remove) == 0:
            return
        print(residues_to_remove)
        residue_to_remove = "not resid " + " and not resid ".join([str(i.resid) for i in residues_to_remove])
        print(residue_to_remove)
        pro_without_residues = pro_wat_uni.select_atoms(residue_to_remove)
        print(len(pro_wat_uni.atoms), len(pro_without_residues.atoms))
        # write the new pdb
        pro_without_residues.write(pdb)
        return
    print(df.columns)
    print('pdbs to check:', len(df))
    cutoff = 12
    cutoff_overlapping = 1
    def mda_check_and_remove_residues(data):
        n, i, cutoff, cutoff_overlapping = data
        print(f"{n+1}/{len(df)} : {i['pdb_id']}")
        try:
            check_pdb(i["proteinhswater_pdb"], i["lig_pdb_hs"])
            check_pdb(i["proteinhswaterother_pdb"], i["lig_pdb_hs"])
        except Exception as e:
            print(f"Error with {i['pdb_id']}: {e}")
        return 
    if n_jobs == 1:
        for n, row in df.iterrows():
            mda_check_and_remove_residues((n, row, cutoff, cutoff_overlapping))
    else:
        Parallel(n_jobs=n_jobs)(delayed(mda_check_and_remove_residues)((n, row, cutoff, cutoff_overlapping)) for n, row in df.iterrows())
    return 
