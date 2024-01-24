import MDAnalysis as mda
import warnings
import pandas as pd
warnings.filterwarnings("ignore")


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


def split_pdb_into_components(pdb_path, pdb_id=None, count=None, verbose=0):
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
