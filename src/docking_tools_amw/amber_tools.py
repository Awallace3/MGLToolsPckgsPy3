import subprocess
import os


def fix_CYS_HIS_cpptraj(pdb_path, output_path=None, identifier="", prepareforleap_args=['nosugar']):
    """
    Fix CYS and HIS residues to have right resname and sulfur/disulfide bridges
    in a PDB file using cpptraj's prepareforleap command.
    """
    if output_path is None or output_path == pdb_path:
        pdb_out_path = pdb_path.replace('.pdb','_pfl.pdb')
    else:
        pdb_out_path = output_path
    cpptraj_cmd = f"""
parm {pdb_path}
loadcrd {pdb_path} name tmp1
prepareforleap crdset tmp1 name tmp2 pdbout {pdb_out_path} {" ".join(prepareforleap_args)}
"""
    with open(f'cpptraj{identifier}.in', 'w') as f:
        f.write(cpptraj_cmd)
    cmd = f'cpptraj -i cpptraj{identifier}.in > cpptraj{identifier}.in'
    out = subprocess.run(cmd, shell=True, check=True)
    if out.returncode != 0:
        raise RuntimeError("prepareforleap failed")
    if output_path == pdb_path:
        os.system(f'mv {pdb_out_path} {pdb_path}')
    os.system('rm cpptra*.in cpptra*.out')
    return


def get_amber_charge(PRO_pdb_path):
    """
    Get the total charge of a protein using tleap from AmberTools.
    """
    def_dir = os.getcwd()
    path_to_pdb = "/".join(PRO_pdb_path.split("/")[:-1])
    os.chdir(path_to_pdb)
    pdb_name = os.path.basename(PRO_pdb_path)
    pdb_name_no_ext = ".".join(pdb_name.split(".")[:-1])
    mol2_path = pdb_name.replace("pdb", "mol2")
    tleap_in = f"""
source leaprc.protein.ff19SB
source leaprc.water.opc
mol = loadPdb {pdb_name}
savemol2 mol {mol2_path} 1
quit
"""
    with open(f"tleap_{pdb_name_no_ext}.in", "w") as f:
        f.write(tleap_in)
    cmd = f"tleap -f tleap.in > tleap_{pdb_name_no_ext}.dat"
    out = subprocess.run(cmd, shell=True, check=True)
    if out.returncode != 0:
        os.chdir(def_dir)
        raise RuntimeError("tleap failed")
    with open(mol2_path, 'rb') as input:
        raw_data = input.read()
    mol2_data = raw_data.decode('utf-8', 'ignore').split('\n')
    # mol2_data = decode_file_into_utf8(mol2_path).split('\n')
    with open(mol2_path, "w") as output:
        output.write("\n".join(mol2_data))
    # print(path_to_pdb, mol2_path, sep="/")
    out_data, write = [], False
    for count, line in enumerate(mol2_data):
        if line.strip() == "@<TRIPOS>ATOM":
            write = True
        if line.strip() == "@<TRIPOS>BOND":
            write = False
        if write == True and line.strip() != "@<TRIPOS>ATOM":
            out_data.append(line)
    charges_fn = f"charges_{pdb_name_no_ext}.txt"

    with open(charges_fn, "w") as output:
        output.write("\n".join(out_data))

    total_charge = 0
    updated_charge = False
    # Perhaps awk for charge column and sum. Grep can filter out non-utf8 characters
    # os.system(f"awk '{{print $3, $9}}' {charges_fn} > charges_{pdb_name_no_ext}_clean.txt")
    with open(charges_fn) as file:
        lines = file.readlines()
        init_resi = lines[0].split()[-3]
        charge = 0
        for line in lines:
            resi_name = line.split()[-3]
            if resi_name == init_resi:
                charge += float(line.split()[-2])
                updated_charge = True
            else:
                total_charge += charge
                init_resi = resi_name
                updated_charge = True
                charge = float(line.split()[-2])

    if updated_charge:
        os.system(f"rm charges*.txt {mol2_path} tleap*.in tleap*.dat tleap*.log")
        total_charge = round(total_charge)
    os.chdir(def_dir)
    return total_charge
