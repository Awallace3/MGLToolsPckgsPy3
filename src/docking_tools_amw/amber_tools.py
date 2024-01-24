import subprocess
import os
from pprint import pprint as pp


def fix_CYS_HIS_cpptraj(pdb_path, output_path=None, identifier="", prepareforleap_args=['nosugar']):
    """
    Fix CYS and HIS residues to have right resname and sulfur/disulfide bridges
    in a PDB file using cpptraj's prepareforleap command.

    NOTE: the identifier is used to name the cpptraj input and output files to
          avoid issue with parallelization.
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
    cpptraj_in = f'cpptraj{identifier}.in'
    cpptraj_out = f'cpptraj{identifier}.out'
    with open(cpptraj_in, 'w') as f:
        f.write(cpptraj_cmd)
    cmd = f'cpptraj -i {cpptraj_in} > {cpptraj_out}'
    out = subprocess.run(cmd, shell=True, check=True)
    if out.returncode != 0:
        raise RuntimeError("prepareforleap failed")
    if output_path == pdb_path:
        os.system(f'mv {pdb_out_path} {pdb_path}')
    os.system(f'rm {cpptraj_in} {cpptraj_out}')
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
    # Run tleap
    tleap_in_fn = f"tleap_{pdb_name_no_ext}.in"
    with open(tleap_in_fn, "w") as f:
        f.write(tleap_in)
    cmd = f"tleap -f {tleap_in_fn} > tleap_{pdb_name_no_ext}.dat"
    out = subprocess.run(cmd, shell=True, check=True)
    # check if tleap failed
    if out.returncode != 0:
        os.chdir(def_dir)
        raise RuntimeError("tleap failed")

    # Extract charges from specific section in mol2 file
    start_linenumber = subprocess.run(f"grep -n '@<TRIPOS>ATOM' {mol2_path} | cut -d: -f1", shell=True, check=True, capture_output=True)
    start_linenumber = int(start_linenumber.stdout.decode('utf-8').strip())
    end_linenumber = subprocess.run(f"grep -n '@<TRIPOS>BOND' {mol2_path} | cut -d: -f1", shell=True, check=True, capture_output=True)
    end_linenumber = int(end_linenumber.stdout.decode('utf-8').strip())
    # Need to get resname and atom amber charge; however, need to be careful of special characters in mol2 file
    cmd = f"sed -n '{start_linenumber},{end_linenumber}p' {mol2_path} | sed '1d;$d' | sed 's/[^a-zA-Z0-9 +-\\*\\.]//g' | awk '{{print $8, $(NF - 1)}}'"
    out = subprocess.run(cmd, shell=True, check=True, capture_output=True)
    resname_charge = out.stdout.decode('utf-8', 'ignore').strip().split("\n")
    # set counters
    total_charge, res_charge, resnum = 0, 0, 0
    resnum_charge_dict = {}
    for n, i in enumerate(resname_charge):
        i = i.strip().split()
        if len(i) == 2:
            resname, charge = i
            if len(resname) < 3:
                resname = init_resname
        elif len(i) == 1:
            charge = i[0]
            resname = init_resname
        else:
            raise ValueError("Error in parsing charges")
        charge = float(charge)
        if n == 0:
            init_resname, res_charge = resname, charge
            resnum += 1
            resname = f"{resname}{resnum}"
            total_charge += charge
        elif resname == init_resname:
            resname = f"{resname}{resnum}"
            total_charge += charge
            res_charge += charge
            resnum_charge_dict[resname] = res_charge
        else:
            init_resname, res_charge = resname, charge
            resname = f"{resname}{resnum}"
            total_charge += charge
            resnum += 1
    os.system(f"rm {mol2_path} {tleap_in_fn} tleap_{pdb_name_no_ext}.dat *.log")
    total_charge = round(total_charge)
    os.chdir(def_dir)
    return total_charge, resnum_charge_dict


def get_charge(PRO_pdb_path):
    def_dir = os.getcwd()
    path_to_pdb = "/".join(PRO_pdb_path.split("/")[:-1])
    os.chdir(path_to_pdb)
    pdb_name = os.path.basename(PRO_pdb_path)
    pdb_name_no_ext = ".".join(pdb_name.split(".")[:-1])
    mol2_path = pdb_name.replace("pdb", "mol2")
    with open(f"tleap_{pdb_name_no_ext}.in", "w") as f:
        f.write("source leaprc.protein.ff19SB\nsource leaprc.water.opc\n\n")
        f.write(f'mol = loadPdb "{pdb_name}"\n')
        f.write(f"savemol2 mol {mol2_path} 1\n")
        f.write("quit")
    cmd = f"tleap -f tleap_{pdb_name_no_ext}.in > tleap_{pdb_name_no_ext}.dat"
    out = subprocess.run(cmd, shell=True, check=True)
    if out.returncode != 0:
        os.chdir(def_dir)
        raise RuntimeError("tleap failed")
    with open(mol2_path, 'rb') as input:
        raw_data = input.read()
    mol2_data = raw_data.decode('utf-8', 'ignore').split('\n')
    with open(mol2_path, "w") as output:
        output.write("\n".join(mol2_data))
    # print(path_to_pdb, mol2_path, sep="/")
    out_data, write = [], False
    for count, line in enumerate(mol2_data):
        if line.strip() == "@<TRIPOS>ATOM":
            write = True
        elif write == True:
            out_data.append(line)
        elif line.strip() == "@<TRIPOS>BOND":
            write = False
    charges_fn = f"charges_{pdb_name_no_ext}.txt"
    print(charges_fn)

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
        # os.system(f"rm {charges_fn} {mol2_path} tleap_{pdb_name_no_ext}.in tleap_{pdb_name_no_ext}.in *.log")
        total_charge = round(total_charge)
    os.chdir(def_dir)
    return total_charge
