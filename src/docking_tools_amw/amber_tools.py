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
    Get the total charge of a protein using tleap from AmberTools. Calls to 
    sed, grep, and awk are used to extract the charges from the mol2 file
    and handle special characters/cases.
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


def prepare_protein_QM(pdb_path, output_pdb_path: str=None):
    if output_pdb_path is None:
        output_pdb_path = pdb_fname.replace(".pdb", "_hfixed.pdb")
    # Need to run 'reduce' to add hydrogens
    cmd = ["reduce", "-FLIP", pdb_fname]
    with open(output_pdb_path, "w") as file, open("error.log", "w") as error_file:
        result = subprocess.run(cmd, stdout=file, stderr=error_file)
    # Need to add missing residues/atoms
    fixer = PDBFixer(filename=output_pdb_path)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=pH)
    PDBFile.writeFile(fixer.topology, fixer.positions, open(output_pdb_path, "w"))
    # Here we need to fix the HIS and CYS residues to have correct resnames for tleap
    fix_CYS_HIS_cpptraj(output_pdb_path, output_pdb_path, identifier=n)
    # Get the charge of the protein with tleap
    total_charge, resnum_charge_dict = get_amber_charge(output_pdb_path)
    return total_charge, resnum_charge_dict

