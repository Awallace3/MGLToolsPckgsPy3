from pymol import cmd, stored
import os

this_file_for_pymol_initialization = os.path.abspath(__file__)

def pocketfrag(monoC = None, cutcap=False):
    name="pocketfrag"
    # Get ligand info, add to system
    cmd.select("lig_A", "bm. hetatm and not sol. and not metals")
    # Now that ligand has been digested, color it
    cmd.color("red", "lig_A")
    # Carve away monomer C
    if monoC is not None:
        # Selecting by distance to ligand?
        if "be." in monoC and "w." not in monoC:
            # Digest monoC declaration
            monoC = monoC.strip('"')
            cutoff = monoC.split()[-1]
            # Make pocket selection
            ## Select full residues or metal cofactors in the pocket
            # CTS: Running command below cuts C (of carbonyl group) and N. Following lines adjust the cut
            cmd.select('sys%s_B' % cutoff, "(byres not lig_A w. %s of lig_A) or (metals w. %s of lig_A) " % (cutoff, cutoff))
            # CTS expand B such that only alpha carbon -- carbon bonds are broken. This only handles when N is on the QM side.
            cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and elem N) xt. 1)" % (cutoff, cutoff)) 
            cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and elem C) xt. 1 and elem O)" % (cutoff, cutoff)) 
            #cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and elem S) xt. 1)" % (cutoff, cutoff)) 
            #cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and elem S) xt. 1)" % (cutoff, cutoff)) 
            cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and elem Se) xt. 1)" % (cutoff, cutoff)) 
            cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and elem Se) xt. 1)" % (cutoff, cutoff)) 
            cmd.select('sys%s_B' % cutoff, "sys%s_B + (sys%s_B (xt. 1 and elem H))" % (cutoff, cutoff)) 
            ## Select everything else, that's monomer C
            cmd.select('mono_C', "not sys%s_B and not lig_A" % cutoff)
            ## Exclude peptide bonds from border residues, i.e., B--C border is cutting across CA--(sidechain) bond
            # CTS expand C such that only alpha carbon -- carbon bonds are broken. This handles when N is on the MM side.
            cmd.select("mono_C", "mono_C + ((mono_C and elem S) xt. 1 and elem S)")
            cmd.select("mono_C", "mono_C + ((mono_C and elem S) xt. 2)")
            cmd.select("mono_C", "mono_C + ((mono_C and name CA) xt. 1 and elem H)")
            cmd.select("mono_C", "mono_C + ((mono_C and name CA) xt. 1 and elem N)")
            cmd.select("mono_C", "mono_C + ((mono_C and elem N) xt. 1)") 
            cmd.select("mono_C", "mono_C + ((mono_C and elem C) xt. 1 and elem O)") 
            # CTS: expand C so that endcaps aren't on fronteir regions
            cmd.select("mono_C", "mono_C + ((mono_C and elem C) xt. 3 and resn ACE)") 
            cmd.select("mono_C", "mono_C + ((mono_C and elem C) xt. 3 and resn NMA)") 
            # CTS now re-specify system B so that there are no overlapping atoms (in both system B and monoC)
            cmd.select("sys%s_B" % cutoff, "not mono_C and not lig_A")
            cmd.hide("sticks", "mono_C")
            cmd.show("sticks", "mono_C and sol.")
            cmd.show("lines", "mono_C and not sol.")
            # CTS identify any MM residues between two QM residues
            stored.Cresis = []
            # Get atoms in mono_C that are directly bound to system B
            cmd.select("boundary_cs", "mono_C and bound_to sys%s_B" % cutoff)
            # Add the residue numbers of these atoms to a list
            cmd.iterate("boundary_cs", "stored.Cresis.append(resi)")
            #print(stored.Cresis)
            # For each residue bound directly to system B
            for r in stored.Cresis:
                resis = []
                up = str(int(r)+1)
                down = str(int(r)-1)
                # get atoms that are in the neighboring residue and part of system B and not C or O
                cmd.select("nextto", "sys%s_B and resi %s and not name C and not name O" % (cutoff, up))
                cmd.select("nextto", "nextto + (sys%s_B and resi %s and not name C and not name O)" % (cutoff, down))
                #cmd.save('%s.pdb' % r, "nextto")
                # find where there is more than one neighboring QM residue
                stored.neighbor_resis = []
                cmd.iterate("nextto", "stored.neighbor_resis.append(resi)")
                neighbor_resis = set(stored.neighbor_resis)
                if len(neighbor_resis) == 2:
                    #print(neighbor_resis)
                    cmd.select('sys%s_B' % cutoff, "sys%s_B + resi %s" % (cutoff, r))
                    cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and resi %s and elem N) xt. 1)" % (cutoff, cutoff, r)) 
                    cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and resi %s and elem C) xt. 1 and elem O)" % (cutoff, cutoff, r)) 
                    cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and resi %s and elem S) xt. 1)" % (cutoff, cutoff, r)) 
                    cmd.select('sys%s_B' % cutoff, "sys%s_B + (br. (sys%s_B and elem S))" % (cutoff, cutoff)) 
                    cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and resi %s and elem Se) xt. 1)" % (cutoff, cutoff, r)) 
                    cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and resi %s and elem Se) xt. 1)" % (cutoff, cutoff, r)) 
                    cmd.select('sys%s_B' % cutoff, "sys%s_B + (sys%s_B (xt. 1 and elem H))" % (cutoff, cutoff)) 
                    cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and elem N) xt. 1)" % (cutoff, cutoff)) 
                    cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and elem C) xt. 1 and elem O)" % (cutoff, cutoff)) 
                    # print('Number of atoms in QM protein:', cmd.count_atoms('sys%s_B' % cutoff))
                    ## Select everything else, that's monomer C
                    cmd.select('mono_C', "not sys%s_B and not lig_A" % cutoff)
                    ## Exclude peptide bonds from border residues, i.e., B--C border is cutting across CA--(sidechain) bond
                    # CTS expand C such that only alpha carbon -- carbon bonds are broken. This handles when N is on the MM side.
                    cmd.select("mono_C", "mono_C + ((mono_C and elem S) xt. 1 and elem S)")
                    cmd.select("mono_C", "mono_C + ((mono_C and elem S) xt. 2)")
                    cmd.select("mono_C", "mono_C + ((mono_C and name CA) xt. 1 and elem H)")
                    cmd.select("mono_C", "mono_C + ((mono_C and name CA) xt. 1 and elem N)")
                    cmd.select("mono_C", "mono_C + ((mono_C and elem N) xt. 1)") 
                    cmd.select("mono_C", "mono_C + ((mono_C and elem C) xt. 1 and elem O)") 
                    # CTS: expand C so that endcaps aren't on fronteir regions
                    cmd.select("mono_C", "mono_C + ((mono_C and elem C) xt. 3 and resn ACE)") 
                    cmd.select("mono_C", "mono_C + ((mono_C and elem C) xt. 3 and resn NMA)") 
                    # CTS now re-specify system B so that there are no overlapping atoms (in both system B and monoC)
                    cmd.select("sys%s_B" % cutoff, "not mono_C and not lig_A")
                    cmd.hide("sticks", "mono_C")
                    cmd.show("sticks", "mono_C and sol.")
                    cmd.show("lines", "mono_C and not sol.")
            # loop through each residue at boundary (in monoC)
            # expand; if expansion returns two atoms in QM
        elif "be." in monoC and "w." in monoC:
            monoC = monoC.strip('"').split()
            B_cutoff = monoC[monoC.index('be.') + 1]
            C_cutoff = monoC[monoC.index('w.') + 1]
            multilevel(B_cutoff, C_cutoff, cutcap=bool(cutcap))
        else:
            cmd.select('mono_C', monoC.strip('"'))
        # Process selected monomer C
        #colorsele("mono_C")
    else:
        B_cutoff = ''
        # Make pocket selection
        cmd.select('sys%s_B' % B_cutoff, "not lig_A")
    
    cmd.select('qm_sys','sys%s_B + lig_A' % cutoff)
    cmd.remove('not qm_sys')
    cmd.save('session.pse')
    cmd.h_add()
    cmd.save(f'{args.PDB_file[:-4]}_QM.pdb')
    # print('Number of atoms in QM protein:', cmd.count_atoms('all' % cutoff))
    return

def protein_truncation(PDB_file, cutoff=5.0, output_pdb_filename="truncated_protein.pdb"):
    """
    Usage:
        num_atoms = pymol_tools.protein_truncation(pdb_file, cutoff=5.0, output_pdb_filename="testing.pdb")
    """
    cmd.load(f'{PDB_file}')
    cmd.do(f'run {this_file_for_pymol_initialization}')
    num_atoms = protein_truncation_pymol_(monoC=f"be. {cutoff}", output_pdb_filename=output_pdb_filename)
    return num_atoms

def protein_truncation_pymol_(monoC=None, output_pdb_filename="truncated_protein.pdb"):
    name="pocketfrag"
    # Get ligand info, add to system
    cmd.select("lig_A", "bm. hetatm and not sol. and not metals")
    # Now that ligand has been digested, color it
    cmd.color("red", "lig_A")
    # Selecting by distance to ligand?
    if "be." in monoC and "w." not in monoC:
        # Digest monoC declaration
        monoC = monoC.strip('"')
        cutoff = monoC.split()[-1]
        # Make pocket selection
        ## Select full residues or metal cofactors in the pocket
        # CTS: Running command below cuts C (of carbonyl group) and N. Following lines adjust the cut
        cmd.select('sys%s_B' % cutoff, "(byres not lig_A w. %s of lig_A) or (metals w. %s of lig_A) " % (cutoff, cutoff))
        # CTS expand B such that only alpha carbon -- carbon bonds are broken. This only handles when N is on the QM side.
        cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and elem N) xt. 1)" % (cutoff, cutoff)) 
        cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and elem C) xt. 1 and elem O)" % (cutoff, cutoff)) 
        cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and elem Se) xt. 1)" % (cutoff, cutoff)) 
        cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and elem Se) xt. 1)" % (cutoff, cutoff)) 
        cmd.select('sys%s_B' % cutoff, "sys%s_B + (sys%s_B (xt. 1 and elem H))" % (cutoff, cutoff)) 
        ## Select everything else, that's monomer C
        cmd.select('mono_C', "not sys%s_B and not lig_A" % cutoff)
        ## Exclude peptide bonds from border residues, i.e., B--C border is cutting across CA--(sidechain) bond
        # CTS expand C such that only alpha carbon -- carbon bonds are broken. This handles when N is on the MM side.
        cmd.select("mono_C", "mono_C + ((mono_C and elem S) xt. 1 and elem S)")
        cmd.select("mono_C", "mono_C + ((mono_C and elem S) xt. 2)")
        cmd.select("mono_C", "mono_C + ((mono_C and name CA) xt. 1 and elem H)")
        cmd.select("mono_C", "mono_C + ((mono_C and name CA) xt. 1 and elem N)")
        cmd.select("mono_C", "mono_C + ((mono_C and elem N) xt. 1)") 
        cmd.select("mono_C", "mono_C + ((mono_C and elem C) xt. 1 and elem O)") 
        # CTS: expand C so that endcaps aren't on fronteir regions
        cmd.select("mono_C", "mono_C + ((mono_C and elem C) xt. 3 and resn ACE)") 
        cmd.select("mono_C", "mono_C + ((mono_C and elem C) xt. 3 and resn NMA)") 
        # CTS now re-specify system B so that there are no overlapping atoms (in both system B and monoC)
        cmd.select("sys%s_B" % cutoff, "not mono_C and not lig_A")
        cmd.hide("sticks", "mono_C")
        cmd.show("sticks", "mono_C and sol.")
        cmd.show("lines", "mono_C and not sol.")
        # CTS identify any MM residues between two QM residues
        stored.Cresis = []
        # Get atoms in mono_C that are directly bound to system B
        cmd.select("boundary_cs", "mono_C and bound_to sys%s_B" % cutoff)
        # Add the residue numbers of these atoms to a list
        cmd.iterate("boundary_cs", "stored.Cresis.append(resi)")
        # For each residue bound directly to system B
        for r in stored.Cresis:
            resis = []
            up = str(int(r)+1)
            down = str(int(r)-1)
            # get atoms that are in the neighboring residue and part of system B and not C or O
            cmd.select("nextto", "sys%s_B and resi %s and not name C and not name O" % (cutoff, up))
            cmd.select("nextto", "nextto + (sys%s_B and resi %s and not name C and not name O)" % (cutoff, down))
            #cmd.save('%s.pdb' % r, "nextto")
            # find where there is more than one neighboring QM residue
            stored.neighbor_resis = []
            cmd.iterate("nextto", "stored.neighbor_resis.append(resi)")
            neighbor_resis = set(stored.neighbor_resis)
            if len(neighbor_resis) == 2:
                #print(neighbor_resis)
                cmd.select('sys%s_B' % cutoff, "sys%s_B + resi %s" % (cutoff, r))
                cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and resi %s and elem N) xt. 1)" % (cutoff, cutoff, r)) 
                cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and resi %s and elem C) xt. 1 and elem O)" % (cutoff, cutoff, r)) 
                cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and resi %s and elem S) xt. 1)" % (cutoff, cutoff, r)) 
                cmd.select('sys%s_B' % cutoff, "sys%s_B + (br. (sys%s_B and elem S))" % (cutoff, cutoff)) 
                cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and resi %s and elem Se) xt. 1)" % (cutoff, cutoff, r)) 
                cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and resi %s and elem Se) xt. 1)" % (cutoff, cutoff, r)) 
                cmd.select('sys%s_B' % cutoff, "sys%s_B + (sys%s_B (xt. 1 and elem H))" % (cutoff, cutoff)) 
                cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and elem N) xt. 1)" % (cutoff, cutoff)) 
                cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and elem C) xt. 1 and elem O)" % (cutoff, cutoff)) 
                # print('Number of atoms in QM protein:', cmd.count_atoms('sys%s_B' % cutoff))
                ## Select everything else, that's monomer C
                cmd.select('mono_C', "not sys%s_B and not lig_A" % cutoff)
                ## Exclude peptide bonds from border residues, i.e., B--C border is cutting across CA--(sidechain) bond
                # CTS expand C such that only alpha carbon -- carbon bonds are broken. This handles when N is on the MM side.
                cmd.select("mono_C", "mono_C + ((mono_C and elem S) xt. 1 and elem S)")
                cmd.select("mono_C", "mono_C + ((mono_C and elem S) xt. 2)")
                cmd.select("mono_C", "mono_C + ((mono_C and name CA) xt. 1 and elem H)")
                cmd.select("mono_C", "mono_C + ((mono_C and name CA) xt. 1 and elem N)")
                cmd.select("mono_C", "mono_C + ((mono_C and elem N) xt. 1)") 
                cmd.select("mono_C", "mono_C + ((mono_C and elem C) xt. 1 and elem O)") 
                # CTS: expand C so that endcaps aren't on fronteir regions
                cmd.select("mono_C", "mono_C + ((mono_C and elem C) xt. 3 and resn ACE)") 
                cmd.select("mono_C", "mono_C + ((mono_C and elem C) xt. 3 and resn NMA)") 
                # CTS now re-specify system B so that there are no overlapping atoms (in both system B and monoC)
                cmd.select("sys%s_B" % cutoff, "not mono_C and not lig_A")
                cmd.hide("sticks", "mono_C")
                cmd.show("sticks", "mono_C and sol.")
                cmd.show("lines", "mono_C and not sol.")
    else:
        cmd.select('mono_C', monoC.strip('"'))

    cmd.select('qm_sys','sys%s_B + lig_A' % cutoff)
    cmd.remove('not qm_sys')
    cmd.save('session.pse')
    cmd.h_add()
    cmd.save(output_pdb_filename)
    number_of_QM_atoms = cmd.count_atoms('all')
    return number_of_QM_atoms

if __name__ == "__main__":
    import argparse
    cmd.extend("pocketfrag", pocketfrag)
    parser=argparse.ArgumentParser()
    parser.add_argument('PDB_file', type=str)
    parser.add_argument('cutoff', type=float)
    args=parser.parse_args()
    cmd.load(f'{args.PDB_file}')
    cmd.do('run cut_protein.py')
    cmd.do(f'pocketfrag split="ss", monoC="be. {args.cutoff}"')
    #cmd.do(f'make_dictionary {args.cutoff}')
