import os


def pdb2pqr_create(data):
    """
    Not finished implementation
    """
    try:
        if len(data) == 3:
            n, i, conn = data
        else:
            n, i = data
        pdb_start = i["pro_pdb_hs"].split("_hfixed")[0]
        pdb_id = pdb_start.split("/")[-1]
        count = pdb_id[4]
        pdb_id = pdb_id[:4]
        lig_resname, chain, resnum = i["lig_name"].split(":")
        print(
            f"{n+1}/{len(df)} : {i['pdb_id']} : {lig_resname} : pH {i['ph']} : pdb_start {pdb_start}"
        )
        # merge protein and water
        pdb_in = i["pro_pdb_hs"]
        pdb_uni = mda.Universe(pdb_in)
        pdb_uni = mda.Merge(pdb_uni.atoms, mda.Universe(i["lig_pdb_hs"]).atoms)
        if i["wat_pdb"] is not None:
            pdb_uni = mda.Merge(pdb_uni.atoms, mda.Universe(i["wat_pdb"]).atoms)
        if i['oth_pdb'] is not None:
            pdb_uni = mda.Merge(pdb_uni.atoms, mda.Universe(i["oth_pdb"]).atoms)
        # print number of water atoms
        print(f"    {len(pdb_uni.select_atoms('resname HOH'))} water atoms")
        pdb_in = f"{pdb_start}_{chain}_{resnum}_merged.pdb"
        pdb_uni.atoms.write(pdb_in)
        cmd = f"sed -i '/^REMARK/d' {pdb_in}"
        subprocess.run(cmd, shell=True, check=True)

        print(f"{pdb_in = }, {lig_resname = }")
        ph = i["ph"]
        pdb_out = f"{pdb_start}_hfixed_pqr.pdb"
        lig_pdb = i["lig_pdb_hs"]
        cmd = f"pdb2pqr --ff={ff} {pdb_in} --neutraln --neutralc --with-ph={ph} --pdb-out {pdb_out} {pdb_start}_{resnum}.pqr --quiet --log-level={log_level}"
        result = subprocess.run(cmd, shell=True)
        if result.returncode != 0:
            raise ValueError(f"Error pdb2pqr exited with {result.returncode}")
        (
            pro_out,
            pro_wat_out,
            pro_oth_out,
            pro_wat_oth_out,
            pro_atom_count,
            pro_wat_atom_count,
            pro_oth_atom_count,
            pro_wat_oth_atom_count,
        ) = docking_tools_amw.mda_tools.pdb_to_componets_mixed(
            pdb_out,
            i['oth_pdb'], # None,
            pdb_start + "_pqr",
        )
        if verbose:
            print(f"{pro_out = }\n{pro_wat_out = }\n{pro_oth_out = }\n{pro_wat_oth_out = }")
        # os.system(f"rm {pdb_start}_{resnum}.pqr {pdb_in}")
        os.system(f"rm {pdb_start}_{resnum}.pqr ")

        if pro_atom_count == pro_wat_atom_count or pro_wat_out is None:
            pro_wat_out = None
            pro_wat_charge = None
        pro_charge, _ = docking_tools_amw.amber_tools.get_amber_charge(pro_out, identifier=resnum)
        try:
            oth_charge, _ = docking_tools_amw.amber_tools.get_amber_charge(pro_oth_out, identifier=resnum)
            wat_oth_charge = oth_charge
        except Exception as e:
            oth_charge = None
            wat_oth_charge = None
            pass
        os.chdir(def_dir)
        pro_wat_charge = pro_charge
        # Update the database
        out = {
            "pl_id": i['pl_id'],
            "pro_out": pro_out,
            "pro_wat_out": pro_wat_out,
            "pro_oth_out": pro_oth_out,
            "pro_wat_oth_out": pro_wat_oth_out,
            "pro_charge": pro_charge,
            "pro_wat_charge": pro_wat_charge,
            "pro_oth_charge": oth_charge,
            "pro_wat_oth_charge": wat_oth_charge,
        }
        print(f"   {out = }")
        if len(data) == 3:
            cur = conn.cursor()
            cur.execute(f"""UPDATE bmoad.protein_ligand SET proteinhs_pdb = '{out['pro_out']}' WHERE pl_id={out['pl_id']};""")
            cur.execute(f"""UPDATE bmoad.protein_ligand SET proteinhs_charge = {out['pro_charge']} WHERE pl_id={out['pl_id']};""")
            if out['pro_wat_out'] is not None:
                cur.execute(f"""UPDATE bmoad.protein_ligand SET proteinhswater_pdb = '{out['pro_wat_out']}' WHERE pl_id={out['pl_id']};""")
            if out['pro_wat_charge'] is not None:
                cur.execute(f"""UPDATE bmoad.protein_ligand SET proteinhswater_charge = {out['pro_wat_charge']} WHERE pl_id={out['pl_id']};""")
            if out['pro_oth_out'] is not None:
                cur.execute(f"""UPDATE bmoad.protein_ligand SET proteinhsother_pdb = '{out['pro_oth_out']}' WHERE pl_id={out['pl_id']};""")
            if out['pro_wat_oth_out'] is not None:
                cur.execute(f"""UPDATE bmoad.protein_ligand SET proteinhswaterother_pdb = '{out['pro_wat_oth_out']}' WHERE pl_id={out['pl_id']};""")
            if out['pro_wat_oth_charge'] is not None:
                cur.execute(f"""UPDATE bmoad.protein_ligand SET proteinhswaterother_charge = {out['pro_wat_oth_charge']} WHERE pl_id={out['pl_id']};""")
            if out['pro_oth_charge'] is not None:
                cur.execute(f"""UPDATE bmoad.protein_ligand SET proteinhsother_charge = {out['pro_oth_charge']} WHERE pl_id={out['pl_id']};""")
            print(f"    {out['pl_id']} updated {out['pro_out']}")
            conn.commit()
        return out
    except Exception as e:
        os.chdir(def_dir)
        print(f"Error with {i['pl_id']}: {e}")
        if len(data) == 3:
            cur = conn.cursor()
            cur.execute(f"""UPDATE bmoad.protein_ligand SET pdb2pqr_parse_error = %s WHERE pl_id={i['pl_id']};""", (str(e),))
            conn.commit()
        return {
            "pl_id": i['pl_id'],
            "pdb2pqr_error": e,
        }
