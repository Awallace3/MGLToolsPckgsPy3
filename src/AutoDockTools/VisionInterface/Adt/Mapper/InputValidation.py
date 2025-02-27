########################################################################
#
#    Vision Node - Python source code - file generated by vision
#    Monday 22 June 2009 13:08:37 
#    
#       The Scripps Research Institute (TSRI)
#       Molecular Graphics Lab
#       La Jolla, CA 92037, USA
#
# Copyright: Daniel Stoffler, Michel Sanner and TSRI
#   
# revision: Guillaume Vareille
#  
#########################################################################
#
# $Header: /mnt/raid/services/cvs/python/packages/share1.5/AutoDockTools/VisionInterface/Adt/Mapper/InputValidation.py,v 1.1 2010/07/26 23:00:19 jren Exp $
#
# $Id: InputValidation.py,v 1.1 2010/07/26 23:00:19 jren Exp $
#

# import node's base class node
from NetworkEditor.items import NetworkNode
class InputValidation(NetworkNode):
    """
    Get inputs for virtual screening, and check for validity of inputs.
    """
    mRequiredTypes = {}
    mRequiredSynonyms = [
    ]

    def __init__(self, constrkw = {},  name='ADTFileNames', **kw):
        kw['constrkw'] = constrkw
        kw['name'] = name
        NetworkNode.__init__(*(self,), **kw)

        ip = self.inputPortsDescr
        ip.append({'name': 'recpetor_obj', 'datatype': 'receptor'})

        op = self.outputPortsDescr
        op.append({'name': 'GPF_template', 'datatype': 'gpf_template'})
        op.append({'name': 'DPF_template', 'datatype': 'dpf_template'})
        op.append({'name': 'result_dir', 'datatype': 'string'})


        code = """def doit(self, receptor_obj):
        import os
        from AutoDockTools.VisionInterface.Adt.receptor import receptor
        from AutoDockTools.VisionInterface.Adt.dpf_template import dpf_template
        from AutoDockTools.VisionInterface.Adt.gpf_template import gpf_template

        receptor_id = receptor_obj.get_id()
        receptor_dir = receptor_obj.get_workdir() 

        gpf_file = receptor_dir + os.sep + receptor_id + '.gpf'
        dpf_file = receptor_dir + os.sep + receptor_id + '.dpf'
        GPF_template = gpf_template(gpf_file)
        DPF_template = dpf_template(dpf_file)
        result_dir = os.path.abspath(receptor_dir + os.sep + '..' + os.sep + receptor_id)

        if not(os.path.exists(gpf_file)):
            print "ERROR: GPF template " + gpf_file + " does not exist!"
            return 'stop'
        elif not(os.path.exists(dpf_file)):
            print "ERROR: DPF template " + dpf_file + " does not exist!"
            return 'stop'

        pdbqt_loc = receptor_obj.get_ext_loc('pdbqt')
        pqr_loc = receptor_obj.get_ext_loc('pqr')
        pdb_loc = receptor_obj.get_ext_loc('pdb')
            
        if pdbqt_loc == None and pqr_loc == None and pdb_loc == None:
            print "ERROR: No valid structure file found, none of the following exist"
            print "    " + pdbqt_loc + ", " + pqr_loc + ", " + pdb_loc
            return 'stop'
        
        print "-------------------------------------------------------"
        print "     INPUTS THAT WILL BE USED FOR VIRTUAL SCREENING    "
        print "GPF Template:                  " + GPF_template.fullpath
        print "DPF Template:                  " + DPF_template.fullpath
        print "Results will be downloaded to: " + result_dir
        print "-------------------------------------------------------"
        
        pass
        self.outputData(GPF_template=GPF_template, DPF_template=DPF_template, result_dir=result_dir)
"""
        self.configure(function=code)


    def beforeAddingToNetwork(self, net):
        try:
            ed = net.getEditor()
        except:
            import traceback; traceback.print_exc()
            print('Warning! Could not import widgets')

