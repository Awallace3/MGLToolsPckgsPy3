#############################################################################
#
# Author: Ruth HUEY, Michel F. SANNER
#
# Copyright: M. Sanner TSRI 2000
#
#############################################################################

#
# $Header: /mnt/raid/services/cvs/python/packages/share1.5/AutoDockTools/autodpfCommands.py,v 1.100.2.5 2016/02/23 18:46:45 annao Exp $
#
# $Id: autodpfCommands.py,v 1.100.2.5 2016/02/23 18:46:45 annao Exp $
#
#
#
#
#
#

"""
This Module facilitates producing a docking parameter file for AutoDock. The steps in this process are:

    * Selecting the macromolecule filename: The user can select the macromolecule for autodpf in three ways: it can be chosen from molecules previously added to the moleculeViewer, it can be picked as a PDB file,  or it can be picked as a MOL2 file:

        o Choose Macromol...

        o Select PDB Macromolecule 

        o Select MOL2 Macromolecule

    * Selecting the small molecule which has been previously formatted by AutoTors: 

        o Via Reading a PDBQ-File which adds the ligand to the viewer

    * The user sets parameters pertaining to the small molecule 

        o Checking that a grid map exists for each of the ligand atom types 

        o Indicating whether a floating grid map exists

        o Setting the initial translation of the small molecule

            - by choosing  the 'random' option which sets a random starting position for the ligand

            - by entering the desired coordinates in the entry

        o Setting the initial quaternion of the small molecule

            - by choosing  the 'random' option which sets a random starting quaternion.

            - by entering the desired initial quaternion -Qx,Qy,Qz,Qw in the entry.  Qx, Qy, Qz define the unit vector in the direction of rigid body rotation and Qw the angle of rotation about this unit vector.

        o Setting the coefficient of the torsional DOF

        o By choosing to set the initial dihedrals for the small molecule or not: If not, AutoDock assumes that the chi1, chi2, chi3 etc are all zero and does not change the initial ligand torsion angles. If the user chooses to set the initial dihedrals, he further chooses:

            - for them to be randomly assigned 

            - an initial relative dihedral angle for each active torsion in the ligand.

        o The user can specify two types of torsion constraints for the ligand:

            -  Gaussian constraints which use an inverted Gaussian bell curve to calculate the energy function input of the constraint.  This type of constraint is specified by two floating point numbers: the perferred angle in the range -180-+180decreeds and the half-width which is the difference between two angles at which the energy is half the barrier PLUS an integer which identifies the torsion according to the list at the top of the AutoTors-generated input ligand PDBQ file. More than one constraint of this type may be specified for a single torsion.

            - Hard torsion constraints may also be specified. These differ from the previous type in that the torsion is never allowed to take values bewond the range defined and in that the second parameter is the full width of the allowed range of torsion angles. Moreover, only one constraint of this type is allowed per torsion.

        o If the user specifies torsion constraints, he may also specify the height of the energy barrier to be applied to these constraints.

        o If the user specifies Gaussian torsion constraints, he may also specify whether to store and output the torsion energies

    * The user sets parameters pertaining to docking algorithm(s) he wishes to use
:
        o Setting Simulated Annealing parameters.

        o Setting Genetic Algorithm parameters (GA).

        o Setting Local Search parameters (LS).

    It is important to remember that any of these may be used alone but only GA and LS may be used together


    * The user adjusts these additional parameters: 
    
        o the step sizes of translation, quaternion rotation and dihedral torsion change.
        o  energy parameters including energy assigned to atoms outside the grid volume, the maximum allowable initial energy and the maximum number of retries.

        o output format parameters including the level of detail for the output, the rms cluster tolerance, the reference file for rms calculations and whether to do symmetry checking in the rms calculations.


    * The user selects which kind of docking parameter file to write : 
    
        o Simulated Annealing 

        o GA

        o LS

        o GALS


    * The results of the previous steps are written to a file. The user selects a filename via a filebrowser.  By convention, the file should have a .dpf extension. If no macromolecule has been selected, it is not possible to write a grid parameter file and the user gets a warning message to that effect. Likewise, the types of the maps to be calculated must be set before the grid parameter file is written and a warning message to this effect appears if the types have not been set.
(A checkbutton, "DONE", allows the user to withdraw the autoTools menuBar)
    
"""

import os, numpy
import Pmw

from ViewerFramework.VFCommand import CommandGUI
from mglutil.gui.InputForm.Tk.gui import InputFormDescr
from mglutil.util.callback import CallBackFunction

from Pmv.mvCommand import MVCommand
from MolKit.tree import TreeNode, TreeNodeSet
from MolKit.molecule import AtomSet
from MolKit import Read
from Pmv.guiTools import MoleculeChooser
from tkinter.simpledialog import SimpleDialog
import types, string, tkinter,os
from .energyConstants import Rij, epsij, SolVol, SolPar, SolCon
from .DockingParameters import DockingParameters, simulated_annealing_list,\
    local_search_list, genetic_algorithm_list, cluster_list,\
    genetic_algorithm_local_search_list
from .DockingParameters import  simulated_annealing_list4,\
    local_search_list4, genetic_algorithm_list4,\
    genetic_algorithm_local_search_list4_with_parameter_file, \
    genetic_algorithm_local_search_list4
from .DockingParameters import  simulated_annealing_list4_1,\
    local_search_list4_1, genetic_algorithm_list4_1,\
    genetic_algorithm_local_search_list4_1_with_parameter_file, \
    genetic_algorithm_local_search_list4_1
from .DockingParameters import  simulated_annealing_list4_2,\
    local_search_list4_2, genetic_algorithm_list4_2,\
    genetic_algorithm_local_search_list4_2_with_parameter_file, \
    genetic_algorithm_local_search_list4_2, epdb_list4_2
from .DockingParameters import ConfigFileMaker
from .autotorsCommands import checkMolCharges



#these are the texts on menubuttons, menu entries etc:
menuText = {}
#menuText['AutoDpfMB'] = ' Set Docking Parameters '
menuText['AutoDpfMB'] = 'Docking'
menuText['ReadDpfMB'] = 'Open DPF...'

menuText['MacromoleculeMB'] = 'Macromolecule'
menuText['ReadMacro'] = 'Set Filename...(AD3)'
menuText['ChooseMacro'] = 'Choose...(AD3)'
menuText['ReadMacro4'] = 'Set Rigid Filename...'
menuText['ReadFlexRes4'] = 'Set Flexible Residues Filename...'
menuText['ChooseMacro4'] = 'Choose...'

menuText['SetLigandParmsMB'] = 'Ligand'
menuText['ReadLigand4'] = 'Open...'
menuText['ChooseLigand4'] = 'Choose...'
menuText['AdjustLigand4'] = 'Ligand Parameters...'
menuText['ReadLigand'] = 'Open...(AD3)'
menuText['ChooseLigand'] = 'Choose...(AD3)'
menuText['AdjustLigand'] = 'Ligand Parameters...(AD3)'

menuText['SetSearchParmsMB'] = 'Search Parameters'
menuText['SA'] = 'Simulated Annealing...'
menuText['GA'] = 'Genetic Algorithm...'
menuText['LS'] = 'Local Search...'

menuText['SetDockingRunParmsMB'] = 'Docking Parameters...'
menuText['OtherOptionsMB'] = 'Other Options...'
menuText['SetAutoDock4Parameters'] = 'AutoDock4 Parameters'
menuText['SetAutoDock41Parameters'] = 'AutoDock4.2 Parameters'

menuText['WriteDpfMB'] = 'Output'
menuText['WriteSA4'] = 'Simulated Annealing...'
menuText['WriteGA4'] = 'Genetic Algorithm...'
menuText['WriteLS4'] = 'Local Search...'
menuText['WriteGALS4'] = 'Lamarckian GA...'
menuText['WriteSA41'] = 'Simulated Annealing(4.2)...'
menuText['WriteGA41'] = 'Genetic Algorithm(4.2)...'
menuText['WriteLS41'] = 'Local Search(4.2)...'
menuText['WriteGALS41'] = 'Lamarckian GA(4.2)...'
menuText['WriteEPDB41'] = 'Evaluate Energy(EPDB)...'
menuText['WriteCONFIG41'] = 'Vina Config(config.txt)...'
menuText['WriteCluster4'] = 'Clustering...'
menuText['WriteSA'] = 'Simulated Annealing...(AD3)'
menuText['WriteGA'] = 'Genetic Algorithm...(AD3)'
menuText['WriteLS'] = 'Local Search...(AD3)'
menuText['WriteGALS'] = 'Lamarckian GA...(AD3)'
menuText['WriteCluster'] = 'Clustering...(AD3)'

menuText['EditDpfMB'] = 'Edit DPF...'



def checkHasDpo(vf):
    if not hasattr(vf, 'dpo'):
        vf.dpo = DockingParameters()
        #the dpo needs to know its vf
        vf.dpo.vf = vf
    if not hasattr(vf, 'vinaDict'):
        vf.vinaDict = {}



class DpfSetDpo(MVCommand):


    def onAddCmdToViewer(self):
        checkHasDpo(self.vf)
        if not hasattr(self.vf, 'vinaDict'):
            self.vf.vinaDict = {}


    def doit(self, *args, **kw):
        if not len(list(kw.items())):
            return 'ERROR'
        for key, val in list(kw.items()):
            self.vf.dpo[key]['value'] = val     
            if key in ['ligand','receptor','flexres','flex','center','out','size']:
                self.vf.vinaDict[key]['value'] = val



class DpfLoadDefaults(MVCommand):
    """ allows user to select a file containing a set of defaults"""


    def onAddCmdToViewer(self):
        checkHasDpo(self.vf)


    def __call__(self, dpffile, **kw):
        """None<-ADdpf_read
        dpffile is name of file whose contents are to be used to set values in the dpo
"""
        #check that dpffile actually exists
        if not os.path.exists(dpffile):
            raise IOError
        self.doitWrapper(*(dpffile,), **kw)


    def doit(self, dpffile):
        self.vf.dpo.read(dpffile)
                

    def guiCallback(self):
        """called each time the 'select defaultList file' button is pressed"""
        dpfFile = self.vf.askFileOpen(types=[('select default filename:', '*.dpf')],
                title = 'DPF Default File:')
        if dpfFile:
            self.doitWrapper(dpfFile, log=1, redraw=0)


DpfLoadDefaultsGUI = CommandGUI()
DpfLoadDefaultsGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'], menuText['ReadDpfMB'])



class DpfMacroSelector(MVCommand):
    """ allows user to select a filename for the macromolecule"""


    def onAddCmdToViewer(self):
        checkHasDpo(self.vf)
        if not hasattr(self.vf, 'vinaDict'): self.vf.vinaDict={} 


    def guiCallback(self):
        """called each time the 'select pdb macromolecule' button is pressed"""
        macroFile = self.vf.askFileOpen(types=[('PDBQS files', '*.pdbqs')],
                title = 'PDBQS Macromolecule File:')
        if macroFile:
            filename=os.path.basename(macroFile)
            self.doitWrapper(macroFile,log=1,redraw=0)


    def __call__ (self, macroFile, **kw):
        """None<-ADdpf_readMacromolecule
        macroFile file containing the receptor
"""
        if not os.path.exists(macroFile):
            raise IOError
        self.doitWrapper(*(macroFile,), **kw)


    def doit(self, macroFile):

        filename=os.path.basename(macroFile)

        ftype = string.split(filename, '.')[-1]
        if ftype != 'pdbqs':
            msgStr="macromolecule must be in pdbqs format!"
            self.vf.warningMsg(msgStr)
            return 'ERROR'

        rnum=string.rfind(filename, '.')
        if rnum<0:
            t= "illegal filename "+ filename
            self.vf.warningMsg(t)
            return
        #setting dpo.molstem->setting 'fld' and 'map'
        self.vf.dpo.molstem = filename[:rnum]
        #call set_receptor to set 'fld'
        self.vf.dpo.set_receptor(filename)



DpfMacroSelectorGUI=CommandGUI()
DpfMacroSelectorGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'],\
menuText['ReadMacro'], cascadeName = menuText['MacromoleculeMB'])



class Dpf4MacroSelector(MVCommand):
    """ allows user to select a filename for the macromolecule"""


    def onAddCmdToViewer(self):
        checkHasDpo(self.vf)
        if not hasattr(self.vf, 'vinaDict'): self.vf.vinaDict={}


    def onRemoveObjFromViewer(self, obj):
        if hasattr(self.vf.dpo, 'molstem') and self.vf.dpo.molstem==obj.name:
            delattr(self.vf.dpo, 'molstem')
            

    def guiCallback(self):
        """called each time the 'select pdbqt macromolecule filename' button is pressed"""
        macroFile = self.vf.askFileOpen(types=[('PDBQT files', '*.pdbqt')],
                title = 'PDBQT Macromolecule File:')
        if macroFile:
            filename=os.path.basename(macroFile)
            self.doitWrapper(macroFile,log=1,redraw=0)


    def __call__ (self, macroFile, **kw):
        """None<-ADdpf4_readMacromolecule
        macroFile file containing the receptor
"""
        if not os.path.exists(macroFile):
            raise IOError
        self.doitWrapper(*(macroFile,), **kw)


    def doit(self, macroFile):
        filename=os.path.basename(macroFile)
        molstem, ftype = os.path.splitext(filename)
        if ftype != '.pdbqt':
            msgStr="macromolecule must be in pdbqt format!"
            self.vf.warningMsg(msgStr)
            return 'ERROR'
        self.vf.dpo.molstem = molstem
        if hasattr(self.vf, 'flexDict') and \
                'rigid_filename' in self.vf.flexDict and \
                self.vf.flexDict['rigid_filename']==macroFile:
            filename = self.vf.flexDict['rigid_filename']
        self.vf.dpo.set_receptor(filename)


Dpf4MacroSelectorGUI=CommandGUI()
Dpf4MacroSelectorGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'],\
menuText['ReadMacro4'], cascadeName = menuText['MacromoleculeMB'])



class Dpf4FlexResSelector(MVCommand):
    """ allows user to select a flexible residue filename"""


    def onAddCmdToViewer(self):
        self.vf.DPF_FLEXRES_TYPES = []
        checkHasDpo(self.vf)
        if not hasattr(self.vf, 'vinaDict'): self.vf.vinaDict={}


    def guiCallback(self):
        """called each time the 'select flexible residue filename' button is pressed"""
        macroFile = self.vf.askFileOpen(types=[('PDBQT files', '*.pdbqt')],
                title = 'PDBQT Macromolecule Flexible Residue File:')
        if macroFile:
            filename=os.path.basename(macroFile)
            self.doitWrapper(macroFile,log=1,redraw=0)


    def __call__ (self, macroFile, **kw):
        """None<-ADdpf4_setMacroFlexibleResidueFilename
        macroFile file containing the receptor
"""
        if not os.path.exists(macroFile):
            raise IOError
        self.doitWrapper(*(macroFile,), **kw)


    def doit(self, macroFile):
        filename = os.path.basename(macroFile)
        molstem, ftype = os.path.splitext(filename)
        if ftype != '.pdbqt':
            msgStr="macromolecule must be in pdbqt format!"
            self.vf.warningMsg(msgStr)
            return 'ERROR'

        if hasattr(self.vf, 'flexDict') and \
                'flex_filename' in self.vf.flexDict and \
                self.vf.flexDict['flex_filename']==macroFile:
            filename = os.path.basename(self.vf.flexDict['flex_filename'])
        self.vf.dpo['flexres']['value'] = filename
        self.vf.dpo['flexres_flag']['value'] = True
        # check for any new 'ligand_types' in flexres  file
        self.vf.dpo.flexres = Read(macroFile)[0]
        d = {}
        for a in self.vf.dpo.flexres.allAtoms:
            d[a.autodock_element] = 1
        flex_types = list(d.keys())
        self.vf.DPF_FLEXRES_TYPES = flex_types
        self.vf.dpo['ndihe']['value'] += self.vf.dpo.flexres.parser.keys.count('BRANCH')
        self.vf.dpo.flexres.ndihe = self.vf.dpo.flexres.parser.keys.count('BRANCH')
        self.vf.vinaDict['flexres'] = filename


Dpf4FlexResSelectorGUI=CommandGUI()
Dpf4FlexResSelectorGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'],\
menuText['ReadFlexRes4'], cascadeName = menuText['MacromoleculeMB'])



class DpfMacroChooser(MVCommand):
    """ allows user to choose a molecule already present for the macromolecule"""
    

    def onAddCmdToViewer(self):
        checkHasDpo(self.vf)
        if not hasattr(self.vf, 'vinaDict'): self.vf.vinaDict={}

    def onRemoveObjFromViewer(self, obj):
        if hasattr(self.vf.dpo, 'molstem') and self.vf.dpo.molstem==obj.name:
            delattr(self.vf.dpo, 'molstem')
            

    def __init__(self, mode='single', title = 'Choose Macromolecule'):
        MVCommand.__init__(self)
        self.mode = mode
        self.title = title


    def chooseMolecule_cb(self, event = None):
        """called each time the 'choose Molecule' button is pressed"""
        mols = self.chooser.getMolSet()
        if mols: self.doitWrapper(mols,redraw=0)
        self.chooser.form.withdraw()


    def guiCallback(self):
        self.chooser = MoleculeChooser(self.vf, self.mode, self.title)
        self.chooser.ipf.append({'name':'Select Button',
                                 'widgetType':tkinter.Button,
                                 'text':'Select Molecule',
                                 'wcfg':{'bd':6},
                                 'gridcfg':{'sticky':tkinter.E+tkinter.W},
                                 'command': self.chooseMolecule_cb})
        self.form = self.chooser.go(modal=0, blocking=0)
        lb = self.chooser.ipf.entryByName['Molecule']['widget'].lb
        lb.bind("<Double-Button-1>",self.chooseMolecule_cb)


    def __call__(self, nodes, **kw):
        self.doitWrapper(*(nodes,), **kw)


    def doit(self, nodes):
        macro = self.vf.expandNodes(nodes)[0]
        errCharge, resList = checkMolCharges(macro, self.vf)
        if hasattr(macro, 'outputfilename'):
            filename = macro.outputfilename
        else:
            filename = os.path.basename(macro.parser.filename)
        ext = os.path.splitext(filename)[1]
        if ext!='.pdbqs':
            msg = filename + " must be a pdbqs file for autodock3"
            self.vf.warningMsg(msg)
            return "ERROR"
        errCharge, resList = checkMolCharges(macro, self.vf)
        #setting dpo.molstem->setting 'fld' and 'map'
        self.vf.dpo.molstem = os.path.splitext(filename)[0]
        #call set_receptor to set 'fld'
        self.vf.dpo.set_receptor(filename)


DpfMacroChooserGUI=CommandGUI()
DpfMacroChooserGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'],\
menuText['ChooseMacro'], cascadeName = menuText['MacromoleculeMB'], separatorAbove=1)



class Dpf4MacroChooser(MVCommand):
    """ allows user to choose a molecule already present for the macromolecule"""
    

    def onAddCmdToViewer(self):
        #self.vf.loadModule('editCommands', 'Pmv')
        checkHasDpo(self.vf)
        if not hasattr(self.vf, 'vinaDict'): self.vf.vinaDict={}


    def __init__(self, mode='single', title = 'Choose Macromolecule'):
        MVCommand.__init__(self)
        self.mode = mode
        self.title = title


    def chooseMolecule_cb(self, event = None):
        """called each time the 'choose Molecule' button is pressed"""
        mols = self.chooser.getMolSet()
        if mols: self.doitWrapper(mols,log=1,redraw=0)
        self.chooser.form.withdraw()


    def guiCallback(self):
        self.chooser = MoleculeChooser(self.vf, self.mode, self.title)
        self.chooser.ipf.append({'name':'Select Button',
                                 'widgetType':tkinter.Button,
                                 'text':'Select Molecule',
                                 'wcfg':{'bd':6},
                                 'gridcfg':{'sticky':tkinter.E+tkinter.W},
                                 'command': self.chooseMolecule_cb})
        self.form = self.chooser.go(modal=0, blocking=0)
        lb = self.chooser.ipf.entryByName['Molecule']['widget'].lb
        lb.bind("<Double-Button-1>",self.chooseMolecule_cb)


    def __call__(self, nodes, **kw):
        self.doitWrapper(*(nodes,), **kw)


    def doit(self, nodes):
        macro = self.vf.expandNodes(nodes)[0]
        if hasattr(self.vf, 'flexDict') and self.vf.flexDict['macroname']==macro.name \
                and 'rigid_filename' in self.vf.flexDict:
            filename = self.vf.flexDict['rigid_filename']
            errCharge, resList = self.vf.checkResCharges(residues)
            self.vf.dpo['receptor_types']['value'] = join(self.vf.flexDict['rigid_types'])
        else:
            filename = os.path.basename(macro.parser.filename)
            errCharge, resList = checkMolCharges(macro, self.vf)
        self.vf.dpo.molstem, ext = os.path.splitext(filename)
        if ext!='.pdbqt':
            msg = filename + " must be a pdbqt file for autodock4"
            self.vf.warningMsg(msg)
            return "ERROR"
        #call set_receptor to set 'fld'
        self.vf.dpo.set_receptor(filename)
        print("set receptor_filename to ", self.vf.dpo.receptor_filename)


Dpf4MacroChooserGUI=CommandGUI()
Dpf4MacroChooserGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'],\
menuText['ChooseMacro4'], cascadeName = menuText['MacromoleculeMB'])



class DpfInitLigand(MVCommand):
    """initializes the ligand"""


    def onAddCmdToViewer(self):
        checkHasDpo(self.vf)
        #self.vf.loadModule('editCommands', 'Pmv')
        #remove torsion constraint stuff?
        self.torsNum=[]
        self.halfWidth=[]
        self.hardTorCts=[]
        self.gaussTorCts=[]
        self.hTCtsnum=0
        self.torsPAngle=[]
        if self.vf.hasGui:
            self.about= tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.about.set('')
            self.types = tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.torsdof= 0
            self.torsdofcoeff = tkinter.StringVar(master=self.vf.GUI.ROOT)

            self.barrier = tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.barVar = tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.specifyTorCts = tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.showtorpen = tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.addTorCts = tkinter.IntVar(master=self.vf.GUI.ROOT)
            #counter for gaussian torsion constraints
            self.tcts = tkinter.IntVar(master=self.vf.GUI.ROOT)
            self.fmap = tkinter.IntVar(master=self.vf.GUI.ROOT)
            self.initTransType =tkinter.IntVar(master=self.vf.GUI.ROOT)
            self.initQuatType = tkinter.IntVar(master=self.vf.GUI.ROOT)
            self.specifyDihe = tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.specifyDihe.set('1')
            self.initDiheType = tkinter.IntVar(master=self.vf.GUI.ROOT)
            self.tran0 = tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.quat0 = tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.dihe0 = tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.ligMsgStr = tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.ligMsgStr.set('Ligand: ')
            self.typesMsgStr = tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.typesMsgStr.set('Ligand Atom Types: ' + self.types.get())
            self.tdofMsgStr = tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.tdofMsgStr.set('No Torsional Degrees of Freedom\nSpecified in Ligand!')
            self.centerMsgStr = tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.centerMsgStr.set('Center of Ligand Molecule:  ' + self.about.get())
            self.ndiheMsgStr=tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.ndiheMsgStr.set('No Torsions Specified in Ligand!')


    def onRemoveObjectFromViewer(self, obj):
        if hasattr(self.vf.dpo, 'ligand') and self.vf.dpo.ligand==obj:
            self.vf.dpo.ligand = None
            delattr(self.vf.dpo, 'ligand')
            #close the form if it is open
            try:
                self.Close_cb()
            except:
                pass
        if hasattr(self, 'old_ligand') and obj==self.old_ligand:
            self.old_ligand = None
            delattr(self, 'old_ligand')


    def buildForm(self):
        #the newstuff:
        ifd = self.ifd = InputFormDescr(title = "AutoDpf Ligand Parameters")
        ifd.append( {'name': 'ligLab',
            'widgetType':tkinter.Label,
            'textvariable': self.ligMsgStr,
            'gridcfg':{'sticky':tkinter.W,'columnspan':3}})
        ifd.append( {'name': 'typesLab',
            'widgetType':tkinter.Label,
            'label': 'Ligand Atom Types:',
            'textvariable': self.typesMsgStr,
            'gridcfg':{'sticky':tkinter.W,'columnspan':6}})
        ifd.append( {'name': 'centerLab',
            'widgetType':tkinter.Label,
            'label': 'Center of Ligand Molecule:',
            'textvariable': self.centerMsgStr,
            'gridcfg':{'sticky':tkinter.W,'columnspan':3}})
        ifd.append( {'name': 'initLab',
            'widgetType':tkinter.Label,
            'text': 'Set Initial State of Ligand:',
            'gridcfg':{'sticky':tkinter.W,'columnspan':3}})
        ifd.append( {'name':'tran0Ent',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'User-Specified Initial Position:',
                'textvariable': self.tran0
            },
            'gridcfg':{'sticky':tkinter.E,'columnspan':3}})
        ifd.append({'name': 'initTransCB',
            'widgetType':tkinter.Checkbutton,
            'text': 'Random ',
            'variable': self.initTransType,
            'command': self.set_initTransType,
            'gridcfg':{'sticky':tkinter.E, 'row':-1,'column':3}})
        ifd.append( {'name': 'quat0Ent',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Initial Relative Dihedral Offset(quat0):',
                'textvariable': self.quat0
            },
            'gridcfg':{'sticky':tkinter.W,'columnspan':3}})
        ifd.append({'name': 'initQuatCB',
            'widgetType':tkinter.Checkbutton,
            'text': 'Random ',
            'variable': self.initQuatType,
            'command': self.set_initQuatType,
            'gridcfg':{'sticky':tkinter.E, 'row':-1,'column':3}})
        ifd.append( {'name': 'ndiheLab',
            'widgetType':tkinter.Label,
            'textvariable': self.ndiheMsgStr,
            'gridcfg':{'sticky':tkinter.W,'columnspan':3}})
        ifd.append( {'name': 'tdofLab',
            'widgetType':tkinter.Label,
            'textvariable': self.tdofMsgStr,
            'gridcfg':{'sticky':tkinter.W,'columnspan':3}})
        ifd.append( {'name': 'torsdofcoeffEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Coefficient of torsional degrees of freedom:',
                'textvariable': self.torsdofcoeff
            },
            'gridcfg':{'sticky':tkinter.W,'columnspan':6}})
        ifd.append( {'name': 'diheChoic3Lab',
            'widgetType':tkinter.Label,
            'text': 'Specify Initial Dihedrals?',
            'gridcfg':{'sticky':tkinter.W, 'columnspan':2}})
        ifd.append({'name': 'yesDiheRB',
            'widgetType':tkinter.Radiobutton,
            'wcfg': {'value':'1'},
            'text': 'Yes',
            'variable': self.specifyDihe,
            'command': self.getDihe,
            'gridcfg':{'sticky':tkinter.W,'row':-1, 'column':3}})
        ifd.append({'name': 'noDiheRB',
            'widgetType':tkinter.Radiobutton,
            'wcfg': {'value':'0'},
            'text': 'No',
            'command': self.getDihe,
            'variable': self.specifyDihe,
            'gridcfg':{'sticky':tkinter.W, 'row':-1, 'column':4}})
        ifd.append( {'name': 'useInitDiheLab',
            'widgetType':tkinter.Label,
            'text': 'User Specified\nInitial Relative Dihedrals:',
            'gridcfg':{'sticky':tkinter.W}})
        ifd.append( {'name': 'dihe0Ent',
            'widgetType':tkinter.Entry,
            'wcfg':{ 'textvariable': self.dihe0 },
            'gridcfg':{'sticky':tkinter.E+tkinter.W,'row':-1, 'column':1, 'columnspan':3}})
        ifd.append({'name': 'initDiheCB',
            'widgetType':tkinter.Checkbutton,
            'text': 'Random ',
            'variable': self.initDiheType,
            'command': self.set_initDiheType,
            'gridcfg':{'sticky':tkinter.E, 'row':-1,'column':4}})
        #ifd.append( {'name': 'torCtsChoic3Lab',
            #'widgetType':Tkinter.Label,
            #'text': 'Specify Ligand Torsion Constraints?',
            #'gridcfg':{'sticky':Tkinter.W, 'columnspan':2}})
        #ifd.append({'name': 'yestorCtsRB',
            #'widgetType':Tkinter.Radiobutton,
            #'wcfg': {'value':'1'},
            #'text': 'Yes',
            #'variable': self.specifyTorCts,
            #'command': self.getTorCts,
            #'gridcfg':{'sticky':Tkinter.W,'row':-1, 'column':3}})
        #ifd.append({'name': 'noTorCtsRB',
            #'widgetType':Tkinter.Radiobutton,
            #'wcfg': {'value':'0'},
            #'text': 'No',
            #'command': self.getTorCts,
            #'variable': self.specifyTorCts,
            #'gridcfg':{'sticky':Tkinter.W, 'row':-1, 'column':4}})
        #ifd.append( {'name': 'barrierLab',
            #'widgetType':Tkinter.Label,
            #'text': 'Height of Energy-barrier applied\nto Constrained Torsions:',
            #'gridcfg':{'sticky':Tkinter.W}})
        #ifd.append( {'name': 'barrierEnt',
            #'widgetType':Tkinter.Entry,
            #'textvariable': self.barrier,
            #'gridcfg':{'sticky':Tkinter.E+Tkinter.W,'row':-1, 'column':1, 'columnspan':3}})
        #ifd.append({'name': 'barrierCB',
            #'widgetType':Tkinter.Checkbutton,
            #'wcfg':{ 'onvalue':'1', 'offvalue':'0'},
            #'text': 'None ',
            #'variable': self.barVar,
            #'command': self.set_barrier,
            #'gridcfg':{'sticky':Tkinter.E, 'row':-1,'column':4}})
        #ifd.append( {'name': 'addGCtsBut',
            #'widgetType':Tkinter.Button,
            #'text': 'Add Gaussian Constraint:\n(>1 possible per torsion)',
            #'command': self.getNewGaussTorCts,
            #'gridcfg':{'sticky':Tkinter.W+Tkinter.E, 'columnspan':2}})
        #ifd.append( {'name': 'addHCtsBut',
            #'widgetType':Tkinter.Button,
            #'text': 'Add A Hard Constraint:\n(only 1 possible per torsion)',
            #'command': self.getNewHardTorCts,
            #'gridcfg':{'sticky':Tkinter.W+Tkinter.E,'row':-1, 'column':2, 'columnspan':4}})
        ifd.append({'name': 'acceptB',
            'widgetType': tkinter.Button,
            'text':'Accept',
            'wcfg':{'bd':4},
            'gridcfg':{'sticky':tkinter.E+tkinter.W, 'columnspan':3},
            'command':self.Accept_cb})
        ifd.append({'name': 'closeB',
            'widgetType': tkinter.Button,
            'text':'Close',
            'wcfg':{'bd':4},
            'gridcfg':{'sticky':tkinter.E+tkinter.W, 'row':-1, 'column':3,'columnspan':3},
            'command':self.Close_cb})

        #initialize tran0, dihe0 and quat0
        self.tran0.set(self.about.get())
        nstr= '0. '*self.ndihe
        self.dihe0.set(nstr)
        self.quat0.set('1.0 0. 0. 0.')
        
    

    def Accept_cb(self, event=None):
        changeVals = {}
        for item in [ 'dihe0Ent', 'tran0Ent', 'quat0Ent']:
            var = self.ifd.entryByName[item]['wcfg']['textvariable']
            #FIX THIS
            if self.vf.dpo[item[:-3]]['value']!= var.get():
                changeVals[item[:-3]] =  var.get()
        oldVal = self.vf.dpo['torsdof']['value']
        newTDC =  float(self.torsdofcoeff.get())
        if oldVal[0]!= self.torsdof or oldVal[1]!= newTDC:
                changeVals['torsdof'] =  [self.torsdof, newTDC]
        if len(list(changeVals.keys()))>0:
            changeVals['topCommand'] = 0
            self.doitWrapper(*(), **changeVals)
        self.form.withdraw()


    def Close_cb(self, event=None):
        self.form.withdraw()


    def doit(self,*args,  **kw):
        self.vf.ADdpf_setDpo(*(), **kw)
        ##3/5:
        #DO SOME LIGAND CHECKING HERE?
        #ligand = self.vf.dpo.ligand
        #errCharge, resList = checkMolCharges(ligand, self.vf)
        #allAts = ligand.allAtoms
        #self.vf.colorByAtomType(ligand, ['lines'], topCommand = 0)
        #aromaticCs = AtomSet(filter(lambda x: x.autodock_element=='A',allAts))
        #if len(aromaticCs): 
            #self.vf.color(aromaticCs,((0.,1.,0.,1.),),['lines'],topCommand=0, redraw=1)


    def getNewGaussTorCts(self, event=None):
        num = self.tcts.get()
        self.torsNum.append(tkinter.StringVar(master=self.vf.GUI.ROOT))
        self.torsPAngle.append(tkinter.StringVar(master=self.vf.GUI.ROOT))
        self.halfWidth.append(tkinter.StringVar(master=self.vf.GUI.ROOT))
        ifd2 = self.ifd2 = InputFormDescr(title = "Gaussian Torsion Constraint Parameters")
        numStr = 'Torsion Number:\n(1-'+str(self.ndihe)+')'
        self.ifd2.append( { 'name': 'tnumEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': numStr,
                'textvariable': self.torsNum[num]
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':2}})
        self.ifd2.append( { 'name': 'pangEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Perferred Angle(degrees):',
                'textvariable': self.torsPAngle[num]
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':2}})
        self.ifd2.append( { 'name': 'hwidthEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Half-Width(degrees):',
                'textvariable': self.halfWidth[num]
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':2}})
        self.ifd2.append({'name': 'showTorpenCB',
            'widgetType':tkinter.Radiobutton,
            'text': 'Store + Output torsion energies',
            'wcfg': {'value':'1'},
            'variable': self.showtorpen,
            'gridcfg':{'sticky':tkinter.W}})
        self.ifd2.append({'name': 'noTorpenCB',
            'widgetType':tkinter.Radiobutton,
            'text': 'Don\'t Store + Output torsion energies',
            'variable': self.showtorpen,
            'wcfg': {'value':'0'},
            'gridcfg':{'sticky':tkinter.W}})
        val = self.vf.getUserInput(self.ifd2)
        if len(val) and val['tnumEnt']=='' or val['pangEnt']=='' or val['hwidthEnt']=='':
            val = None
        if not val:
            del self.torsNum[num]
            del self.torsPAngle[num]
            del self.halfWidth[num]
        else:
            newStr=[]
            newStr.append(self.torsNum[num].get())
            newStr.append(self.torsPAngle[num].get())
            newStr.append(self.halfWidth[num].get())
            self.gaussTorCts.append(newStr)
            self.tcts.set(num+1)
            
        
    def getNewHardTorCts(self, event=None):
        torsNum =tkinter.StringVar(master=self.vf.GUI.ROOT)
        torsPAngle=tkinter.StringVar(master=self.vf.GUI.ROOT)
        halfWidth=tkinter.StringVar(master=self.vf.GUI.ROOT)
        ifd2 = self.ifd2 = InputFormDescr(title = "Hard Torsion Constraint Parameters")
        numStr = 'Torsion Number:\n(1-'+str(self.ndihe)+')'
        self.ifd2.append( { 'name': 'tnumEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': numStr,
                'textvariable': torsNum
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':2}})
        self.ifd2.append( { 'name': 'pangEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Perferred Relative Angle(degrees):',
                'textvariable': torsPAngle
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':2}})
        self.ifd2.append( { 'name': 'hwidthEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Full Width of allowed range(degrees):',
                'textvariable': halfWidth
            },
            'gridcfg':{'sticky':tkinter.W, 'columnspan':4}})
        val = self.vf.getUserInput(self.ifd2)
        if val:
            num = int(val['tnumEnt'])
            newEnt=[]
            newEnt.append(val['tnumEnt'])
            newEnt.append(val['pangEnt'])
            newEnt.append(val['hwidthEnt'])
            self.hardTorCts[num]=newEnt


    def getTorCts(self, event=None):
        w=self.ifd.entryByName['barrierLab']
        w1=self.ifd.entryByName['barrierEnt']
        w2=self.ifd.entryByName['barrierCB']
        w3=self.ifd.entryByName['addGCtsBut']
        w4=self.ifd.entryByName['addHCtsBut']
        if self.specifyTorCts.get()=='1':
            w['widget'].grid(w['gridcfg'])
            w1['widget'].grid(w1['gridcfg'])
            w2['widget'].grid(w2['gridcfg'])
            w3['widget'].grid(w3['gridcfg'])
            w4['widget'].grid(w4['gridcfg'])
        else:
            w['widget'].grid_forget()
            w1['widget'].grid_forget()
            w2['widget'].grid_forget()
            w3['widget'].grid_forget()
            w4['widget'].grid_forget()
                

    def getDihe(self, event=None):
        w=self.ifd.entryByName['useInitDiheLab']
        w1=self.ifd.entryByName['dihe0Ent']
        w2=self.ifd.entryByName['initDiheCB']
        if self.specifyDihe.get()=='1':
            w['widget'].grid(w['gridcfg'])    
            w1['widget'].grid(w1['gridcfg'])    
            w2['widget'].grid(w2['gridcfg'])    
        else:
            w['widget'].grid_forget()
            w1['widget'].grid_forget()
            w2['widget'].grid_forget()


    def set_barrier(self):
        if self.barVar.get()=='1':
            self.barrier.set('')


    def set_initDiheType(self):
        if self.initDiheType.get()==1:
            self.dihe0.set('random')
        else:
            nstr= '0. '*self.ndihe
            self.dihe0.set(nstr)


    def set_initTransType(self):
        w=self.ifd.entryByName['tran0Ent']['widget']
        if self.initTransType.get()==1:
            self.tran0.set('random')
        else:
            self.tran0.set(self.about.get())

        
    def set_initQuatType(self):
        if self.initQuatType.get()==1:
            self.quat0.set('random')
        else:
            self.quat0.set('1.0 0. 0. 0.')


    def guiCallback(self):
        ligandfilename =self.vf.dpo['move']['value']
        ligand = self.vf.dpo.ligand
        ligstring = ligand.types
        if hasattr(ligand, 'ndihe'):
            self.ndihe = ligand.ndihe
        elif hasattr(ligand, 'torscount'):
            self.ndihe = ligand.torscount
        else:
            msgStr="Selected Ligand not formatted!!"
            self.vf.warningMsg(msgStr)
            return
        self.torsdof = ligand.TORSDOF
        self.types.set(ligand.types)
        if ligand==None:
            msgStr="No Ligand Selected!"
            self.vf.warningMsg(msgStr)
            return
        if not hasattr(self, 'form'):
            self.buildForm()
            self.form = self.vf.getUserInput(self.ifd, modal=0, blocking=0)
            self.form.root.protocol('WM_DELETE_WINDOW',self.Close_cb)
            #self.barVar.set('0')
            self.initTransType.set(1)
            self.initDiheType.set(1)
            self.initQuatType.set(1)
            self.fmap.set(0)
            #this overwrites previous values
            self.set_initTransType()
            self.set_initQuatType()
            self.set_initDiheType()
        else:
            self.form.root.deiconify()
        self.torsdofcoeff.set(str(self.vf.dpo['torsdof']['value'][1]))
        self.vf.dpo['torsdof']['value'][0] = self.torsdof
        #self.fmap.set(self.vf.dpo['fmap']['value'])
        #self.tcts.set(self.vf.dpfDict['tcts'])
        #self.barrier.set(self.vf.dpfDict['barrier'])
        #self.getTorCts()
        #self.specifyDihe.set(self.vf.dpfDict['specifyDihe'])
        #self.getDihe()
        #self.specifyTorCts.set(str(self.vf.dpfDict['specifyTorCts']))
        #self.showtorpen.set(str(self.vf.dpfDict['showtorpen']))
        if not hasattr(self, 'old_ligand') or self.old_ligand!=ligand:
            v = self.vf.dpo['about']['value']
            strCenter = "%6.3f %6.3f %6.3f" %(v[0], v[1], v[2])
            self.about.set(strCenter)
            self.old_ligand=ligand
            self.ligMsgStr.set( 'Ligand:   '+ ligandfilename)
            if self.ndihe:
                self.ndiheMsgStr.set('Number of Active Torsions in Ligand:   '+ str(self.ndihe))
                #self.hardTorCts has one possible entry 
                #per active torsion,kept in hardTorCts list
                #for i in range(int(self.ndihe)):
                    #self.hardTorCts.append([])
            else:
                self.ndiheMsgStr.set('No Torsions Specified in Ligand!')
            if self.torsdof:
                self.tdofMsgStr.set('Number of Torsional Degrees of Freedom(torsdof) in Ligand:   '+ str(self.torsdof))
            else:
                self.tdofMsgStr.set('No Torsional Degrees of\nFreedom in Ligand!  ')
            self.typesMsgStr.set('Ligand Atom Types: ' + self.types.get())
            self.centerMsgStr.set('Center of Ligand Molecule:  ' + self.about.get())


DpfInitLigandGUI=CommandGUI()
DpfInitLigandGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'], menuText['AdjustLigand'], cascadeName = menuText['SetLigandParmsMB'])



sortKeyList =  ['C','A','N','O','S','H','P','n','f','F','c','b','I','M']
def setDpoFields(ligand, dpo):
    #first the strings:
    if hasattr(ligand, 'outputfile'):
        filename = os.path.basename(ligand.outputfile)
    else:
        filename = os.path.basename(ligand.parser.filename)
    dpo['move']['value'] = filename
    dpo['rmsref']['value'] =  filename
    #next a list:
    l = ligand.allAtoms.autodock_element
    sortedL = []
    for t in sortKeyList:
        if t in l: sortedL.append(t)
    ligand.types = string.join(sortedL,'')
    dpo['types']['value'] =  ligand.types

    #next a list of floats:
    if not hasattr(ligand, 'dpf_center'):
        #to find center of ligand atoms only:
        #   check whether there are BEGIN_RES/END_RES records
        if hasattr(ligand.parser, 'keys') and 'BEGIN_' in ligand.parser.keys:
            key_list = ligand.parser.keys
            ind = key_list.index('BEGIN_')
            k_list = key_list[:ind]
            atom_ct = k_list.count("ATOM") + k_list.count("HETATM")
            ats = ligand.allAtoms[:atom_ct]
            center = list(numpy.sum(ats.coords, 0)/(len(ats)*1.0))
            for i in range(3):
                center[i] = round(center[i], 4)
            ligand.center = center
        else:
            ligand.getCenter()
        ligand.dpf_center = ligand.center
    dpo['about']['value'] =  [round(ligand.center[0],3), round(ligand.center[1],3),\
            round(ligand.center[2],3)]

    #the last integers:
    dpo['ndihe']['value'] =  ligand.ndihe
    dpo['torsdof']['value'][0] =  ligand.TORSDOF
    dpo['torsdof4']['value'][0] =  ligand.ndihe
    dpo.ligand = ligand
    dpo.ligand_filename = filename



class DpfLigandChooser(MVCommand):
    """ allows user to choose a molecule already present for the ligand"""


    def onAddCmdToViewer(self):
        checkHasDpo(self.vf)
        #self.vf.loadModule('autotorsCommands', 'AutoDockTools')
        if not hasattr(self.vf, 'atorsDict'):
            self.vf.atorsDict = {}
        self.vf.loadModule('displayCommands')
        self.vf.loadModule('bondsCommands')
        self.vf.loadModule('fileCommands')


    def chooseLigand_cb(self, event = None):
        """called each time the 'choose Ligand' button is pressed"""
        mol = self.chooser.getMolSet()
        if mol: 
            ask = 1
            self.doitWrapper(mol, ask, log=1, redraw=1)
        self.chooser.form.withdraw()


    def guiCallback(self):
        self.mode = 'single'
        self.title= 'Choose Ligand'
        self.chooser = MoleculeChooser(self.vf, self.mode, self.title)
        self.chooser.ipf.append({'name':'Select Button',
                                 'widgetType':tkinter.Button,
                                 'text':'Select Ligand',
                                 'wcfg':{'bd':6},
                                 'gridcfg':{'sticky':tkinter.E+tkinter.W},
                                 'command': self.chooseLigand_cb})
        self.form = self.chooser.go(modal=0, blocking=0)
        lb = self.chooser.ipf.entryByName['Molecule']['widget'].lb
        lb.bind("<Double-Button-1>",self.chooseLigand_cb)


    def __call__(self, nodes, ask=1, **kw):
        self.doitWrapper(*(nodes, ask), **kw)


    def doit(self, nodes, ask):
        lig = self.vf.expandNodes(nodes)[0]
        #what about if >1 of these are true?
        #NB this will change with AutoDockLigand class is implemented

        isAtorsMol = 0
        isDpoLigand = 0

        hasTorTree = hasattr(lig, 'torTree')

        d = self.vf.atorsDict
        if 'molecule' in d and d['molecule']==lig:
            isAtorsMol = 1
            lig.ndihe = lig.torscount
            #lig.ndihe = d['torscount']
            #lig.TORSDOF = d['TORSDOF']

        if hasattr(self.vf, 'dpo'):
            obj = self.vf.dpo
            if  hasattr(obj, 'ligand') and obj.ligand == lig:
                isDpoLigand = 1

        if hasTorTree or isAtorsMol or isDpoLigand:
            setDpoFields(lig, self.vf.dpo)
        else:
            msgStr = 'can only selectLigand\n preformatted by autotors\
            (and Written to a file)\n for this option!'
            self.vf.warningMsg(msgStr)
            return 'ERROR'

        if self.vf.hasGui: 
            self.vf.colorByAtomType(lig, topCommand = 0, redraw=1)
            aromaticCs = AtomSet([x for x in lig.allAtoms if x.autodock_element=='A'])
            if len(aromaticCs):
                self.vf.color(aromaticCs,((0.,1.,0.),),['lines'],topCommand=0, redraw=1)
            if ask:
                self.vf.ADdpf_initLigand.guiCallback()
            else:
                d = {}
                if hasattr(lig, 'outputfile'):
                    d['move'] = lig.outputfile
                    d['torsdof'] = [lig.torsdof, self.vf.dpo['torsdof']['value'][1]]
                    d['about'] = lig.center
                    d['ndihe'] = lig.ndihe
                    d['types'] = lig.types
                    self.vf.ADdpf_setDpo(*(), **d)


DpfLigandChooserGUI=CommandGUI()
DpfLigandChooserGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'],\
menuText['ChooseLigand'], cascadeName = menuText['SetLigandParmsMB'], separatorAbove=1)



class DpfLigPDBQReader(MVCommand):
    """ allows user to choose a PDBQ file for the ligand"""

    def onAddCmdToViewer(self):
        if self.vf.hasGui and not hasattr(self.vf,'readMolecule'):
            self.vf.loadCommand('fileCommands', 'readMolecule', 'Pmv')
        checkHasDpo(self.vf)


    def guiCallback(self):
        """called each time the 'select ligand' button is pressed"""
        ligFile = self.vf.askFileOpen(types=[('Autotors-format convention:', '*.out.pdbq'),('PDBQ files:', '*.pdbq')],
                title = 'Formatted Ligand File:')
        if ligFile:
            self.doitWrapper(ligFile, log=1, redraw=1)


    def __call__(self, ligFile, ask = 1, **kw):
        """None<-ADdpf_readFormattedLigand
        ligFile: file containing the ligand 
        ask: flag for whether to update geometry 
"""
        if not os.path.exists(ligFile):
            raise IOError
        kw['ask'] = ask
        self.doitWrapper(*(ligFile,), **kw)


    def doit(self, ligFile, ask=1):
        ligand = self.vf.readMolecule(ligFile)[0]
        ligFile = os.path.basename(ligFile)
        #ligand = self.vf.readPDBQ(ligFile)[0]
        assert hasattr(ligand, 'TORSDOF'), 'ligand must be preformatted with autotors'
        setDpoFields(ligand, self.vf.dpo)
        #color and check charges here
        if self.vf.hasGui and hasattr(ligand,'ndihe') and ask:
            self.vf.colorByAtomType(ligand, topCommand = 0, redraw=1)
            aromaticCs = AtomSet([x for x in ligand.allAtoms if x.autodock_element=='A'])
            if len(aromaticCs):
                self.vf.color(aromaticCs,((0.,1.,0.),),['lines'],topCommand=0, redraw=1)
            self.vf.ADdpf_initLigand.guiCallback()


        
DpfLigPDBQReaderGUI = CommandGUI()
DpfLigPDBQReaderGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'], menuText['ReadLigand'],  cascadeName = menuText['SetLigandParmsMB'])

   

class Dpf4InitLigand(MVCommand):
    """initializes the ligand for AutoDock4  FIX THIS TO SET TYPES APPROPRIATELY"""


    def onAddCmdToViewer(self):
        checkHasDpo(self.vf)
        self.vf.DPF_LIGAND_TYPES = []
        #self.vf.loadModule('editCommands', 'Pmv')
        #remove torsion constraint stuff?
        self.torsNum=[]
        self.halfWidth=[]
        self.hardTorCts=[]
        self.gaussTorCts=[]
        self.hTCtsnum=0
        self.torsPAngle=[]
        if self.vf.hasGui:
            self.about= tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.about.set('')
            self.types = tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.torsdof = 0
            self.torsdofcoeff = tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.torsdofcoeff.set("")
            self.barrier = tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.barVar = tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.specifyTorCts = tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.showtorpen = tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.addTorCts = tkinter.IntVar(master=self.vf.GUI.ROOT)
            #counter for gaussian torsion constraints
            self.tcts = tkinter.IntVar(master=self.vf.GUI.ROOT)
            self.fmap = tkinter.IntVar(master=self.vf.GUI.ROOT)
            self.initTransType =tkinter.IntVar(master=self.vf.GUI.ROOT)
            self.initQuatType = tkinter.IntVar(master=self.vf.GUI.ROOT)
            self.specifyDihe = tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.specifyDihe.set('1')
            self.initDiheType = tkinter.IntVar(master=self.vf.GUI.ROOT)
            self.tran0 = tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.quat0 = tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.dihe0 = tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.ligMsgStr = tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.ligMsgStr.set('Ligand: ')
            self.typesMsgStr = tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.typesMsgStr.set('Ligand Atom Types: ' + self.types.get())
            self.tdofMsgStr = tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.tdofMsgStr.set('No Torsional Degrees of Freedom\nSpecified in Ligand!')
            self.centerMsgStr = tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.centerMsgStr.set('Center of Ligand Molecule:  ' + self.about.get())
            self.ndiheMsgStr=tkinter.StringVar(master=self.vf.GUI.ROOT)
            self.ndiheMsgStr.set('No Torsions Specified in Ligand!')


    def onRemoveObjectFromViewer(self, obj):
        if hasattr(self.vf.dpo, 'ligand') and self.vf.dpo.ligand==obj:
            self.vf.dpo.ligand = None
            delattr(self.vf.dpo, 'ligand')
            #close the form if it is open
            try:
                self.Close_cb()
            except:
                pass
        if hasattr(self, 'old_ligand') and obj==self.old_ligand:
            self.old_ligand = None
            delattr(self, 'old_ligand')


    def buildForm(self):
        #the newstuff:
        ifd = self.ifd = InputFormDescr(title = "AutoDpf4 Ligand Parameters")
        ifd.append( {'name': 'ligLab',
            'widgetType':tkinter.Label,
            'textvariable': self.ligMsgStr,
            'gridcfg':{'sticky':tkinter.W,'columnspan':3}})
        ifd.append( {'name': 'typesLab',
            'widgetType':tkinter.Label,
            'label': 'Ligand Atom Types:',
            'textvariable': self.typesMsgStr,
            'gridcfg':{'sticky':tkinter.W,'columnspan':6}})
        ifd.append( {'name': 'centerLab',
            'widgetType':tkinter.Label,
            'label': 'Center of Ligand Molecule:',
            'textvariable': self.centerMsgStr,
            'gridcfg':{'sticky':tkinter.W,'columnspan':3}})
        ifd.append( {'name': 'initLab',
            'widgetType':tkinter.Label,
            'text': 'Set Initial State of Ligand:',
            'gridcfg':{'sticky':tkinter.W,'columnspan':3}})
        ifd.append( {'name':'tran0Ent',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'User-Specified Initial Position:',
                'textvariable': self.tran0
            },
            'gridcfg':{'sticky':tkinter.E,'columnspan':3}})
        ifd.append({'name': 'initTransCB',
            'widgetType':tkinter.Checkbutton,
            'text': 'Random ',
            'variable': self.initTransType,
            'command': self.set_initTransType,
            'gridcfg':{'sticky':tkinter.E, 'row':-1,'column':3}})
        ifd.append( {'name': 'quat0Ent',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Initial Relative Dihedral Offset(quat0):',
                'textvariable': self.quat0
            },
            'gridcfg':{'sticky':tkinter.E,'columnspan':3}})
        ifd.append({'name': 'initQuatCB',
            'widgetType':tkinter.Checkbutton,
            'text': 'Random ',
            'variable': self.initQuatType,
            'command': self.set_initQuatType,
            'gridcfg':{'sticky':tkinter.E, 'row':-1,'column':3}})
        ifd.append( {'name': 'ndiheLab',
            'widgetType':tkinter.Label,
            'textvariable': self.ndiheMsgStr,
            'gridcfg':{'sticky':tkinter.W,'columnspan':3}})
        ifd.append( {'name': 'tdofLab',
            'widgetType':tkinter.Label,
            'textvariable': self.tdofMsgStr,
            'gridcfg':{'sticky':tkinter.W,'columnspan':3}})
        #ifd.append( {'name': 'torsdofcoeffEnt',
        #    'widgetType':Tkinter.Entry,
        #    'wcfg':{
        #        'label': 'Coefficient of torsional degrees of freedom:',
        #        'textvariable': self.torsdofcoeff
        #    },
        #    'gridcfg':{'sticky':Tkinter.W,'columnspan':6}})
        ifd.append( {'name': 'diheChoic3Lab',
            'widgetType':tkinter.Label,
            'text': 'Specify Initial Dihedrals?',
            'gridcfg':{'sticky':tkinter.W, 'columnspan':2}})
        ifd.append({'name': 'yesDiheRB',
            'widgetType':tkinter.Radiobutton,
            'wcfg': {'value':'1'},
            'text': 'Yes',
            'variable': self.specifyDihe,
            'command': self.getDihe,
            'gridcfg':{'sticky':tkinter.W,'row':-1, 'column':3}})
        ifd.append({'name': 'noDiheRB',
            'widgetType':tkinter.Radiobutton,
            'wcfg': {'value':'0'},
            'text': 'No',
            'command': self.getDihe,
            'variable': self.specifyDihe,
            'gridcfg':{'sticky':tkinter.W, 'row':-1, 'column':4}})
        ifd.append( {'name': 'useInitDiheLab',
            'widgetType':tkinter.Label,
            'text': 'User Specified\nInitial Relative Dihedrals:',
            'gridcfg':{'sticky':tkinter.W}})
        ifd.append( {'name': 'dihe0Ent',
            'widgetType':tkinter.Entry,
            'wcfg':{ 'textvariable': self.dihe0 },
            'gridcfg':{'sticky':tkinter.E+tkinter.W,'row':-1, 'column':1, 'columnspan':3}})
        ifd.append({'name': 'initDiheCB',
            'widgetType':tkinter.Checkbutton,
            'text': 'Random ',
            'variable': self.initDiheType,
            'command': self.set_initDiheType,
            'gridcfg':{'sticky':tkinter.E, 'row':-1,'column':4}})
        #ifd.append( {'name': 'torCtsChoic3Lab',
            #'widgetType':Tkinter.Label,
            #'text': 'Specify Ligand Torsion Constraints?',
            #'gridcfg':{'sticky':Tkinter.W, 'columnspan':2}})
        #ifd.append({'name': 'yestorCtsRB',
            #'widgetType':Tkinter.Radiobutton,
            #'wcfg': {'value':'1'},
            #'text': 'Yes',
            #'variable': self.specifyTorCts,
            #'command': self.getTorCts,
            #'gridcfg':{'sticky':Tkinter.W,'row':-1, 'column':3}})
        #ifd.append({'name': 'noTorCtsRB',
            #'widgetType':Tkinter.Radiobutton,
            #'wcfg': {'value':'0'},
            #'text': 'No',
            #'command': self.getTorCts,
            #'variable': self.specifyTorCts,
            #'gridcfg':{'sticky':Tkinter.W, 'row':-1, 'column':4}})
        #ifd.append( {'name': 'barrierLab',
            #'widgetType':Tkinter.Label,
            #'text': 'Height of Energy-barrier applied\nto Constrained Torsions:',
            #'gridcfg':{'sticky':Tkinter.W}})
        #ifd.append( {'name': 'barrierEnt',
            #'widgetType':Tkinter.Entry,
            #'textvariable': self.barrier,
            #'gridcfg':{'sticky':Tkinter.E+Tkinter.W,'row':-1, 'column':1, 'columnspan':3}})
        #ifd.append({'name': 'barrierCB',
            #'widgetType':Tkinter.Checkbutton,
            #'wcfg':{ 'onvalue':'1', 'offvalue':'0'},
            #'text': 'None ',
            #'variable': self.barVar,
            #'command': self.set_barrier,
            #'gridcfg':{'sticky':Tkinter.E, 'row':-1,'column':4}})
        #ifd.append( {'name': 'addGCtsBut',
            #'widgetType':Tkinter.Button,
            #'text': 'Add Gaussian Constraint:\n(>1 possible per torsion)',
            #'command': self.getNewGaussTorCts,
            #'gridcfg':{'sticky':Tkinter.W+Tkinter.E, 'columnspan':2}})
        #ifd.append( {'name': 'addHCtsBut',
            #'widgetType':Tkinter.Button,
            #'text': 'Add A Hard Constraint:\n(only 1 possible per torsion)',
            #'command': self.getNewHardTorCts,
            #'gridcfg':{'sticky':Tkinter.W+Tkinter.E,'row':-1, 'column':2, 'columnspan':4}})
        ifd.append({'name': 'acceptB',
            'widgetType': tkinter.Button,
            'text':'Accept',
            'wcfg':{'bd':4},
            'gridcfg':{'sticky':tkinter.E+tkinter.W, 'columnspan':3},
            'command':self.Accept_cb})
        ifd.append({'name': 'closeB',
            'widgetType': tkinter.Button,
            'text':'Close',
            'wcfg':{'bd':4},
            'gridcfg':{'sticky':tkinter.E+tkinter.W, 'row':-1, 'column':3,'columnspan':3},
            'command':self.Close_cb})

        #initialize tran0, dihe0 and quat0
        self.tran0.set(self.about.get())
        nstr= '0. '*self.ndihe
        self.dihe0.set(nstr)
        self.quat0.set('1.0 0. 0. 0.')
        
    

    def Accept_cb(self, event=None):
        changeVals = {}
        for item in [ 'dihe0Ent', 'tran0Ent', 'quat0Ent']:
            var = self.ifd.entryByName[item]['wcfg']['textvariable']
            #FIX THIS
            if self.vf.dpo[item[:-3]]['value']!= var.get():
                changeVals[item[:-3]] =  var.get()
        oldVal = self.vf.dpo['torsdof4']['value']
        if oldVal[0]!= self.torsdof:
                changeVals['torsdof4'] =  [self.torsdof]
        if len(list(changeVals.keys()))>0:
            changeVals['topCommand'] = 0
            self.doitWrapper(*(), **changeVals)
        else:
            self.form.withdraw()


    def Close_cb(self, event=None):
        self.form.withdraw()


    def doit(self,*args,  **kw):
        ligand = self.vf.dpo.ligand
        d = {}
        for a in ligand.allAtoms:
            d[a.autodock_element] = 1
        self.vf.DPF_LIGAND_TYPES = list(d.keys())
        self.vf.ADdpf_setDpo(*(), **kw)


    def getNewGaussTorCts(self, event=None):
        num = self.tcts.get()
        self.torsNum.append(tkinter.StringVar(master=self.vf.GUI.ROOT))
        self.torsPAngle.append(tkinter.StringVar(master=self.vf.GUI.ROOT))
        self.halfWidth.append(tkinter.StringVar(master=self.vf.GUI.ROOT))
        ifd2 = self.ifd2 = InputFormDescr(title = "Gaussian Torsion Constraint Parameters")
        numStr = 'Torsion Number:\n(1-'+str(self.ndihe)+')'
        self.ifd2.append( { 'name': 'tnumEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': numStr,
                'textvariable': self.torsNum[num]
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':2}})
        self.ifd2.append( { 'name': 'pangEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Perferred Angle(degrees):',
                'textvariable': self.torsPAngle[num]
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':2}})
        self.ifd2.append( { 'name': 'hwidthEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Half-Width(degrees):',
                'textvariable': self.halfWidth[num]
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':2}})
        self.ifd2.append({'name': 'showTorpenCB',
            'widgetType':tkinter.Radiobutton,
            'text': 'Store + Output torsion energies',
            'wcfg': {'value':'1'},
            'variable': self.showtorpen,
            'gridcfg':{'sticky':tkinter.W}})
        self.ifd2.append({'name': 'noTorpenCB',
            'widgetType':tkinter.Radiobutton,
            'text': 'Don\'t Store + Output torsion energies',
            'variable': self.showtorpen,
            'wcfg': {'value':'0'},
            'gridcfg':{'sticky':tkinter.W}})
        val = self.vf.getUserInput(self.ifd2)
        if len(val) and val['tnumEnt']=='' or val['pangEnt']=='' or val['hwidthEnt']=='':
            val = None
        if not val:
            del self.torsNum[num]
            del self.torsPAngle[num]
            del self.halfWidth[num]
        else:
            newStr=[]
            newStr.append(self.torsNum[num].get())
            newStr.append(self.torsPAngle[num].get())
            newStr.append(self.halfWidth[num].get())
            self.gaussTorCts.append(newStr)
            self.tcts.set(num+1)
            
        
    def getNewHardTorCts(self, event=None):
        torsNum =tkinter.StringVar(master=self.vf.GUI.ROOT)
        torsPAngle=tkinter.StringVar(master=self.vf.GUI.ROOT)
        halfWidth=tkinter.StringVar(master=self.vf.GUI.ROOT)
        ifd2 = self.ifd2 = InputFormDescr(title = "Hard Torsion Constraint Parameters")
        numStr = 'Torsion Number:\n(1-'+str(self.ndihe)+')'
        self.ifd2.append( { 'name': 'tnumEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': numStr,
                'textvariable': torsNum
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':2}})
        self.ifd2.append( { 'name': 'pangEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Perferred Relative Angle(degrees):',
                'textvariable': torsPAngle
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':2}})
        self.ifd2.append( { 'name': 'hwidthEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Full Width of allowed range(degrees):',
                'textvariable': halfWidth
            },
            'gridcfg':{'sticky':tkinter.W, 'columnspan':4}})
        val = self.vf.getUserInput(self.ifd2)
        if val:
            num = int(val['tnumEnt'])
            newEnt=[]
            newEnt.append(val['tnumEnt'])
            newEnt.append(val['pangEnt'])
            newEnt.append(val['hwidthEnt'])
            self.hardTorCts[num]=newEnt


    def getTorCts(self, event=None):
        w=self.ifd.entryByName['barrierLab']
        w1=self.ifd.entryByName['barrierEnt']
        w2=self.ifd.entryByName['barrierCB']
        w3=self.ifd.entryByName['addGCtsBut']
        w4=self.ifd.entryByName['addHCtsBut']
        if self.specifyTorCts.get()=='1':
            w['widget'].grid(w['gridcfg'])
            w1['widget'].grid(w1['gridcfg'])
            w2['widget'].grid(w2['gridcfg'])
            w3['widget'].grid(w3['gridcfg'])
            w4['widget'].grid(w4['gridcfg'])
        else:
            w['widget'].grid_forget()
            w1['widget'].grid_forget()
            w2['widget'].grid_forget()
            w3['widget'].grid_forget()
            w4['widget'].grid_forget()
                

    def getDihe(self, event=None):
        w=self.ifd.entryByName['useInitDiheLab']
        w1=self.ifd.entryByName['dihe0Ent']
        w2=self.ifd.entryByName['initDiheCB']
        if self.specifyDihe.get()=='1':
            w['widget'].grid(w['gridcfg'])    
            w1['widget'].grid(w1['gridcfg'])    
            w2['widget'].grid(w2['gridcfg'])    
        else:
            w['widget'].grid_forget()
            w1['widget'].grid_forget()
            w2['widget'].grid_forget()


    def set_barrier(self):
        if self.barVar.get()=='1':
            self.barrier.set('')


    def set_initDiheType(self):
        if self.initDiheType.get()==1:
            self.dihe0.set('random')
        else:
            nstr= '0. '*self.ndihe
            self.dihe0.set(nstr)


    def set_initTransType(self):
        w=self.ifd.entryByName['tran0Ent']['widget']
        if self.initTransType.get()==1:
            self.tran0.set('random')
        else:
            self.tran0.set(self.about.get())

        
    def set_initQuatType(self):
        if self.initQuatType.get()==1:
            self.quat0.set('random')
        else:
            self.quat0.set('1.0 0. 0. 0.')


    def guiCallback(self):
        #Dpf4InitLigand
        if not hasattr(self.vf.dpo,'ligand'):
            msgStr="Choose/Open Ligand first."
            self.vf.warningMsg(msgStr)
            return
        
        ligandfilename =self.vf.dpo['move']['value']
        ligand = self.vf.dpo.ligand
        ligstring = ligand.types
        if hasattr(ligand, 'ndihe'):
            self.ndihe = ligand.ndihe
        elif hasattr(ligand, 'torscount'):
            self.ndihe = ligand.torscount
        else:
            msgStr="Selected Ligand not formatted!!"
            self.vf.warningMsg(msgStr)
            return
        self.torsdof = ligand.TORSDOF
        d = {}
        for a in ligand.allAtoms:
            d[a.autodock_element] = 1
        autodock_types = list(d.keys())
        for t in self.vf.DPF_FLEXRES_TYPES:
            if t not in autodock_types:
                autodock_types.append(t)
        autodock_types.sort()
            
        autodock_types_str = autodock_types[0]
        for t in autodock_types[1:]:
            autodock_types_str = autodock_types_str + " " + t
        self.types.set(autodock_types_str)
        self.vf.dpo['ligand_types']['value'] = autodock_types_str
        ligand.autodock_types = autodock_types_str
        if ligand==None:
            msgStr="No Ligand Selected!"
            self.vf.warningMsg(msgStr)
            return
        if not hasattr(self, 'form'):
            self.buildForm()
            self.form = self.vf.getUserInput(self.ifd, modal=0, blocking=0)
            self.form.root.protocol('WM_DELETE_WINDOW',self.Close_cb)
            self.initTransType.set(1)
            self.initDiheType.set(1)
            self.initQuatType.set(1)
            self.fmap.set(0)
            #this overwrites previous values
            self.set_initTransType()
            self.set_initQuatType()
            self.set_initDiheType()
        else:
            self.form.root.deiconify()
        self.torsdofcoeff.set("")
        self.vf.dpo['torsdof4']['value'][0] = ligand.TORSDOF
        if not hasattr(self, 'old_ligand') or self.old_ligand!=ligand:
            v = self.vf.dpo['about']['value']
            strCenter = "%6.3f %6.3f %6.3f" %(v[0], v[1], v[2])
            self.about.set(strCenter)
            self.old_ligand=ligand
            self.ligMsgStr.set( 'Ligand:   '+ ligandfilename)
            if self.ndihe:
                self.ndiheMsgStr.set('Number of Active Torsions in Ligand:   '+ str(self.ndihe))
            else:
                self.ndiheMsgStr.set('No Torsions Specified in Ligand!')
            if self.torsdof:
                self.tdofMsgStr.set('Number of Torsional Degrees of Freedom(torsdof) in Ligand:   '+ str(self.torsdof))
            else:
                self.tdofMsgStr.set('No Torsional Degrees of\nFreedom in Ligand!  ')
            self.typesMsgStr.set('Ligand Atom Types: ' + self.types.get())
            self.centerMsgStr.set('Center of Ligand Molecule:  ' + self.about.get())


Dpf4InitLigandGUI=CommandGUI()
Dpf4InitLigandGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'], menuText['AdjustLigand4'], cascadeName = menuText['SetLigandParmsMB'])



class Dpf4LigandChooser(MVCommand):
    """ allows user to choose a molecule already present for the ligand"""


    def onAddCmdToViewer(self):
        checkHasDpo(self.vf)
        #self.vf.loadModule('autotorsCommands', 'AutoDockTools')
        if not hasattr(self.vf, 'atorsDict'):
            self.vf.atorsDict = {}
        self.vf.loadModule('displayCommands')
        self.vf.loadModule('bondsCommands')
        self.vf.loadModule('fileCommands')


    def chooseLigand_cb(self, event = None):
        """called each time the 'choose Ligand' button is pressed"""
        mol = self.chooser.getMolSet()
        if mol: 
            ask = 1
            self.doitWrapper(mol, ask, log=1, redraw=1)
        self.chooser.form.withdraw()


    def guiCallback(self):
        self.mode = 'single'
        self.title= 'Choose Ligand'
        self.chooser = MoleculeChooser(self.vf, self.mode, self.title)
        self.chooser.ipf.append({'name':'Select Button',
                                 'widgetType':tkinter.Button,
                                 'text':'Select Ligand',
                                 'wcfg':{'bd':6},
                                 'gridcfg':{'sticky':tkinter.E+tkinter.W},
                                 'command': self.chooseLigand_cb})
        self.form = self.chooser.go(modal=0, blocking=0)
        lb = self.chooser.ipf.entryByName['Molecule']['widget'].lb
        lb.bind("<Double-Button-1>",self.chooseLigand_cb)


    def __call__(self, nodes, ask=1, **kw):
        self.doitWrapper(*(nodes, ask), **kw)


    def doit(self, nodes, ask):
        lig = self.vf.expandNodes(nodes)[0]
        d = {}
        for a in lig.allAtoms:
            d[a.autodock_element] = 1
        self.vf.DPF_LIGAND_TYPES = list(d.keys())
        isAtorsMol = 0
        isDpoLigand = 0

        hasTorTree = hasattr(lig, 'torTree')

        d = self.vf.atorsDict
        if 'molecule' in d and d['molecule']==lig:
            isAtorsMol = 1
            lig.ndihe = lig.torscount
        if hasattr(self.vf, 'dpo'):
            obj = self.vf.dpo
            if  hasattr(obj, 'ligand') and obj.ligand == lig:
                isDpoLigand = 1

        if hasTorTree or isAtorsMol or isDpoLigand:
            setDpoFields(lig, self.vf.dpo)
        else:
            msgStr = 'can only selectLigand\n preformatted by autotors\
            (and Written to a file)\n for this option!'
            self.vf.warningMsg(msgStr)
            return 'ERROR'

        if self.vf.hasGui: 
            self.vf.colorByAtomType(lig, topCommand = 0, redraw=1)
            aromaticCs = AtomSet([x for x in lig.allAtoms if x.autodock_element=='A'])
            if len(aromaticCs):
                self.vf.color(aromaticCs,((0.,1.,0.),),['lines'],topCommand=0, redraw=1)
            if ask:
                self.vf.ADdpf4_initLigand.guiCallback()
            else:
                d = {}
                if hasattr(lig, 'outputfile'):
                    d['move'] = lig.outputfile
                    d['torsdof4'] = [lig.torsdof, self.vf.dpo['torsdof4']['value'][1]]
                    d['about'] = lig.center
                    d['ndihe'] = lig.ndihe
                    d['types'] = lig.types
                    w = {}
                    for a in lig.allAtoms:
                        w[a.autodock_element] = 1
                    ligtypes = list(w.keys())
                    ligtypes.sort()
                    ligtypestr = ligtypes[0]
                    for t in ligtypes[1:]:
                        lig_type_str = lig_type_str + " " + t
                    d['ligand_types'] = lig_type_str
                    self.vf.ADdpf_setDpo(*(), **d)


Dpf4LigandChooserGUI=CommandGUI()
Dpf4LigandChooserGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'],\
menuText['ChooseLigand4'], cascadeName = menuText['SetLigandParmsMB'])



class Dpf4LigPDBQReader(MVCommand):
    """ allows user to choose a PDBQ file for the ligand"""

    def onAddCmdToViewer(self):
        if self.vf.hasGui and not hasattr(self.vf,'readMolecule'):
            self.vf.loadCommand('fileCommands', 'readMolecule', 'Pmv')
        checkHasDpo(self.vf)


    def guiCallback(self):
        """called each time the 'select ligand' button is pressed"""
        ligFile = self.vf.askFileOpen(types=[ ('PDBQT files:', '*.pdbqt'), 
                        ('Autotors-format convention:', '*.out.pdbqt')],
                        title = 'Formatted Ligand File:')
        if ligFile:
            self.doitWrapper(ligFile, log=1, redraw=1)


    def __call__(self, ligFile, ask = 1, **kw):
        """None<-ADdpf4_readFormattedLigand
        ligFile: file containing the ligand 
        ask: flag for whether to update geometry 
"""
        if not os.path.exists(ligFile):
            raise IOError
        kw['ask'] = ask
        self.doitWrapper(*(ligFile,), **kw)


    def doit(self, ligFile, ask=1):
        ligand = self.vf.readMolecule(ligFile)[0]
        ligFile = os.path.basename(ligFile)
        assert hasattr(ligand, 'TORSDOF'), 'ligand must be preformatted with autotors'
        setDpoFields(ligand, self.vf.dpo)
        #color and check charges here
        d = {}
        for a in ligand.allAtoms:
            d[a.autodock_element] = 1
        self.vf.DPF_LIGAND_TYPES = list(d.keys())
        if self.vf.hasGui and hasattr(ligand,'ndihe') and ask:
            self.vf.colorByAtomType(ligand, topCommand = 0, redraw=1)
            aromaticCs = AtomSet([x for x in ligand.allAtoms if x.autodock_element=='A'])
            if len(aromaticCs):
                self.vf.color(aromaticCs,((0.,1.,0.),),['lines'],topCommand=0, redraw=1)
            self.vf.ADdpf4_initLigand.guiCallback()

        
Dpf4LigPDBQReaderGUI = CommandGUI()
Dpf4LigPDBQReaderGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'], menuText['ReadLigand4'],  cascadeName = menuText['SetLigandParmsMB'])



class DpfEditor(MVCommand):
    """ allows user to edit current output file and write it"""


    def onAddCmdToViewer(self):
        checkHasDpo(self.vf)


    def __call__(self,  **kw):
        self.doitWrapper(*(), **kw)


    def doit(self):
        if self.vf.dpo.dpf_written_filename:
            filename = self.vf.dpo.dpf_written_filename
            fptr=open(filename,'r')
            allLines=fptr.readlines()
            ss=''
            for item in allLines:
                ss=ss+item
        else:
            ss = ''
            filename=''
        titleStr = 'Edit  dpf'
        self.ifd = ifd  = InputFormDescr(title = titleStr)
        ifd.append({'name': 'dpfText',
            'size':[80,30],
            'label':filename,
            'defaultValue':ss,
            'widgetType':'ScrolledText',
            'writeFileType':[('Docking Parameter Files','*.dpf'),('Vina Config Files', '*.txt')],
            'readFileType':[('Docking Parameter Files','*.dpf'),('Vina Config Files', '*.txt')],
            'readButton':1,'writeButton':1})
        vals = self.vf.getUserInput(ifd)


    def guiCallback(self):
        self.doitWrapper(log=1,redraw=0)


DpfEditorGUI=CommandGUI()
DpfEditorGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'], menuText['EditDpfMB'])



class DpfSAWriter(MVCommand):
    """ allows user to choose an output filename and write simulated annealing parameters"""

    def onAddCmdToViewer(self):
        checkHasDpo(self.vf)


    def __call__(self, outfile, **kw):
        self.doitWrapper(*(outfile,), **kw)


    def doit(self, outfile):
        #to remove pickle problem, assume dpo is current self.vf.dpo
        if not len(self.vf.dpo.receptor_stem):
            self.vf.warningMsg("You must choose a macromolecule before writing dpf")
            return 'ERROR'

        #if a rms reference file has been specified, write it to dpf
        if self.vf.dpo['rmsref']['value']!=self.vf.dpo['move']['value']:
            l = simulated_annealing_list
            ind = l.index('rmstol') + 1
            l.insert(ind, 'rmsref')
            
        self.vf.dpo.write(outfile, simulated_annealing_list)
        #this is set when the dpo is written
        #self.vf.dpo.dpf_filename = outfile


    def guiCallback(self):
        if not len(self.vf.dpo.receptor_stem):
            self.vf.warningMsg("You must choose a macromolecule before writing dpf")
            return 'ERROR'
        outfile = self.vf.askFileSave(types=[('dpf file', '*.dpf')],
                title = 'SA Docking Parameter Output File:')
        if outfile: 
            self.doitWrapper(outfile, log=1,redraw=0)


DpfSAWriterGUI=CommandGUI()
DpfSAWriterGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'],\
menuText['WriteSA'], cascadeName = menuText['WriteDpfMB'])



class Dpf41EPDBWriter(MVCommand):
    """ allows user to choose an output filename and write it"""


    def onAddCmdToViewer(self):
        checkHasDpo(self.vf)


    def __call__(self, outfile, **kw):
        self.doitWrapper(*(outfile,), **kw)


    def doit(self, outfile):
        if not len(self.vf.dpo.receptor_stem):
            self.vf.warningMsg("You must choose a macromolecule before writing dpf")
            return 'ERROR'
        orig_epdb_flag = self.vf.dpo['epdb_flag']['value']
        self.vf.dpo['epdb_flag']['value'] = 1
        self.vf.dpo.write42(outfile, epdb_list4_2)
        #msg = "restoring dpo epdb_flag value from current "+ str(self.vf.dpo['epdb_flag']['value'])
        self.vf.dpo['epdb_flag']['value'] = orig_epdb_flag
        #msg= msg + " to " + str(self.vf.dpo['epdb_flag']['value'])
        #self.vf.warningMsg(msg)


        
    def guiCallback(self):
        if not len(self.vf.dpo.receptor_stem):
            self.vf.warningMsg("You must choose a macromolecule before writing dpf")
            return 'ERROR'
        outfile = self.vf.askFileSave(types=[('dpf file', '*.dpf')],
                title = 'LS Docking Parameter Output File:')
        if outfile: 
            self.doitWrapper(outfile, log=1,redraw=0)


Dpf41EPDBWriterGUI=CommandGUI()
Dpf41EPDBWriterGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'], menuText['WriteEPDB41'],\
cascadeName = menuText['WriteDpfMB'])



class Dpf41CONFIGWriter(MVCommand):
    """ allows user to write a new config file for AutoDockVina"""


    def onAddCmdToViewer(self):
        if not hasattr(self.vf, 'vinaDict'): self.vf.vinaDict={}
        if not hasattr(self.vf, 'readMolecule'):
            self.vf.loadCommand('fileCommands', 'readMolecule', 'Pmv')
        self.vina_keys = ['receptor', 'ligand', 'flexres',\
                          'center_x', 'center_y', 'center_z', \
                          'size_x', 'size_y', 'size_z', 'out']


    def browseRec(self, event=None):
        rec = self.vf.askFileOpen(types=[('PDBQT files', '*.pdbqt')],
                title = 'Receptor File(PDBQT):')
        if rec: 
            pfname=os.path.split(rec)[-1]
            self.receptor.set(pfname) #StringVar
            self.ifd.entryByName['receptorEnt']['wcfg']['textvariable'].set(pfname)
            #msg = "set receptor to " + pfname
            #self.vf.warningMsg(msg)
            self.vf.vinaDict['receptor'] = pfname


    def browseLig(self, event=None):
        lig = self.vf.askFileOpen(types=[('PDBQT files', '*.pdbqt')],
                title = 'Ligand File(PDBQT):')
        if lig: 
            #don't strip off the pathname??
            pfname=os.path.split(lig)[-1]
            self.ligand.set(pfname)
            self.ifd.entryByName['ligandEnt']['wcfg']['textvariable'].set(pfname)
            self.vf.vinaDict['ligand'] = pfname


    def browseOut(self, event=None):
        Out = self.vf.askFileSave(types=[('PDBQT files', '*.pdbqt')],
                title = 'Vina Result File(PDBQT):')
        if Out: 
            #don't strip off the pathname??
            pfname=os.path.split(Out)[-1]
            self.out.set(pfname)
            self.ifd.entryByName['out']['wcfg']['textvariable'].set(pfname)
            self.vf.vinaDict['out'] = pfname


    def browseFlexRes(self, event=None):
        flexres = self.vf.askFileOpen(types=[('PDBQT files', '*.pdbqt')],
                title = 'Flexible Residues File(PDBQT):')
        if flexres: 
            #don't strip off the pathname??
            pfname=os.path.split(flexres)[-1]
            self.flexres.set(pfname)
            self.ifd.entryByName['flexEnt']['wcfg']['textvariable'].set(pfname)
            self.vf.vinaDict['flexres'] = pfname


    def browseRefConfig(self, event=None):
        refconfig = self.vf.askFileOpen(types=[('txt files', '*.txt')],
                title = 'Vina Configuration File(txt):')
        if refconfig: 
            #don't strip off the pathname??
            pfname=os.path.split(refconfig)[-1]
            self.ifd.entryByName['refconfigEnt']['wcfg']['textvariable'].set(pfname)
            #config.set(pfname)
            self.vf.vinaDict['refconfigEnt'] = pfname
            fptr = open(pfname)
            lines = fptr.readlines()
            fptr.close()
            for ll in lines:
                llist = ll.split()
                if len(llist):
                    print(llist)
                    key = llist[0]
                    value = llist[-1]
                    if key in ['ligand', 'receptor','energy_range','center_x', 'center_y', 'center_z', 
                               'log','outputconfig','flexres','exhaustiveness',
                               'size_x','size_y','size_z','num_modes','out','refconfig','numCPUs','seed']:
                        key = key + "Ent"
                    self.ifd.entryByName[key]['wcfg']['textvariable'].set(value)


    def browseOutputConfig(self, event=None):
        output = self.vf.askFileSave(types=[('txt files', '*.txt')],
                title = 'Vina Config Filename(.txt):')
        if output: 
            #don't strip off the pathname??

            pfname=os.path.split(output)[-1]
            self.ifd.entryByName['outputconfig']['wcfg']['textvariable'].set(pfname)
            #config.set(pfname)
            self.vf.vinaDict['outputconfig'] = pfname


    def updateReceptor(self, event=None):
        self.vf.warningMsg("in updateReceptor")
        self.vf.warningMsg( self.receptor.get()) #StringVar
        self.ifd.entryByName['receptor']['wcfg']['textvariable'].set(self.receptor.get()) #StringVar


    def updateLigand(self, event=None):
        self.vf.warningMsg("in updateLigand")
        self.vf.warningMsg( self.ligand.get())
        self.ifd.entryByName['ligand']['wcfg']['textvariable'].set(self.ligand.get())


    def showVinaOpts(self, event=None):
        self.show_VinaOpts.set(1)
        self.ifd.entryByName['scoreOnlyCB']['widget'].grid(self.ifd.entryByName['scoreOnlyCB']['gridcfg'])
        self.ifd.entryByName['localOnlyCB']['widget'].grid(self.ifd.entryByName['localOnlyCB']['gridcfg'])
        self.ifd.entryByName['randomizeOnlyCB']['widget'].grid(self.ifd.entryByName['randomizeOnlyCB']['gridcfg'])
        #more options: numCPUs,seed,exhaustiveness, num_modes, energy_range
        #row 8
        self.ifd.entryByName['numCPUsLab']['widget'].grid(self.ifd.entryByName['numCPUsLab']['gridcfg'])
        self.ifd.entryByName['numCPUsEnt']['widget'].grid(self.ifd.entryByName['numCPUsEnt']['gridcfg'])
        self.ifd.entryByName['num_modesLab']['widget'].grid(self.ifd.entryByName['num_modesLab']['gridcfg'])
        self.ifd.entryByName['num_modesEnt']['widget'].grid(self.ifd.entryByName['num_modesEnt']['gridcfg'])
        self.ifd.entryByName['exhaustivenessLab']['widget'].grid(self.ifd.entryByName['exhaustivenessLab']['gridcfg'])
        self.ifd.entryByName['exhaustivenessEnt']['widget'].grid(self.ifd.entryByName['exhaustivenessEnt']['gridcfg'])
        #row 9
        self.ifd.entryByName['seedLab']['widget'].grid(self.ifd.entryByName['seedLab']['gridcfg'])
        self.ifd.entryByName['seedEnt']['widget'].grid(self.ifd.entryByName['seedEnt']['gridcfg'])
        self.ifd.entryByName['energy_rangeLab']['widget'].grid(self.ifd.entryByName['energy_rangeLab']['gridcfg'])
        self.ifd.entryByName['energy_rangeEnt']['widget'].grid(self.ifd.entryByName['energy_rangeEnt']['gridcfg'])
        self.form.autoSize()
 ###        #??weights
 ###        print "--weight_gauss1" #-0.035579
 ###        print "--weight_gauss2" #-0.005156
 ###        print "--weight_repulsion" #0.840245
 ###        print "--weight_hydrophobic"#-0.035069
 ###        print "--weight_hydrogen"#-0.587439
 ###        print "--weight_rot" #-0.05845999999
 ###        #misc
 ###        print "--cpu" #try to detect number of CPUs; defaults to 1
 ###        print "--seed" #explicit random seed
 ###        print "--exhaustiveness" #roughly proportional to time: 1+
 ###        print "--num_modes" #9: maximum number of binding modes to generate
 ###        print "--energy_range" #3: maximum energy difference between the best
 ###                               # binding mode and the worst one displayed
 ###                               # (kcal/mol
 
    def hideVinaOpts(self, event=None):
        self.show_VinaOpts.set(0)
        self.ifd.entryByName['scoreOnlyCB']['widget'].grid_forget()
        self.ifd.entryByName['localOnlyCB']['widget'].grid_forget()
        self.ifd.entryByName['randomizeOnlyCB']['widget'].grid_forget()
        #more options: numCPUs,seed,exhaustiveness, num_modes, energy_range
        #row 8
        self.ifd.entryByName['numCPUsLab']['widget'].grid_forget()
        self.ifd.entryByName['numCPUsEnt']['widget'].grid_forget()
        self.ifd.entryByName['num_modesLab']['widget'].grid_forget()
        self.ifd.entryByName['num_modesEnt']['widget'].grid_forget()
        self.ifd.entryByName['exhaustivenessLab']['widget'].grid_forget()
        self.ifd.entryByName['exhaustivenessEnt']['widget'].grid_forget()
        #row 8
        self.ifd.entryByName['seedLab']['widget'].grid_forget()
        self.ifd.entryByName['seedEnt']['widget'].grid_forget()
        self.ifd.entryByName['energy_rangeLab']['widget'].grid_forget()
        self.ifd.entryByName['energy_rangeEnt']['widget'].grid_forget()
        self.form.autoSize()
     
 
    def clearVinaOpts(self, event=None):
        self.scoreOnly.set(0)
        self.localOnly.set(0)
        self.randomizeOnly.set(0)
     
 
    def readFldFile(self, event=None):
        fldFile = self.vf.askFileOpen(types=[('select autogrid fld filename:', '*.fld')],
                 title = 'AutoGrid fld File:')
        if fldFile:
             fptr = open(fldFile)
             lines = fptr.readlines()
             fptr.close()
             SPACING=-1
             NELEMENTS=[]
             CENTER=[]
             found_spacing=False
             found_nelements=False
             found_ncenter=False
             #lines from fld files:
             #SPACING 0.375
             #NELEMENTS 60 60 60
             #CENTER 2.500 6.500 -7.500
             for ll in lines:
                 if ll.find("SPACING")>-1:
                     SPACING = float(ll.strip().split()[1])
                     print("SPACING set to ", SPACING)
                     found_spacing=True
                 if ll.find("NELEMENTS")>-1:
                     NELEMENTS = list(map(float,ll.strip().split()[1:4]))
                     print("NELEMENTS set to ", NELEMENTS)
                     found_nelements=True
                 if ll.find("CENTER")>-1:
                     CENTER = list(map(float,ll.strip().split()[1:4]))
                     print("CENTER set to ", CENTER)
                     found_center=True
             if found_center and found_spacing and found_nelements:
                 self.center_x.set(str(CENTER[0]))
                 self.center_y.set(str(CENTER[1]))
                 self.center_z.set(str(CENTER[2]))
                 size_x = SPACING*(NELEMENTS[0]-1)
                 size_y = SPACING*(NELEMENTS[1]-1)
                 size_z = SPACING*(NELEMENTS[2]-1)
                 self.size_x.set(str(size_x))
                 self.size_y.set(str(size_y))
                 self.size_z.set(str(size_z))
 

    def set_scoreOnly(self, event=None):
         #score_only: can omit search space
        self.scoreOnly.set(1)
        self.localOnly.set(0)
        self.randomizeOnly.set(0)
 

    def set_localOnly(self, event=None):
        #local_only: only local search
        self.localOnly.set(1)
        self.scoreOnly.set(0)
        self.randomizeOnly.set(0)
 

    def set_randomizeOnly(self, event=None):
        #randomize_only, attempting to avoid clashes 
        self.randomizeOnly.set(1)
        self.localOnly.set(0)
        self.scoreOnly.set(0)


    def buildForm(self):
        #vina input variables:
        self.receptor=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.ligand=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.flexres=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.center_x=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.center_y=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.center_z=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.size_x=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.size_y=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.size_z=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.outputconfig=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.refconfig=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.out=tkinter.StringVar(master=self.vf.GUI.ROOT) #RESULT FILE!
        self.log=tkinter.StringVar(master=self.vf.GUI.ROOT) #1
        self.log.set("")
        self.cpu_number=tkinter.StringVar(master=self.vf.GUI.ROOT) #2
        self.seed=tkinter.StringVar(master=self.vf.GUI.ROOT) #3
        self.exhaustiveness=tkinter.StringVar(master=self.vf.GUI.ROOT) #4
        self.exhaustiveness.set("8")
        self.num_modes=tkinter.StringVar(master=self.vf.GUI.ROOT) #5
        self.num_modes.set("9")
        self.energy_range=tkinter.StringVar(master=self.vf.GUI.ROOT) #6
        self.energy_range.set("3")
        self.scoreOnly= tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.localOnly= tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.randomizeOnly= tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.show_VinaOpts = tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.scoreOnly.set(0)
        self.localOnly.set(0)
        self.randomizeOnly.set(0)
        self.show_VinaOpts.set(0)
        self.outputconfig.set("config.txt")
        ifd = self.ifd = InputFormDescr(title = "Vina Input Parameters")
        ifd.append( {'name': 'receptorLab', #R0_C0-2
            'widgetType':tkinter.Label,
            'tooltip': 'input receptor pdbqt filename',
            'wcfg':{'text': 'Receptor:' },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':3 }})
        ifd.append({'name': 'receptorButton', #R0_C3
            'widgetType': tkinter.Button,
            'tooltip': 'open file browser to select receptor file (pdbqt)',
            'text':'Browse',
            'wcfg':{'bd':1},
            'command':self.browseRec,
            'gridcfg':{'sticky':tkinter.E,'row':-1, 'column':3}})
        ifd.append( {'name': 'ligandLab', #R0_C4-6
            'widgetType':tkinter.Label,
            'tooltip': 'input ligand pdbqt filename',
            'wcfg':{'text': '   Ligand: ' },
            'gridcfg':{'sticky':tkinter.E, 'row':-1, 'column':4, 'columnspan':3 }})
        ifd.append({'name': 'ligandButton', #R0_C7
            'widgetType': tkinter.Button,
            'tooltip': 'open file browser to select ligand file (pdbqt)',
            'text':'Browse',
            'wcfg':{'bd':1},
            'command':self.browseLig,
            'gridcfg':{'sticky':tkinter.W,'row':-1, 'column':7}})
        ifd.append( {'name': 'flexLab', #R0_C8-10
            'widgetType':tkinter.Label,
            'tooltip': '[optional] input flexible residue pdbqt filename',
            'wcfg':{'text': '   Flex: '},
            'gridcfg':{'sticky':tkinter.E,'row':-1, 'column':8, 'columnspan':3}})
        ifd.append({'name': 'flexButton', #R0_C11
            'widgetType': tkinter.Button,
            'text':'Browse',
            'tooltip': 'open file browser to select [optional] flexible residue file (pdbqt)',
            'wcfg':{'bd':1},
            'command':self.browseFlexRes,
            'gridcfg':{'sticky':tkinter.W,'row':-1, 'column':11}})
        # receptor_entry  ligand_entry  flex_entry
        ifd.append( {'name': 'receptorEnt',
            'widgetType':tkinter.Entry, #R1_C0-3
            'tooltip': 'receptor pdbqt filename',
            'wcfg':{ 'textvariable': self.receptor,'width':4 },
            'gridcfg':{'sticky':tkinter.W+tkinter.E, 'columnspan':4}})
        ifd.append( {'name': 'ligandEnt', #R1_C4-7
            'widgetType':tkinter.Entry,
            'tooltip': 'ligand pdbqt filename',
            'wcfg':{ 'textvariable': self.ligand,'width':4 },
            'gridcfg':{'sticky':tkinter.W+tkinter.E, 'row':-1, 'columnspan':4, 'column':4}})
        ifd.append( {'name': 'flexEnt', #R1_C8-11
            'widgetType':tkinter.Entry,
            'tooltip': '[optional] flexible residue pdbqt filename',
            'wcfg':{ 'textvariable': self.flexres,'width':4},
            'gridcfg':{'sticky':tkinter.W+tkinter.E, 'row':-1, 'columnspan':4, 'column':8}})
        #center
        ifd.append( {'name': 'centerLab', #R2_C0-2
            'widgetType':tkinter.Label,
            'tooltip': 'xyz center of search space coords',
            'wcfg':{'text': 'Center:'},
            'gridcfg':{'sticky':tkinter.E, 'columnspan':3}})
        ifd.append( {'name': 'center_xEnt', #R2_C3
            'widgetType':tkinter.Entry,
            'tooltip': 'x coordinate of search space center',
            'wcfg':{ 'textvariable': self.center_x, 'width':3 },
            'gridcfg':{'sticky':tkinter.E+tkinter.W, 'row':-1, 'column':3, 'columnspan':1}})
        ifd.append( {'name': 'center_yEnt', #R2_C4
            'widgetType':tkinter.Entry,
            'wcfg':{ 'textvariable': self.center_y, 'width':3 },
            'tooltip': 'y coordinate of search space center',
            'gridcfg':{'sticky':tkinter.E+tkinter.W, 'row':-1, 'column':4, 'columnspan':1}})
        ifd.append( {'name': 'center_zEnt', #R2_C5
            'widgetType':tkinter.Entry,
            'tooltip': 'z coordinate of search space center',
            'wcfg':{ 'textvariable': self.center_z, 'width':3 },
            'gridcfg':{'sticky':tkinter.E+tkinter.W, 'row':-1, 'column':5, 'columnspan':1}})
        #size
        ifd.append( {'name': 'sizeLab', #R2_C6-8
            'widgetType':tkinter.Label,
            'tooltip': 'Search space dimensions(Angstrom)',
            'wcfg':{ 'text': 'Size:' },
            'gridcfg':{'sticky':tkinter.E+tkinter.W, 'row':-1, 'column':6,'columnspan':3}})
        ifd.append( {'name': 'size_xEnt',#R2_C9
            'widgetType':tkinter.Entry,
            'tooltip': 'x dimension of search space in Angstrom',
            'wcfg':{ 'textvariable': self.size_x, 'width':5 },
            'gridcfg':{'sticky':tkinter.E+tkinter.W, 'row':-1, 'column':9, 'columnspan':1}})
        ifd.append( {'name': 'size_yEnt',#R2_C10
            'widgetType':tkinter.Entry,
            'tooltip': 'y dimension of search space in Angstrom',
            'wcfg':{ 'textvariable': self.size_y, 'width':5 },
            'gridcfg':{'sticky':tkinter.E+tkinter.W, 'row':-1, 'column':10, 'columnspan':1}})
        ifd.append( {'name': 'size_zEnt',#R2_C11
            'widgetType':tkinter.Entry,
            'tooltip': 'z dimension of search space in Angstrom',
            'wcfg':{ 'textvariable': self.size_z, 'width':5 },
            'gridcfg':{'sticky':tkinter.E+tkinter.W, 'row':-1, 'column':11, 'columnspan':1}})
        
        #------------------------------------------------------------------
        #output result filename
        #------------------------------------------------------------------
        ifd.append( {'name': 'outEnt',  #R3_C0-3 @@ ZERO-BASED!
            'widgetType':tkinter.Entry,
            'tooltip': '[optional] output filename for vina results.\nDefault is ligandname_out.pdbqt',
            'wcfg':{'label': 'Out:',
                'textvariable': self.out,
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':4}})
        ifd.append({'name': 'outButton', #R3_C4-5
            'widgetType': tkinter.Button,
            'text':'Browse',
            'wcfg':{'bd':2},
            'command':self.browseOut,
            'gridcfg':{'sticky':tkinter.W+tkinter.E,'row':-1, 'column':4, 'columnspan':2}})
        ifd.append( {'name': 'logEnt',  #R3_C6-11 @@ ZERO-BASED!
            'widgetType':tkinter.Entry,
            'tooltip': '[optional] write log file to this filename',
            'wcfg':{'label': ' Log:',
                'textvariable': self.log,
            },
            'gridcfg':{'sticky':tkinter.W+tkinter.E, 'row':-1, 'column':6, 'columnspan':6}})
        #output config filename
        ifd.append( {'name': 'refconfigEnt',#R4_C0-5
            'widgetType':tkinter.Entry,
            'tooltip': '[optional] load values from reference filename containing vina input parameters',
            'wcfg':{'label': '   Refconfig filename:', 
                    'textvariable': self.refconfig,
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':6}})
        ifd.append({'name': 'refconfigButton', #R4_C6-11
            'widgetType': tkinter.Button,
            'tooltip': '[optional] use parameters from reference config file',
            'text':'Browse',
            'wcfg':{'bd':2},
            'command':self.browseRefConfig,
            'gridcfg':{'sticky':tkinter.W+tkinter.E,'row':-1, 'column':6, 'columnspan':6}})
        #output config filename
        ifd.append( {'name': 'outputconfigEnt', #R5_C0-5
            'widgetType':tkinter.Entry,
            'tooltip': '[optional] specify filename for vina configuration file instead of default "config.txt"',
            'wcfg':{'label': 'Output filename:',
                'textvariable': self.outputconfig,
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':6}})
        #output config filename
        ifd.append({'name': 'outputButton', #R5_C6-11
            'widgetType': tkinter.Button,
            'tooltip': "[optional] select filename for new config file.\nDefault is 'config.txt'",
            'text':'Browse',
            'wcfg':{'bd':2},
            'command':self.browseOutputConfig,
            'gridcfg':{'sticky':tkinter.W+tkinter.E,'row':-1, 'column':6, 'columnspan':6}})
##############################################################
#  10.1.2012: @@ NOW-> 12 columns in each row
##############################################################
        #show/hide other vina options buttons
        ifd.append({'name': 'displayVinaOptsButton', #R6_C0-3 vina:score_only,local_only, randomize_only
             'widgetType':tkinter.Button,
             'text':'Show options',
             'tooltip': 'set vina options: score_only, local_only, randomize_only',
             'command':self.showVinaOpts,
             'wcfg':{'bd':2},
             'gridcfg':{'sticky':tkinter.W+tkinter.E, 'columnspan':4}})
        ifd.append({'name': 'hideVinaOptsButton', #R6_C4-7
             'widgetType': tkinter.Button,
             'text':'Hide options',
             'tooltip': 'hide vina options: score_only, local_only, randomize_only',
             'wcfg':{'bd':2},
             'command':self.hideVinaOpts,
             'gridcfg':{'sticky':tkinter.W+tkinter.E,'row':-1, 'columnspan':4, 'column':4}})
        ifd.append({'name': 'clearVinaOptsButton', #R6_C8-11
             'widgetType': tkinter.Button,
             'text':'Clear',
             'tooltip': 'turn off vina options: score_only, local_only, randomize_only',
             'wcfg':{'bd':2},
             'command':self.clearVinaOpts,
             'gridcfg':{'sticky':tkinter.W+tkinter.E,'row':-1, 'columnspan':4, 'column':8}})
        ifd.append({'name': 'scoreOnlyCB', #R7_C0-3
             'widgetType':tkinter.Checkbutton,
             'text': 'score_only ',
             'variable': self.scoreOnly,
             'command': self.set_scoreOnly,
             'gridcfg':{'sticky':tkinter.W, 'column':0, 'columnspan':4 }})
        ifd.append({'name': 'localOnlyCB',#R7_C4-7
             'widgetType':tkinter.Checkbutton,
             'text': 'local_only ',
             'variable': self.localOnly,
             'command': self.set_localOnly,
             'gridcfg':{'sticky':tkinter.W, 'row':-1, 'column':4, 'columnspan':4 }})
        ifd.append({'name': 'randomizeOnlyCB',#R7_C8-11
             'widgetType':tkinter.Checkbutton,
             'text': 'randomize_only ',
             'variable': self.randomizeOnly,
             'command': self.set_randomizeOnly,
             'gridcfg':{'sticky':tkinter.W, 'row':-1, 'column':8, 'columnspan':4 }})
        #more options: cpu,seed,exhaustiveness,num_modes, energy_range
        #cpus, exhaustiveness, modes
        #------------------------------------------------------------------
        # shorter strings: cpus, modes, exhaustiveness (4 each)
        #------------------------------------------------------------------
        ifd.append( {'name': 'numCPUsLab',#R8_C0-1
            'widgetType':tkinter.Label,
            'wcfg':{'text': 'cpus' },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':2}})
        ifd.append( {'name': 'numCPUsEnt',#R8_C2-3
            'widgetType':tkinter.Entry,
            'tooltip': 'number of CPUs to use (default is to try to detect number of CPUs or, failing that, use 1)',
            'wcfg':{ 'textvariable': self.cpu_number },
            'gridcfg':{'sticky':tkinter.E, 'row':-1, 'column':2, 'columnspan':2}})
        ifd.append( {'name': 'num_modesLab',#R8_C4-5
            'widgetType':tkinter.Label,
            'wcfg':{'text': 'modes' },
            'gridcfg':{'sticky':tkinter.E, 'row':-1, 'column':4,'columnspan':2}})
        ifd.append( {'name': 'num_modesEnt',#R8_C6-7
            'widgetType':tkinter.Entry,
            'tooltip': 'maximum number of binding modes to generate\n(default is 9)',
            'wcfg':{ 'textvariable': self.num_modes},
            'gridcfg':{'sticky':tkinter.E, 'row':-1, 'column':6, 'columnspan':2}})
        ifd.append( {'name': 'exhaustivenessLab',#R8_C8-9
            'widgetType':tkinter.Label,
            'wcfg':{'text': 'exh.' },
            'gridcfg':{'sticky':tkinter.E, 'row':-1, 'column':8,'columnspan':2}})
        ifd.append( {'name': 'exhaustivenessEnt',#R8_C10-11
            'widgetType':tkinter.Entry,
            'tooltip': 'exhaustiveness of the global search\n(roughly proportional to time): 1+\n(default is 8)',
            'wcfg':{ 'textvariable': self.exhaustiveness, },
            'gridcfg':{'sticky':tkinter.E, 'row':-1, 'column':10, 'columnspan':2}})
        #------------------------------------------------------------------
        # longer strings: seed and energy range (?) (6 each)
        #------------------------------------------------------------------
        ifd.append( {'name': 'seedLab',#R9_C0-2 
            'widgetType':tkinter.Label,
            'wcfg':{'text': 'seed' },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':3}})
        ifd.append( {'name': 'seedEnt',#R9_C3-5
            'widgetType':tkinter.Entry,
            'tooltip': 'explicit random seed',
            'wcfg':{ 'textvariable': self.seed },
            'gridcfg':{'sticky':tkinter.E, 'row':-1, 'column':3, 'columnspan':3}})
        ifd.append( {'name': 'energy_rangeLab',#R9_C6-8
            'widgetType':tkinter.Label,
            'wcfg':{'text': 'energy_range' },
            'gridcfg':{'sticky':tkinter.E, 'row':-1, 'column':6, 'columnspan':3}})
        ifd.append( {'name': 'energy_rangeEnt',#R9_C8-11
            'widgetType':tkinter.Entry,
            'tooltip': 'maximum energy difference between the best binding\nmode and the worst one displayed (kcal/mol)\n(default is 3)',
            'wcfg':{ 'textvariable': self.energy_range },
            'gridcfg':{'sticky':tkinter.E, 'row':-1, 'column':9, 'columnspan':3}})
##############################################################
#       launch/cancel buttons
##############################################################
        ifd.append({'widgetType': tkinter.Button, #R3_C0-5/R8_C0-5
            'text':'Save',
            'tooltip': 'Write new vina config file',
            'wcfg':{'bd':2},
            'gridcfg':{'sticky':tkinter.W+tkinter.E, 'columnspan':6, 'column':0},
            'command':self.Accept_cb})
        ifd.append({'widgetType': tkinter.Button, #R3_C6-11/R8_C6-11
            'text':'Cancel',
            'tooltip': "close 'Vina Input Parameters' panel",
            'wcfg':{'bd':2},
            'gridcfg':{'sticky':tkinter.W+tkinter.E, 'column': 6, 'row':-1,'columnspan':6},
            'command':self.Close_cb})
        self.form = self.vf.getUserInput(self.ifd, scrolledFrame=1, modal=0, blocking=0)
        self.form.root.protocol('WM_DELETE_WINDOW',self.Close_cb)    
    

    def guiCallback(self):
        """called each time 'Docking->Output-> Vina Config(config.txt)' is selected"""
        if self.vf.hasGui: 
            if not hasattr(self, 'ifd'):
                self.buildForm()
                self.topLevel = self.form.root
                self.hideVinaOpts()
            else:
                self.form.deiconify()
            if hasattr(self.vf.gpo, 'ligand_filename') and self.vf.gpo.ligand_filename!="":
                self.ligand.set(self.vf.gpo.ligand_filename)
                self.vf.vinaDict['ligand']=self.vf.gpo.ligand_filename
                ligand_filename = self.vf.gpo.ligand_filename
                #default output = name + "_vina.pdbqt"
                default_out_stem = os.path.basename(ligand_filename).split('.')[0]
                default_out = default_out_stem + "_out.pdbqt"
                self.vf.vinaDict['out']=default_out
                self.out.set(default_out)
            if hasattr(self.vf.gpo, 'receptor_filename') and self.vf.gpo.receptor_filename!="":
                self.receptor.set(self.vf.gpo.receptor_filename)
                self.vf.vinaDict['receptor']=self.vf.gpo.receptor_filename
            if 'npts' in self.vf.gpo: #and self.vf.gpo['npts']['value']!=[0,0,0]:?@@?
                npts = self.vf.gpo['npts']['value']
                spacing = self.vf.gpo['spacing']['value']
                size_x = npts[0]*spacing # npts[0] -1 ??
                size_y = npts[1]*spacing
                size_z = npts[2]*spacing
                self.vf.vinaDict['size_x']=size_x
                self.vf.vinaDict['size_y']=size_y
                self.vf.vinaDict['size_z']=size_z           
                self.size_x.set(size_x)
                self.size_y.set(size_y)
                self.size_z.set(size_z)
            if hasattr(self.vf.dpo, 'flexres') and self.vf.dpo.flexres!="":
                flex_filename = self.vf.dpo.flexres.parser.filename
                self.vf.vinaDict['flex']=flex_filename
                self.vf.vinaDict['flexres']=flex_filename
                self.flexres.set(flex_filename)
            filenames = []
            for m in self.vf.Mols: 
                filenames.append((m, m.parser.filename))
            center_x = ""
            center_y = ""
            center_z = ""
            if 'center_x' in self.vf.vinaDict and 'center_y' in self.vf.vinaDict and 'center_z' in self.vf.vinaDict:
                center_x = self.vf.vinaDict['center_x']
                center_y = self.vf.vinaDict['center_y']
                center_z = self.vf.vinaDict['center_z']
            elif self.vf.gpo['gridcenter']['value']=='auto' and hasattr(self, 'receptor'): #StringVar
                for m, filename in filenames:
                     if filename == self.receptor.get():
                         center_x, center_y, center_z = m.getCenter()
                         break
            elif self.vf.gpo['gridcenter']['value']!='auto':
                center_x, center_y, center_z = self.vf.gpo['gridcenter']['value']
            self.center_x.set(center_x)
            self.center_y.set(center_y)
            self.center_z.set(center_z)
            self.vf.vinaDict['center_x']=center_x
            self.vf.vinaDict['center_y']=center_y
            self.vf.vinaDict['center_z']=center_z
            #??? WOULD THIS BE REASONABLE??? ?@@?
            ###if hasattr(self.vf.dpo, 'receptor_filename') and len(self.vf.dpo.receptor_filename):
            ###    self.receptor.set(self.vf.dpo.receptor_filename)
            ###    if self.vf.vinaDict['receptor']!=self.vf.dpo.receptor_filename:
            ###        msg = "replacing " + self.vf.vinaDict['receptor']+ " with dpo value " + self.vf.dpo.receptor_filename
            ###        self.warningMsg(msg)
            ###        self.vf.vinaDict['receptor']=self.vf.dpo.receptor_filename
                


    def Accept_cb(self, event=None):
        kw = {}
        if self.scoreOnly.get():
            kw['score_only'] = 1
        elif self.localOnly.get():
            kw['local_only'] = 1
        elif self.randomizeOnly.get():
            kw['randomize_only'] = 1
        new_value = self.ifd.entryByName['refconfigEnt']['wcfg']['textvariable'].get()
        if len(new_value):
            kw['refconfig'] = new_value
        else:
            for item in ['receptorEnt','flexEnt','ligandEnt','center_xEnt','center_yEnt','center_zEnt','size_xEnt','size_yEnt','size_zEnt', 'outEnt']:
                new_value = self.ifd.entryByName[item]['wcfg']['textvariable'].get()
                if len(new_value): kw[item[:-3]] = new_value
        outputconfigFile = self.outputconfig.get()
        self.doitWrapper(*(outputconfigFile,), **kw) 
        self.Close_cb()


    def Close_cb(self):
        self.form.root.withdraw()


    def doit(self, outputconfigFile="", **kw):
        """None<-self.ADdpf_writeCONFIG
        outputconfigFile is name of file to write which will contain vina input parameters
        self.workingDirPath = Tkinter.StringVar(master=self.vf.GUI.ROOT, value = os.getcwd() ) # N3P
"""
        if not hasattr(self.vf, 'vinaDict'):
            self.vf.vinaDict = {}
        if 1:
            receptor_filename = kw.pop('receptor')  
            ligand_filename = kw.pop('ligand')
            flex_filename = ""
            if 'flexres' in self.vf.vinaDict: 
                flex_filename = self.vf.vinaDict['flexres']
                if len(flex_filename)==0: 
                    flex_filename = ""
            elif 'flex' in self.vf.vinaDict:
                flex_filename = self.vf.vinaDict['flex']
                if len(flex_filename)==0: 
                    flex_filename = ""
            elif 'flexres' in self.vf.vinaDict:
                flex_filename = self.vf.vinaDict['flexres']
                if len(flex_filename)==0: 
                    flex_filename = ""
                #assert len(flex_filename)>0
                #self.vf.vinaDict['flexres']=flex_filename
                #flexres_filename = kw.pop('flex')  
            #except:
            #    flex_filename = self.vf.vinaDict['flex']
                #flex_filename = flexres_filename
            cm = ConfigFileMaker(receptor=receptor_filename, ligand=ligand_filename, flexres=flex_filename,**kw)
            if not len(outputconfigFile): outputconfigFile = 'config.txt'
            cm.write(outputconfigFile)
            msg =  "wrote "+ outputconfigFile+ " line 2623"
            #self.warningMsg(msg)
            #print "used kws:"
            #for k,v in kw.items(): print k,'~', v
        else:
            if outputconfigFile =="": 
                outputconfigFile = kw.get('outputconfigFile','')
                if outputconfigFile=="":
                    outputconfigFile='config.txt'
            ptr = open(outputconfigFile,'w')
            kw_keys = list(kw.keys())
            options = []
            if len(self.flexres.get()):
                self.vina_keys.insert(2,'flexres')
                #[receptor, ligand, flexres, center_x, center_y, center_z, size_x,size_y,size_z,out]
            for k in self.vina_keys:
                #[receptor, ligand, center_x, center_y, center_z, size_x,size_y,size_z,out]
                self.vf.vinaDict[k] = kw[k]
                line = "%s = %s\n" %(k, kw.pop(k))
                #line = "%s = %s\n" %(k, str(kw[k]))
                ptr.write(line)
            for k in list(kw.keys()):
                self.vf.vinaDict[k] = kw[k]
                if k not in self.vina_keys: 
                    options.append(k)
                    #update vinaDict
            if self.scoreOnly.get() and 'score_only' not in options:
                options.append('score_only')
            if self.localOnly.get() and 'local_only' not in options:
                options.append('local_only')
            if self.randomizeOnly.get() and 'randomize_only' not in options:
                options.append('randomize_only')
            if len(options):
                for k in options: 
                    try:
                        line = "%s = %s\n" %(k,str(kw[k]))
                        print("1: line=", line)
                    except:
                        line = "%s = True\n" %(k)
                        print("2: line=", line)
                    ptr.write(line)
            ptr.close()
            if 'flex' in self.vina_keys:
                self.vina_keys.remove('flex')
    
        if not hasattr(self.vf,'vinaDict'): self.vf.vinaDict={}
        self.vf.vinaDict['configFile'] = outputconfigFile        
            

Dpf41CONFIGWriterGUI=CommandGUI()
Dpf41CONFIGWriterGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'], menuText['WriteCONFIG41'],\
cascadeName = menuText['WriteDpfMB'])



class DpfGAWriter(MVCommand):
    """ allows user to choose an output filename and write it"""


    def onAddCmdToViewer(self):
        checkHasDpo(self.vf)


    def __call__(self, outfile, **kw):
        self.doitWrapper(*(outfile,), **kw)


    def doit(self, outfile):
        if not len(self.vf.dpo.receptor_stem):
            self.vf.warningMsg("You must choose a macromolecule before writing dpf")
            return 'ERROR'

        #if a rms reference file has been specified, write it to dpf
        if self.vf.dpo['rmsref']['value']!=self.vf.dpo['move']['value']:
            l = genetic_algorithm_list
            ind = l.index('rmstol') + 1
            l.insert(ind, 'rmsref')

        self.vf.dpo.write(outfile, genetic_algorithm_list)
        

    def guiCallback(self):
        if not len(self.vf.dpo.receptor_stem):
            self.vf.warningMsg("You must choose a macromolecule before writing dpf")
            return 'ERROR'
        outfile = self.vf.askFileSave(types=[('dpf file', '*.dpf')],
                title = 'GA Docking Parameter Output File:')
        if outfile: 
            self.doitWrapper(outfile, log=1,redraw=0)


DpfGAWriterGUI=CommandGUI()
DpfGAWriterGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'], menuText['WriteGA'],\
cascadeName = menuText['WriteDpfMB'])



class DpfLSWriter(MVCommand):
    """ allows user to choose an output filename and write it"""


    def onAddCmdToViewer(self):
        checkHasDpo(self.vf)


    def __call__(self, outfile, **kw):
        self.doitWrapper(*(outfile,), **kw)


    def doit(self, outfile):
        if not len(self.vf.dpo.receptor_stem):
            self.vf.warningMsg("You must choose a macromolecule before writing dpf")
            return 'ERROR'

        #if a rms reference file has been specified, write it to dpf
        if self.vf.dpo['rmsref']['value']!=self.vf.dpo['move']['value']:
            l = local_search_list
            ind = l.index('rmstol') + 1
            l.insert(ind, 'rmsref')

        #self.vf.dpo.write(outfile, local_search_list)
        oldval = self.vf.dpo['set_sw1']['value']
        if oldval==0:
            self.vf.dpo['set_sw1']['value']=1
            self.vf.dpo['set_psw1']['value']=0
        self.vf.dpo.write(outfile, local_search_list)
        if oldval==0:
            self.vf.dpo['set_sw1']['value']=0
            self.vf.dpo['set_psw1']['value']=1


        
    def guiCallback(self):
        if not len(self.vf.dpo.receptor_stem):
            self.vf.warningMsg("You must choose a macromolecule before writing dpf")
            return 'ERROR'
        outfile = self.vf.askFileSave(types=[('dpf file', '*.dpf')],
                title = 'LS Docking Parameter Output File:')
        if outfile: 
            self.doitWrapper(outfile, log=1,redraw=0)


DpfLSWriterGUI=CommandGUI()
DpfLSWriterGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'], menuText['WriteLS'],\
cascadeName = menuText['WriteDpfMB'])



class DpfGALSWriter(MVCommand):
    """ allows user to choose an output filename and write it"""


    def onAddCmdToViewer(self):
        checkHasDpo(self.vf)


    def __call__(self, outfile, **kw):
        self.doitWrapper(*(outfile,), **kw)


    def doit(self, outfile):
        if not len(self.vf.dpo.receptor_stem):
            self.vf.warningMsg("You must choose a macromolecule before writing dpf")
            return 'ERROR'

        #if a rms reference file has been specified, write it to dpf
        if self.vf.dpo['rmsref']['value']!=self.vf.dpo['move']['value']:
            l = genetic_algorithm_local_search_list
            ind = l.index('rmstol') + 1
            l.insert(ind, 'rmsref')

        oldval = self.vf.dpo['set_sw1']['value']
        if oldval==0:
            self.vf.dpo['set_sw1']['value']=1
            self.vf.dpo['set_psw1']['value']=0
        self.vf.dpo.write(outfile, genetic_algorithm_local_search_list)
        if oldval==0:
            self.vf.dpo['set_sw1']['value']=0
            self.vf.dpo['set_psw1']['value']=1

        

    def guiCallback(self):
        if not len(self.vf.dpo.receptor_stem):
            self.vf.warningMsg("You must choose a macromolecule before writing dpf")
            return 'ERROR'
        outfile = self.vf.askFileSave(types=[('dpf file', '*.dpf')],
                title = ' GALS Docking Parameter Output File:')
        if outfile: 
            self.doitWrapper(outfile, log=1,redraw=0)


DpfGALSWriterGUI=CommandGUI()
DpfGALSWriterGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'],\
menuText['WriteGALS'], cascadeName = menuText['WriteDpfMB'], separatorAbove=1)



class DpfClusterWriter(MVCommand):
    """ allows user to write a dpf to cluster many dockings """


    def onAddCmdToViewer(self):
        checkHasDpo(self.vf)


    def __call__(self, outfile, **kw):
        self.doitWrapper(*(outfile,), **kw)


    def doit(self, outfile):
        if not len(self.vf.dpo.receptor_stem):
            self.vf.warningMsg("You must choose a macromolecule before writing dpf")
            return 'ERROR'
        self.vf.dpo.write(outfile, cluster_list)


    def guiCallback(self):
        if not len(self.vf.dpo.receptor_stem):
            self.vf.warningMsg("You must choose a macromolecule before writing dpf")
            return 'ERROR'
        outfile = self.vf.askFileSave(types=[('dpf file', '*.dpf')],
                title = ' Cluster Docking Parameter File Output:')
        if outfile: 
            self.doitWrapper(outfile, self.vf.dpo, log=1,redraw=0)


DpfClusterWriterGUI=CommandGUI()
DpfClusterWriterGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'],\
menuText['WriteCluster'], cascadeName = menuText['WriteDpfMB'])



class Dpf4SAWriter(MVCommand):
    """ allows user to choose an output filename and write simulated annealing parameters"""


    def onAddCmdToViewer(self):
        checkHasDpo(self.vf)


    def __call__(self, outfile, **kw):
        self.doitWrapper(*(outfile,), **kw)


    def doit(self, outfile):
        #set self.vf.dpo['ligand_types']['value'] here
        # from current values of self.vf.DPF_LIGAND_TYPES 
        #  AND
        # from current values of self.vf.DPF_FLEXRES_TYPES
        autodock_types = self.vf.DPF_LIGAND_TYPES[:]
        for t in self.vf.DPF_FLEXRES_TYPES:
            if t not in autodock_types:
                autodock_types.append(t)
        if len(autodock_types)==0:
            self.vf.warningMsg("You must choose a ligand before writing dpf")
            return 'ERROR'
        autodock_type_str = str(autodock_types[0])
        if len(autodock_types)>1:
            for t in autodock_types[1:]:
                autodock_type_str  = autodock_type_str + " " + str(t)
        self.vf.dpo['ligand_types']['value'] = autodock_type_str
        #to remove pickle problem, assume dpo is current self.vf.dpo
        if not len(self.vf.dpo.receptor_stem):
            self.vf.warningMsg("You must choose a macromolecule before writing dpf")
            return 'ERROR'

        #if a rms reference file has been specified, write it to dpf
        if self.vf.dpo['rmsref']['value']!=self.vf.dpo['move']['value']:
            l = simulated_annealing_list4
            ind = l.index('rmstol') + 1
            l.insert(ind, 'rmsref')
            
        self.vf.dpo.write4(outfile, simulated_annealing_list4)


    def guiCallback(self):
        if not len(self.vf.dpo.receptor_stem):
            self.vf.warningMsg("You must choose a macromolecule before writing dpf")
            return 'ERROR'
        outfile = self.vf.askFileSave(types=[('dpf file', '*.dpf')],
                title = 'AutoDock4 SA Docking Parameter Output File:')
        if outfile: 
            self.doitWrapper(outfile, log=1,redraw=0)


Dpf4SAWriterGUI=CommandGUI()
Dpf4SAWriterGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'],\
menuText['WriteSA4'], cascadeName = menuText['WriteDpfMB'])



class Dpf4GAWriter(MVCommand):
    """ allows user to choose an output filename and write it"""


    def onAddCmdToViewer(self):
        checkHasDpo(self.vf)


    def __call__(self, outfile, **kw):
        self.doitWrapper(*(outfile,), **kw)


    def doit(self, outfile):
        if not len(self.vf.dpo.receptor_stem):
            self.vf.warningMsg("You must choose a macromolecule before writing dpf")
            return 'ERROR'
        #set self.vf.dpo['ligand_types']['value'] here
        # from current values of self.vf.DPF_LIGAND_TYPES 
        #  AND
        # from current values of self.vf.DPF_FLEXRES_TYPES
        autodock_types = self.vf.DPF_LIGAND_TYPES[:]
        for t in self.vf.DPF_FLEXRES_TYPES:
            if t not in autodock_types:
                autodock_types.append(t)
        if len(autodock_types)==0:
            self.vf.warningMsg("You must choose a ligand before writing dpf")
            return 'ERROR'
        autodock_type_str = str(autodock_types[0])
        if len(autodock_types)>1:
            for t in autodock_types[1:]:
                autodock_type_str  = autodock_type_str + " " + str(t)
        self.vf.dpo['ligand_types']['value'] = autodock_type_str

        #if a rms reference file has been specified, write it to dpf
        if self.vf.dpo['rmsref']['value']!=self.vf.dpo['move']['value']:
            l = genetic_algorithm_list4
            ind = l.index('rmstol') + 1
            l.insert(ind, 'rmsref')

        self.vf.dpo.write4(outfile, genetic_algorithm_list4)
        #this is set when the dpo is written
        #self.vf.dpo.dpf_filename = outfile
        

    def guiCallback(self):
        if not len(self.vf.dpo.receptor_stem):
            self.vf.warningMsg("You must choose a macromolecule before writing dpf")
            return 'ERROR'
        outfile = self.vf.askFileSave(types=[('dpf file', '*.dpf')],
                title = 'AutoDock4 GA Docking Parameter Output File:')
        if outfile: 
            self.doitWrapper(outfile, log=1,redraw=0)


Dpf4GAWriterGUI=CommandGUI()
Dpf4GAWriterGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'], menuText['WriteGA4'],\
cascadeName = menuText['WriteDpfMB'])



class Dpf4LSWriter(MVCommand):
    """ allows user to choose an output filename and write it"""


    def onAddCmdToViewer(self):
        checkHasDpo(self.vf)


    def __call__(self, outfile, **kw):
        self.doitWrapper(*(outfile,), **kw)


    def doit(self, outfile):
        if not len(self.vf.dpo.receptor_stem):
            self.vf.warningMsg("You must choose a macromolecule before writing dpf")
            return 'ERROR'
        #set self.vf.dpo['ligand_types']['value'] here
        # from current values of self.vf.DPF_LIGAND_TYPES 
        #  AND
        # from current values of self.vf.DPF_FLEXRES_TYPES
        autodock_types = self.vf.DPF_LIGAND_TYPES[:]
        for t in self.vf.DPF_FLEXRES_TYPES:
            if t not in autodock_types:
                autodock_types.append(t)
        if len(autodock_types)==0:
            self.vf.warningMsg("You must choose a ligand before writing dpf")
            return 'ERROR'
        autodock_type_str = str(autodock_types[0])
        if len(autodock_types)>1:
            for t in autodock_types[1:]:
                autodock_type_str  = autodock_type_str + " " + str(t)
        self.vf.dpo['ligand_types']['value'] = autodock_type_str

        #if a rms reference file has been specified, write it to dpf
        if self.vf.dpo['rmsref']['value']!=self.vf.dpo['move']['value']:
            l = local_search_list4
            ind = l.index('rmstol') + 1
            l.insert(ind, 'rmsref')

        self.vf.dpo.write4(outfile, local_search_list4)
        #this is set when the dpo is written
        #self.vf.dpo.dpf_filename = outfile

        
    def guiCallback(self):
        if not len(self.vf.dpo.receptor_stem):
            self.vf.warningMsg("You must choose a macromolecule before writing dpf")
            return 'ERROR'
        outfile = self.vf.askFileSave(types=[('dpf file', '*.dpf')],
                title = 'AutoDock4 LS Docking Parameter Output File:')
        if outfile: 
            self.doitWrapper(outfile, log=1,redraw=0)


Dpf4LSWriterGUI=CommandGUI()
Dpf4LSWriterGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'], menuText['WriteLS4'],\
cascadeName = menuText['WriteDpfMB'])



class Dpf4GALSWriter(MVCommand):
    """ allows user to choose an output filename and write it"""


    def onAddCmdToViewer(self):
        checkHasDpo(self.vf)


    def __call__(self, outfile, **kw):
        self.doitWrapper(*(outfile,), **kw)


    def doit(self, outfile):
        if not len(self.vf.dpo.receptor_stem):
            self.vf.warningMsg("You must choose a macromolecule before writing dpf")
            return 'ERROR'
        #set self.vf.dpo['ligand_types']['value'] here
        # from current values of self.vf.DPF_LIGAND_TYPES 
        #  AND
        # from current values of self.vf.DPF_FLEXRES_TYPES
        autodock_types = self.vf.DPF_LIGAND_TYPES[:]
        for t in self.vf.DPF_FLEXRES_TYPES:
            if t not in autodock_types:
                autodock_types.append(t)
        if len(autodock_types)==0:
            self.vf.warningMsg("You must choose a ligand before writing dpf")
            return 'ERROR'
        autodock_type_str = str(autodock_types[0])
        if len(autodock_types)>1:
            for t in autodock_types[1:]:
                autodock_type_str  = autodock_type_str + " " + str(t)
        self.vf.dpo['ligand_types']['value'] = autodock_type_str
        #if a rms reference file has been specified, write it to dpf
        if self.vf.dpo['rmsref']['value']!=self.vf.dpo['move']['value']:
            l = genetic_algorithm_local_search_list4
            ind = l.index('rmstol') + 1
            l.insert(ind, 'rmsref')

        self.vf.dpo.write4(outfile, genetic_algorithm_local_search_list4)
        #this is set when the dpo is written
        #self.vf.dpo.dpf_filename = outfile
        

    def guiCallback(self):
        if not len(self.vf.dpo.receptor_stem):
            self.vf.warningMsg("You must choose a macromolecule before writing dpf")
            return 'ERROR'
        outfile = self.vf.askFileSave(types=[('dpf file', '*.dpf')],
                title = 'AutoDock4 GALS Docking Parameter Output File:')
        if outfile: 
            self.doitWrapper(outfile, log=1,redraw=0)


Dpf4GALSWriterGUI=CommandGUI()
Dpf4GALSWriterGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'],\
menuText['WriteGALS4'], cascadeName = menuText['WriteDpfMB'])



class Dpf41SAWriter(MVCommand):
    """ allows user to choose an output filename and write simulated annealing parameters"""


    def onAddCmdToViewer(self):
        checkHasDpo(self.vf)


    def __call__(self, outfile, **kw):
        self.doitWrapper(*(outfile,), **kw)


    def doit(self, outfile):
        #set self.vf.dpo['ligand_types']['value'] here
        # from current values of self.vf.DPF_LIGAND_TYPES 
        #  AND
        # from current values of self.vf.DPF_FLEXRES_TYPES
        autodock_types = self.vf.DPF_LIGAND_TYPES[:]
        for t in self.vf.DPF_FLEXRES_TYPES:
            if t not in autodock_types:
                autodock_types.append(t)
        if len(autodock_types)==0:
            self.vf.warningMsg("You must choose a ligand before writing dpf")
            return 'ERROR'
        autodock_type_str = str(autodock_types[0])
        if len(autodock_types)>1:
            for t in autodock_types[1:]:
                autodock_type_str  = autodock_type_str + " " + str(t)
        self.vf.dpo['ligand_types']['value'] = autodock_type_str
        #to remove pickle problem, assume dpo is current self.vf.dpo
        if not len(self.vf.dpo.receptor_stem):
            self.vf.warningMsg("You must choose a macromolecule before writing dpf")
            return 'ERROR'

        #if a rms reference file has been specified, write it to dpf
        if self.vf.dpo['rmsref']['value']!=self.vf.dpo['move']['value']:
            l = simulated_annealing_list4
            ind = l.index('rmstol') + 1
            l.insert(ind, 'rmsref')

        #if self.vf.dpo['parameter_file']['value']=="" or self.vf.dpo['parameter_file']['value']=="AD4_parameters.dat":
        if self.vf.dpo['parameter_file']['value']=="AD4_parameters.dat":
            self.vf.dpo['parameter_file']['value']="AD4.1_bound.dat"
            self.vf.dpo['custom_parameter_file']['value']= 1
        self.vf.dpo.write42(outfile, simulated_annealing_list4_2)


    def guiCallback(self):
        if not len(self.vf.dpo.receptor_stem):
            self.vf.warningMsg("You must choose a macromolecule before writing dpf")
            return 'ERROR'
        outfile = self.vf.askFileSave(types=[('dpf file', '*.dpf')],
                title = 'AutoDock4 SA Docking Parameter Output File:')
        if outfile: 
            self.doitWrapper(outfile, log=1,redraw=0)


Dpf41SAWriterGUI=CommandGUI()
Dpf41SAWriterGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'],\
menuText['WriteSA41'], cascadeName = menuText['WriteDpfMB'])



class Dpf41GAWriter(MVCommand):
    """ allows user to choose an output filename and write it"""


    def onAddCmdToViewer(self):
        checkHasDpo(self.vf)


    def __call__(self, outfile, **kw):
        self.doitWrapper(*(outfile,), **kw)


    def doit(self, outfile):
        if not len(self.vf.dpo.receptor_stem):
            self.vf.warningMsg("You must choose a macromolecule before writing dpf")
            return 'ERROR'
        #set self.vf.dpo['ligand_types']['value'] here
        # from current values of self.vf.DPF_LIGAND_TYPES 
        #  AND
        # from current values of self.vf.DPF_FLEXRES_TYPES
        autodock_types = self.vf.DPF_LIGAND_TYPES[:]
        for t in self.vf.DPF_FLEXRES_TYPES:
            if t not in autodock_types:
                autodock_types.append(t)
        if len(autodock_types)==0:
            self.vf.warningMsg("You must choose a ligand before writing dpf")
            return 'ERROR'
        autodock_type_str = str(autodock_types[0])
        if len(autodock_types)>1:
            for t in autodock_types[1:]:
                autodock_type_str  = autodock_type_str + " " + str(t)
        self.vf.dpo['ligand_types']['value'] = autodock_type_str

        #if a rms reference file has been specified, write it to dpf
        if self.vf.dpo['rmsref']['value']!=self.vf.dpo['move']['value']:
            l = genetic_algorithm_list4_2
            ind = l.index('rmstol') + 1
            l.insert(ind, 'rmsref')

        #if self.vf.dpo['parameter_file']['value']=="" or self.vf.dpo['parameter_file']['value']=="AD4_parameters.dat":
        if self.vf.dpo['parameter_file']['value']=="AD4_parameters.dat":
            self.vf.dpo['parameter_file']['value']="AD4.1_bound.dat"
            self.vf.dpo['custom_parameter_file']['value']= 1
        self.vf.dpo.write42(outfile, genetic_algorithm_list4_2)
        #this is set when the dpo is written
        #self.vf.dpo.dpf_filename = outfile
        

    def guiCallback(self):
        if not len(self.vf.dpo.receptor_stem):
            self.vf.warningMsg("You must choose a macromolecule before writing dpf")
            return 'ERROR'
        outfile = self.vf.askFileSave(types=[('dpf file', '*.dpf')],
                title = 'AutoDock4 GA Docking Parameter Output File:')
        if outfile: 
            self.doitWrapper(outfile, log=1,redraw=0)


Dpf41GAWriterGUI=CommandGUI()
Dpf41GAWriterGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'], menuText['WriteGA41'],\
cascadeName = menuText['WriteDpfMB'])



class Dpf41LSWriter(MVCommand):
    """ allows user to choose an output filename and write it"""


    def onAddCmdToViewer(self):
        checkHasDpo(self.vf)


    def __call__(self, outfile, **kw):
        self.doitWrapper(*(outfile,), **kw)


    def doit(self, outfile):
        if not len(self.vf.dpo.receptor_stem):
            self.vf.warningMsg("You must choose a macromolecule before writing dpf")
            return 'ERROR'
        #set self.vf.dpo['ligand_types']['value'] here
        # from current values of self.vf.DPF_LIGAND_TYPES 
        #  AND
        # from current values of self.vf.DPF_FLEXRES_TYPES
        autodock_types = self.vf.DPF_LIGAND_TYPES[:]
        for t in self.vf.DPF_FLEXRES_TYPES:
            if t not in autodock_types:
                autodock_types.append(t)
        if len(autodock_types)==0:
            self.vf.warningMsg("You must choose a ligand before writing dpf")
            return 'ERROR'
        autodock_type_str = str(autodock_types[0])
        if len(autodock_types)>1:
            for t in autodock_types[1:]:
                autodock_type_str  = autodock_type_str + " " + str(t)
        self.vf.dpo['ligand_types']['value'] = autodock_type_str

        #if a rms reference file has been specified, write it to dpf
        if self.vf.dpo['rmsref']['value']!=self.vf.dpo['move']['value']:
            l = local_search_list4
            ind = l.index('rmstol') + 1
            l.insert(ind, 'rmsref')

        if self.vf.dpo['parameter_file']['value']=="" or self.vf.dpo['parameter_file']['value']=="AD4_parameters.dat":
            self.vf.dpo['parameter_file']['value']="AD4.1_bound.dat"
            self.vf.dpo['custom_parameter_file']['value']= 1
        self.vf.dpo.write42(outfile, local_search_list4_2)
        #this is set when the dpo is written
        #self.vf.dpo.dpf_filename = outfile

        
    def guiCallback(self):
        if not len(self.vf.dpo.receptor_stem):
            self.vf.warningMsg("You must choose a macromolecule before writing dpf")
            return 'ERROR'
        outfile = self.vf.askFileSave(types=[('dpf file', '*.dpf')],
                title = 'AutoDock4 LS Docking Parameter Output File:')
        if outfile: 
            self.doitWrapper(outfile, log=1,redraw=0)


Dpf41LSWriterGUI=CommandGUI()
Dpf41LSWriterGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'], menuText['WriteLS41'],\
cascadeName = menuText['WriteDpfMB'])



class Dpf41GALSWriter(MVCommand):
    """ allows user to choose an output filename and write an autodock4.1 dpf to it"""


    def onAddCmdToViewer(self):
        checkHasDpo(self.vf)


    def __call__(self, outfile, **kw):
        self.doitWrapper(*(outfile,), **kw)


    def doit(self, outfile):
        if not len(self.vf.dpo.receptor_stem):
            self.vf.warningMsg("You must choose a macromolecule before writing dpf")
            return 'ERROR'
        autodock_types = self.vf.DPF_LIGAND_TYPES[:]
        for t in self.vf.DPF_FLEXRES_TYPES:
            if t not in autodock_types:
                autodock_types.append(t)
        if len(autodock_types)==0:
            self.vf.warningMsg("You must choose a ligand before writing dpf")
            return 'ERROR'
        autodock_type_str = str(autodock_types[0])
        if len(autodock_types)>1:
            for t in autodock_types[1:]:
                autodock_type_str  = autodock_type_str + " " + str(t)
        self.vf.dpo['ligand_types']['value'] = autodock_type_str

        #if a rms reference file has been specified, write it to dpf
        if self.vf.dpo['rmsref']['value']!=self.vf.dpo['move']['value']:
            l = genetic_algorithm_local_search_list4_1
            ind = l.index('rmstol') + 1
            l.insert(ind, 'rmsref')

        if self.vf.dpo['parameter_file']['value']=="" or self.vf.dpo['parameter_file']['value']=="AD4_parameters.dat":
            self.vf.dpo['parameter_file']['value']="AD4.1_bound.dat"
            self.vf.dpo['custom_parameter_file']['value']= 1
        self.vf.dpo.write42(outfile, genetic_algorithm_local_search_list4_2)
        #this is set when the dpo is written
        #self.vf.dpo.dpf_filename = outfile
        

    def guiCallback(self):
        if not len(self.vf.dpo.receptor_stem):
            self.vf.warningMsg("You must choose a macromolecule before writing dpf")
            return 'ERROR'
        outfile = self.vf.askFileSave(types=[('dpf file', '*.dpf')],
                title = 'AutoDock4.2 GALS Docking Parameter Output File:')
        if outfile: 
            self.doitWrapper(outfile, log=1,redraw=0)


Dpf41GALSWriterGUI=CommandGUI()
Dpf41GALSWriterGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'],\
menuText['WriteGALS41'], cascadeName = menuText['WriteDpfMB'])



class Dpf4ClusterWriter(MVCommand):
    """ allows user to write a dpf to cluster many dockings """


    def onAddCmdToViewer(self):
        checkHasDpo(self.vf)


    def __call__(self, outfile, **kw):
        self.doitWrapper(*(outfile,), **kw)


    def doit(self, outfile):
        if not len(self.vf.dpo.receptor_stem):
            self.vf.warningMsg("You must choose a macromolecule before writing dpf")
            return 'ERROR'
        self.vf.dpo.write4(outfile, cluster_list4)


    def guiCallback(self):
        if not len(self.vf.dpo.receptor_stem):
            self.vf.warningMsg("You must choose a macromolecule before writing dpf")
            return 'ERROR'
        outfile = self.vf.askFileSave(types=[('dpf file', '*.dpf')],
                title = ' AutoDock4 Cluster Docking Parameter File Output:')
        if outfile: 
            self.doitWrapper(outfile,  log=1,redraw=0)


Dpf4ClusterWriterGUI=CommandGUI()
Dpf4ClusterWriterGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'],\
menuText['WriteCluster4'], cascadeName = menuText['WriteDpfMB'])



class SimAnneal(MVCommand):
    """ allows user to set necessary parameters for simulated annealing-based autodock job"""

    def guiCallback(self):
        """called each time the 'set other options ' button is pressed"""
        if not hasattr(self, 'form'):
            self.buildForm()
            self.form = self.vf.getUserInput(self.ifd, modal=0, blocking=0)
            self.form.root.protocol('WM_DELETE_WINDOW',self.Close_cb)
        else:
            self.form.root.deiconify()

        itemList = ['accs', 'cycles', 'rejs', 'runs', 'linear_schedule', 'select']
        varList = [self.accs, self.cycles, self.rejs, self.runs,\
                    self.linear_schedule, self.select]
        for i in range(len(itemList)):
            varList[i].set(self.vf.dpo[itemList[i]]['value'])

        itemList2 = ['dihrf','quarf', 'rt0', 'rtrf', 'trnrf']
        varList2 = [self.dihrf, self.quarf, self.rt0, self.rtrf, self.trnrf]
        for i in range(len(itemList2)):
            varList2[i].set(str(self.vf.dpo[itemList2[i]]['value']))

        ##first set the ints:
        #for item in ['accs', 'cycles', 'rejs', 'runs']:
            ##setattr(self,item,self.vf.dpo[item]['value'])
            #exec('self.'+item+'.set(self.vf.dpo[\''+item+"\']['value'])")
#
        ##next set the floats:
        #for item in ['dihrf','quarf', 'rt0', 'rtrf', 'trnrf']:
            #exec('self.'+item+'.set(str(self.vf.dpo[\''+item+"\']['value']))")
            ##setattr(self,item,self.vf.dpo[item]['value'])
#
        ##next set the booleans:
        #for item in ['linear_schedule']:
            #exec('self.'+item+'.set(self.vf.dpo[\''+item+"\']['value'])")
            ##setattr(self,item,self.vf.dpo[item]['value'])
#
        ##next set the strings:
        #for item in ['select']:
            #exec('self.'+item+'.set(self.vf.dpo[\''+item+"\']['value'])")
            #setattr(self,item,self.vf.dpo[item]['value'])

        #self.getTraj()

    def buildForm(self):
        # reduction factor variables:
        self.trnrf = tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.quarf =  tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.dihrf = tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.rtrf = tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.rt0 = tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.linear_schedule= tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.runs= tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.cycles= tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.accs= tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.rejs= tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.select= tkinter.StringVar(master=self.vf.GUI.ROOT)

        ## variables for trajectory output:
        #self.trajVar = Tkinter.StringVar()
        #self.trajVar.set('0')
        #self.trjbeg =Tkinter.StringVar()
        #self.trjbeg.set('45')
        #self.trjend =Tkinter.StringVar()
        #self.trjend.set('50')
        #self.trjfrq =Tkinter.StringVar()
        #self.trjend.set('7500')
        #self.trjout =Tkinter.StringVar()
        #self.trjend.set('')
        #self.trjsel =Tkinter.StringVar()
        #self.trjsel.set('E')
        #self.watch=Tkinter.StringVar()
        #self.watch.set('')

        ifd = self.ifd = InputFormDescr(title = "Simulated Annealing Parameters")
        ifd.append( {'name': 'numLab',
            'widgetType':tkinter.Label,
            'text': 'NUMBER OF:',
            'gridcfg':{'sticky':tkinter.W+ tkinter.E, 'columnspan':6}})
        ifd.append( {'name': 'runsEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Runs:',
                'textvariable': self.runs
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':3}})
        ifd.append( {'name': 'accsEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Accepted steps/cycle:',
                'textvariable': self.accs
            },
            'gridcfg':{'sticky':tkinter.E,'row':-1, 'column':3, 'columnspan':3}})
        ifd.append( {'name': 'cyclesEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Cycles:',
                'textvariable': self.cycles
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':3}})
        ifd.append( {'name': 'rejsEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Rejected steps/cycle:',
                'textvariable': self.rejs
            },
            'gridcfg':{'sticky':tkinter.E,'row':-1, 'column':3, 'columnspan':3}})
        ifd.append({'widgetType': tkinter.Label,
             'text':'_______________________________________',
             'wcfg':{'bd':6},
             'gridcfg':{'sticky':tkinter.E+tkinter.W, 'columnspan':6}})
        ifd.append( {'name': 'nextcycleLab',
            'widgetType':tkinter.Label,
            'text': 'TO BEGIN NEXT CYCLE, USE:',
            'gridcfg':{'sticky':tkinter.W + tkinter.E, 'columnspan':6}})
        ifd.append( {'name': 'minLab',
            'widgetType':tkinter.Label,
            'text': 'minimum state:',
            #'gridcfg':{'sticky':Tkinter.E}})
            'gridcfg':{'sticky':tkinter.E, 'columnspan':2}})
        ifd.append({'name': 'minChoice',
            'widgetType':tkinter.Radiobutton,
            'wcfg': {'value':'m'},
            'variable': self.select,
            'gridcfg':{'sticky':tkinter.W, 'row':-1,'column':2}})
        ifd.append( {'name': 'lastLab',
            'widgetType':tkinter.Label,
            'text': 'last state:',
            'gridcfg':{'sticky':tkinter.E, 'row':-1, 'column':3}})
        ifd.append({'name': 'lastChoice',
            'widgetType':tkinter.Radiobutton,
            'wcfg': {'value':'l'},
            'variable': self.select,
            'gridcfg':{'sticky':tkinter.W, 'row':-1,'column':4}})
        ifd.append({'widgetType': tkinter.Label,
             'text':'_______________________________________',
             'wcfg':{'bd':6},
             'gridcfg':{'sticky':tkinter.E+tkinter.W, 'columnspan':6}})
        ifd.append( {'name': 'schLab',
            'widgetType':tkinter.Label,
            'text': 'REDUCTION SCHEDULE TYPE:',
            'gridcfg':{'sticky':tkinter.W + tkinter.E, 'columnspan':6}})
        ifd.append({'name': 'schLinChoice',
            'widgetType':tkinter.Radiobutton,
            'wcfg': {'value':1},
            'text': 'Linear',
            'variable': self.linear_schedule,
            'gridcfg':{'sticky':tkinter.W,'columnspan':3}})
        ifd.append({'name': 'schGeomChoice',
            'widgetType':tkinter.Radiobutton,
            'wcfg': {'value':0},
            'text': 'Geometric',
            'variable': self.linear_schedule,
            'gridcfg':{'sticky':tkinter.W, 'row':-1, 'column':3}})
        ifd.append({'widgetType': tkinter.Label,
             'text':'_______________________________________',
             'wcfg':{'bd':6},
             'gridcfg':{'sticky':tkinter.E+tkinter.W, 'columnspan':6}})
        ifd.append( {'name': 'rfpcLab',
            'widgetType':tkinter.Label,
            'text': 'REDUCTION FACTORS PER CYCLE:',
            'gridcfg':{'sticky':tkinter.W+tkinter.E, 'columnspan':6}})
        ifd.append( {'name': 'trnrfEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Translation',
                'textvariable': self.trnrf,
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':3}})
        ifd.append( {'name': 'quarfEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Quaternion',
                'textvariable': self.quarf
            },
            'gridcfg':{'sticky':tkinter.E, 'row':-1, 'column':3, 'columnspan':3}})
        ifd.append( {'name': 'dihrfEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Dihedral',
                'textvariable': self.dihrf
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':3}})
        ifd.append( {'name': 'rtrfEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Temperature',
                'textvariable': self.rtrf
            },
            #'textvariable': self.temprf,
            'gridcfg':{'sticky':tkinter.E, 'row':-1, 'column':3, 'columnspan':3}})
        ifd.append({'widgetType': tkinter.Label,
             'text':'_______________________________________',
             'wcfg':{'bd':6},
             'gridcfg':{'sticky':tkinter.E+tkinter.W, 'columnspan':6}})
        ifd.append( {'name': 'rt0Ent',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Initial Temperature(Degrees):',
                'textvariable': self.rt0
            },
            'gridcfg':{'sticky':tkinter.W +tkinter.E, 'columnspan':6}})
        #ifd.append({'widgetType': Tkinter.Label,
                #'text':'_______________________________________',
                #'wcfg':{'bd':6},
                #'gridcfg':{'sticky':Tkinter.E+Tkinter.W, 'columnspan':6}})
        #ifd.append({'widgetType': Tkinter.Label,
            #'text':'STEP SIZE:',
            #'wcfg':{'bd':6},
            #'gridcfg':{'sticky':Tkinter.W+Tkinter.E, 'columnspan':6}})
        #ifd.append( {'name': 'transStepEnt',
            #'widgetType':Tkinter.Entry,
            #'label': 'Translation (Angstrom/step):\nYou can enter values for first and last cycles\nto have AutoDock calculate trnrf',
            #'textvariable': self.tstep,
            #'gridcfg':{'sticky':Tkinter.E, 'columnspan':6}})
        #ifd.append( {'name': 'qStepEnt',
            #'widgetType':Tkinter.Entry,
            #'label': 'Quaternion (Degree/step):',
            #'textvariable': self.qstep,
            #'gridcfg':{'sticky':Tkinter.E, 'columnspan':6}})
        #ifd.append( {'name': 'dStepEnt',
            #'widgetType':Tkinter.Entry,
            #'label': 'Torsion (Degree/step):',
            #'textvariable': self.dstep,
            #'gridcfg':{'sticky':Tkinter.E, 'columnspan':6}})
        #ifd.append({'widgetType': Tkinter.Label,
            #'text':'_______________________________________',
            #'wcfg':{'bd':6},
            #'gridcfg':{'sticky':Tkinter.E+Tkinter.W, 'columnspan':6}})
        #form parts to manage optional trajectory output
        #ifd.append( {'name': 'trajectChoiceLab',
            #'widgetType':Tkinter.Label,
            #'text': 'OUTPUT TRAJECTORY ?',
            #'gridcfg':{'sticky':Tkinter.W + Tkinter.E, 'row':16,'column':0,'columnspan':6}})
        #ifd.append( {'name': 'trajYesLab',
            #'widgetType':Tkinter.Label,
            #'text': 'Yes:',
            #'gridcfg':{'sticky':Tkinter.E, 'columnspan':2}})
        #ifd.append({'name': 'trajYes',
            #'widgetType':Tkinter.Radiobutton,
            #'wcfg': {'value':'1'},
            #'variable': self.trajVar,
            #'gridcfg':{'sticky':Tkinter.W, 'row':-1,'column':2},
            #'command': self.getTraj })
        #ifd.append( {'name': 'trajNoLab',
            #'widgetType':Tkinter.Label,
            #'text': 'No: ',
            #'gridcfg':{'sticky':Tkinter.E, 'row': -1, 'column':3}})
        #ifd.append({'name': 'trajNo',
            #'widgetType':Tkinter.Radiobutton,
            #'wcfg': {'value':'0'},
            #'variable': self.trajVar,
            #'gridcfg':{'sticky':Tkinter.W, 'row':-1,'column':4},
            #'command': self.getTraj })
        # the next entries are governed by getTraj
        #ifd.append( {'name': 'trajbegLab',
            #'widgetType':Tkinter.Label,
            #'text': 'Begin at cycle number:',
            #'gridcfg':{'sticky':Tkinter.E, 'columnspan':3}})
        #ifd.append( {'name': 'trajbeg',
            #'widgetType':Tkinter.Entry,
            #'textvariable': self.trjbeg,
            #'gridcfg':{'sticky':Tkinter.W, 'row':-1, 'column':3, 'columnspan':3}})
        #ifd.append( {'name': 'trajendLab',
            #'widgetType':Tkinter.Label,
            #'text': 'End at cycle number:',
            #'gridcfg':{'sticky':Tkinter.E, 'columnspan':3}})
        #ifd.append( {'name': 'trajend',
            #'widgetType':Tkinter.Entry,
            #'textvariable': self.trjend,
            #'gridcfg':{'sticky':Tkinter.W, 'row':-1, 'column':3,'columnspan':3}})
        #ifd.append( {'name': 'trajfreqLab',
            #'widgetType':Tkinter.Label,
            #'text': 'Output Frequency :',
            #'gridcfg':{'sticky':Tkinter.E, 'columnspan':3}})
        #ifd.append( {'name': 'trajfreq',
            #'widgetType':Tkinter.Entry,
            #'textvariable': self.trjfrq,
            #'gridcfg':{'sticky':Tkinter.W, 'row':-1, 'column':3,'columnspan':3}})
        #ifd.append( {'name': 'trajectOutTypeLab',
            #'widgetType':Tkinter.Label,
            #'text': 'FOR TRAJECTORY,  OUTPUT:',
            #'gridcfg':{'sticky':Tkinter.W+Tkinter.E,'columnspan':6}})
        #ifd.append( {'name': 'trajacconlyLab',
            #'widgetType':Tkinter.Label,
            #'text': 'accepted steps only:',
            #'gridcfg':{'sticky':Tkinter.E, 'columnspan':2}})
        #ifd.append({'name': 'trajacconly',
            #'widgetType':Tkinter.Radiobutton,
            #'wcfg': {'value':'A'},
            #'variable': self.trjsel,
            #'gridcfg':{'sticky':Tkinter.W, 'row':-1,'column':2}})
        #ifd.append( {'name': 'trajaccrejLab',
            #'widgetType':Tkinter.Label,
            #'text': 'accepted + rejected steps:',
            #'gridcfg':{'sticky':Tkinter.E,'row':-1,'column':3,  'columnspan':2}})
        #ifd.append({'name': 'trajaccrej',
            #'widgetType':Tkinter.Radiobutton,
            #'wcfg': {'value':'E'},
            #'variable': self.trjsel,
            #'gridcfg':{'sticky':Tkinter.W, 'row':-1,'column':5}})
        #ifd.append( {'name': 'trajFileEntLab',
            #'widgetType':Tkinter.Label,
            #'text': 'Trajectory output filename:',
            #'gridcfg':{'sticky':Tkinter.W,'columnspan':2}})
        #ifd.append( {'name': 'trajFileEnt',
            #'widgetType':Tkinter.Entry,
            #'textvariable': self.trjout,
            #'gridcfg':{'sticky':Tkinter.E,'row': -1, 'column':2}})
        #ifd.append( {'name': 'watchFileEntLab',
            #'widgetType':Tkinter.Label,
            #'text': 'Watch filename:',
            #'gridcfg':{'sticky':Tkinter.E, 'row': -1, 'column':3,'columnspan':2}})
        #ifd.append( {'name': 'watchFileEnt',
            #'widgetType':Tkinter.Entry,
            #'textvariable': self.watch,
            #'gridcfg':{'sticky':Tkinter.W, 'row':-1, 'column':5}})
        ifd.append({'widgetType': tkinter.Button,
            'text':'Accept',
            'wcfg':{'bd':6},
            'gridcfg':{'sticky':tkinter.E+tkinter.W, 'columnspan':3, 'row':26, 'column':0},
            'command':self.Accept_cb})
        ifd.append({'widgetType': tkinter.Button,
            'text':'Close',
            'wcfg':{'bd':6},
            'gridcfg':{'sticky':tkinter.E+tkinter.W, 'columnspan':3, 'row':-1, 'column':3},
            'command':self.Close_cb})

    
    def Accept_cb(self, event=None):
        changeVals = {}
        for item in [ 'runs', 'accs', 'cycles', 'rejs', 'select',\
            'linear_schedule','trnrf', 'quarf', 'dihrf', 'rtrf', 'rt0']:
            var = eval('self.'+item)
            if self.vf.dpo[item]['value']!= var.get():
                changeVals[item] =  var.get()
        if len(list(changeVals.keys()))>0:
            changeVals['topCommand'] = 0
            self.doitWrapper(*(), **changeVals)
        self.form.withdraw()


    def doit(self, *args, **kw):
        self.vf.ADdpf_setDpo(*(), **kw)


    def Close_cb(self, event=None):
        self.form.withdraw()


    def getTraj(self):
        print('in getTraj')
        wlist = []
        for item in ['trajbegLab','trajbeg','trajendLab', 'trajend','trajfreqLab', 'trajfreq', 'trajectOutTypeLab','trajacconlyLab', 'trajacconly', 'trajaccrejLab', 'trajaccrej', 'trajFileEntLab', 'trajFileEnt','watchFileEntLab', 'watchFileEnt']:
            wlist.append(self.ifd.entryByName[item])
        if self.trajVar.get()=='0':
            for item in wlist:
                item['widget'].grid_forget()
            #possibly do something else here
        else:
            for item in wlist:
                item['widget'].grid(item['gridcfg'])


SimAnnealGUI=CommandGUI()
SimAnnealGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'],menuText['SA'], cascadeName = menuText['SetSearchParmsMB'])



class GA(MVCommand):
    """ allows user to set necessary parameters for genetic algorithm-based autodock job"""


    def guiCallback(self):
        """called each time the 'set other options ' button is pressed"""
        if not hasattr(self, 'form'):
            self.buildForm()
            self.form = self.vf.getUserInput(self.ifd, modal=0, blocking=0)
            self.form.root.protocol('WM_DELETE_WINDOW',self.Close_cb)
        else:
            self.form.root.deiconify()


        itemList = ['ga_run', 'ga_pop_size', 'ga_num_evals', 'ga_num_generations',\
                'ga_elitism', 'ga_window_size', 'do_global_only']
        varList = [self.ga_run, self.ga_pop_size, self.ga_num_evals,\
                   self.ga_num_generations, self.ga_elitism,\
                   self.ga_window_size, self.do_global_only]
        for i in range(len(itemList)):
            varList[i].set(self.vf.dpo[itemList[i]]['value'])

        itemList2 = [ 'ga_mutation_rate', 'ga_cauchy_alpha',\
            'ga_cauchy_beta','ga_crossover_rate']
        varList2 = [self.ga_mutation_rate, self.ga_cauchy_alpha,\
                    self.ga_cauchy_beta, self.ga_crossover_rate]
        for i in range(len(itemList2)):
            varList2[i].set(str(self.vf.dpo[itemList2[i]]['value']))

        ##first set the ints:
        #for item in ['ga_run', 'ga_pop_size', 'ga_num_evals', 'ga_num_generations',\
                #'ga_elitism', 'ga_window_size', 'do_global_only']:
            ##setattr(self,item,self.vf.dpo[item]['value'])
            #exec('self.'+item+'.set(self.vf.dpo[\''+item+"\']['value'])")

        ##next set the floats:
        #for item in ['ga_mutation_rate', 'ga_cauchy_alpha',\
            #'ga_cauchy_beta','ga_crossover_rate']:
            #exec('self.'+item+'.set(str(self.vf.dpo[\''+item+"\']['value']))")
            ##setattr(self,item,self.vf.dpo[item]['value'])


    def buildForm(self):
        #genetic algorithm variables:
        #ga_run is for GALS, do_global_only for GA
        self.ga_run=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.do_global_only=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.ga_pop_size=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.ga_num_evals=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.runType=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.runType.set('medium')
        self.runDict = {}
        self.runDict['short'] = '250000'
        self.runDict['medium'] = '2500000'
        self.runDict['long'] = '25000000'
        self.runTypeList = list(self.runDict.keys())
        self.crossover_modeType=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.crossover_modeType.set('twopt')
        self.crossover_modeTypeList = ['twopt', 'arithmetic', 'uniform']
        self.ga_num_generations=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.ga_elitism=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.ga_mutation_rate=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.ga_crossover_rate=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.ga_window_size=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.ga_cauchy_alpha=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.ga_cauchy_beta=tkinter.StringVar(master=self.vf.GUI.ROOT)
        ifd = self.ifd = InputFormDescr(title = "Genetic Algorithm Parameters")
        ifd.append( {'name': 'ga_runEnt',
            'widgetType':tkinter.Entry,
            'tooltip': 'Each run will result in one docked conformation.',
            'wcfg':{'label': 'Number of GA Runs:',
                'textvariable': self.ga_run,
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':6}})
        ifd.append( {'name': 'ga_pop_sizeEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Population Size:',
                'textvariable': self.ga_pop_size
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':6}})
        ifd.append({'widgetType':tkinter.Label,
                    'tooltip': 'Set maximum number of energy evaluations for each run.\nWhen the maximum is reached, the current generation\nwill finish but a new one will not be started.\nSelect type of run to set default number of evals...',
                    'wcfg': {'text':"Maximum Number of evals:"},
                    'gridcfg':{'sticky':'w'} })
        ifd.append({'widgetType':Pmw.ComboBox,
            'name':'run_type',
            'wcfg':{'entryfield_value':self.runType.get(),
                    'labelpos':'w',
                    'listheight':'60',
                    'entryfield_entry_width':7,
                    'scrolledlist_items': self.runTypeList,
                    'selectioncommand': self.setNumEvals,
                    },
            'gridcfg':{'sticky':'w', 'row':-1, 'column':1}})
        ifd.append( {'name': 'ga_num_evalsEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'textvariable': self.ga_num_evals
            },
            'gridcfg':{'sticky':tkinter.E, 'row':-1, 'column':2, 'columnspan':2}})
        ifd.append( {'name': 'ga_num_generationsEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Maximum Number of generations:',
                'textvariable': self.ga_num_generations
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':6}})
        ifd.append( {'name': 'ga_elitismEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Maximum Number of top individuals \nthat automatically survive:',
                'textvariable': self.ga_elitism
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':6}})
        ifd.append( {'name': 'ga_mutation_rateEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Rate of Gene Mutation:',
                'textvariable': self.ga_mutation_rate
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':6}})
        ifd.append( {'name': 'ga_crossover_rateEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Rate of Crossover:',
                'textvariable': self.ga_crossover_rate
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':6}})
        ifd.append({'name': 'ga_crossovermodeLabel',
                    'widgetType':tkinter.Label,
                    'tooltip': "Options for crossover mode are:\n'two_pt' (default for autodock3 and autodock4)\n'arithmetic' (default for autodock4.1)\n'uniform (experimental)'",
                    'wcfg': {'text':"GA Crossover mode:"},
                    'gridcfg':{'sticky':tkinter.E, 'columnspan':2 } })
        ifd.append({'name': 'ga_crossovermodeChoices',
                    'widgetType':Pmw.ComboBox,
                    'wcfg':{'entryfield_value':self.crossover_modeType.get(),
                    #'labelpos':'e',
                    'listheight':'60',
                    'entryfield_entry_width':17,
                    'scrolledlist_items': self.crossover_modeTypeList,
                    'selectioncommand': self.set_crossover_mode
                    },
            'gridcfg':{'sticky':'e', 'row':-1, 'column':1, 'columnspan':4}})
        ifd.append( {'name': 'ga_cauchy_alphaEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Mean of Cauchy distribution for\ngene mutation:',
                'textvariable': self.ga_cauchy_alpha
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':6}})
        ifd.append( {'name': 'ga_cauchy_betaEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Variance of Cauchy distribution for\ngene mutation:',
                'textvariable': self.ga_cauchy_beta
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':6}})
        ifd.append( {'name': 'ga_window_sizeEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Number of generations for picking\n worst individual:',
                'textvariable': self.ga_window_size
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':6}})
        ifd.append({'widgetType': tkinter.Button,
            'text':'Accept',
            'wcfg':{'bd':6},
            'gridcfg':{'sticky':tkinter.E+tkinter.W, 'columnspan':3, 'column':0},
            'command':self.Accept_cb})
        ifd.append({'widgetType': tkinter.Button,
            'text':'Close',
            'wcfg':{'bd':6},
            'gridcfg':{'sticky':tkinter.E+tkinter.W, 'column': 3, 'row':-1,'columnspan':3},
            'command':self.Close_cb})


    def setNumEvals(self, event=None):
        t = self.ifd.entryByName['run_type']['widget'].get()
        self.ga_num_evals.set(self.runDict[t])
        self.runType.set(t)


    def set_crossover_mode(self, event=None):
        t = self.ifd.entryByName['ga_crossovermodeChoices']['widget'].get()
        self.crossover_modeType.set(t)


    def Accept_cb(self, event=None):
        changeVals = {}
        for item in ['ga_run', 'ga_pop_size', 'ga_num_evals',\
            'ga_num_generations', 'ga_elitism', 'ga_mutation_rate', \
            'ga_crossover_rate', 'ga_cauchy_alpha', 'ga_cauchy_beta',\
            'ga_window_size']:
            var = eval('self.'+item)
            if self.vf.dpo[item]['value']!= var.get():
                changeVals[item] =  var.get()
        if self.crossover_modeType.get()!=self.vf.dpo['ga_crossover_mode']['value']:
            print("changing ga_crossover_mode from ", self.vf.dpo['ga_crossover_mode']['value'], ' to ', self.crossover_modeType.get())
            changeVals['ga_crossover_mode_flag']=1
            changeVals['ga_crossover_mode']= self.crossover_modeType.get()
        if len(list(changeVals.keys()))>0:
            changeVals['topCommand'] = 0
            self.doitWrapper(*(), **changeVals)
        self.form.withdraw()


    def Close_cb(self, event=None):
        self.form.withdraw()


    def doit(self, *args, **kw):
        self.vf.ADdpf_setDpo(*(), **kw)


GAGUI=CommandGUI()
GAGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'],menuText['GA'], cascadeName = menuText['SetSearchParmsMB'])



class LS(MVCommand):
    """ allows user to set necessary parameters for local search-based autodock job"""
    def guiCallback(self):
        """called each time the 'set other options ' button is pressed"""
        if not hasattr(self, 'form'):
            self.buildForm()
            self.form = self.vf.getUserInput(self.ifd, modal=0, blocking=0)
            self.form.root.protocol('WM_DELETE_WINDOW',self.Close_cb)
        else:
            self.form.root.deiconify()


        itemList = ['sw_max_its','sw_max_succ','sw_max_fail','do_local_only']
        varList = [self.sw_max_its, self.sw_max_succ, self.sw_max_fail,\
                    self.do_local_only]
        for i in range(len(itemList)):
            varList[i].set(self.vf.dpo[itemList[i]]['value'])

        itemList2 = ['sw_rho', 'sw_lb_rho', 'ls_search_freq', 'set_psw1','set_sw1']
        varList2 = [self.sw_rho, self.sw_lb_rho, self.ls_search_freq,\
                    self.set_psw1, self.set_sw1]
        for i in range(len(itemList2)):
            varList2[i].set(str(self.vf.dpo[itemList2[i]]['value']))

        ##first set the ints:
        #for item in ['sw_max_its','sw_max_succ','sw_max_fail','do_local_only']:
            ##setattr(self,item,self.vf.dpo[item]['value'])
            #exec('self.'+item+'.set(self.vf.dpo[\''+item+"\']['value'])")
            ##exec('self.'+item+'.set(self.vf.dpo[\''+item+"\']['value'])")
#
        ##next set the floats:
        #for item in ['sw_rho', 'sw_lb_rho', 'ls_search_freq']:
            #exec('self.'+item+'.set(str(self.vf.dpo[\''+item+"\']['value']))")
            ##exec('self.'+item+'.set(str(self.vf.dpo[\''+item+"\']['value']))")
            ##setattr(self,item,self.vf.dpo[item]['value'])
#
        ##last set the booleans:
        #for item in ['set_psw1','set_sw1']:
            #exec('self.'+item+'.set(str(self.vf.dpo[\''+item+"\']['value']))")
            ##exec('self.'+item+'.set(str(self.vf.dpo[\''+item+"\']['value']))")
            #setattr(self,item,self.vf.dpo[item]['value'])



    def buildForm(self):
        self.do_local_only=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.sw_max_its=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.sw_max_succ=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.sw_max_fail=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.sw_rho=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.sw_lb_rho=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.ls_search_freq=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.set_psw1=tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.set_sw1=tkinter.StringVar(master=self.vf.GUI.ROOT)
        ifd = self.ifd = InputFormDescr(title = "Local Search Parameters:")
        ifd.append( {'name': 'do_local_onlyEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Number of LS Runs:',
                'textvariable': self.do_local_only,
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':2}})
        ifd.append( {'name': 'sw_max_itsEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Maximum Number of iterations:',
                'textvariable': self.sw_max_its,
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':2}})
        ifd.append( {'name': 'sw_max_succEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Maximum Number of successes in a row\nbefore changing rho:',
                'textvariable': self.sw_max_succ,
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':2}})
        ifd.append( {'name': 'sw_max_failEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Maximum Number of failures in a row\nbefore changing rho:',
                'textvariable': self.sw_max_fail,
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':2}})
        ifd.append( {'name': 'sw_rhoEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Solis&Wets parameter defining initial variance\nand size of local space to sample (rho):',
                'textvariable': self.sw_rho,
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':2}})
        ifd.append( {'name': 'sw_lb_rhoEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Lower bound on rho:',
                'textvariable': self.sw_lb_rho,
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':2}})
        ifd.append( {'name': 'ls_search_freqEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'label': 'Probability of any particular phenotype being\nsubjected to local search:',
                'textvariable': self.ls_search_freq,
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':2}})
        ifd.append( {'name': 'lsChoiceLab',
            'widgetType':tkinter.Label,
            'text': 'FOR LOCAL SEARCH, USE: ',
            'gridcfg':{'sticky':tkinter.W + tkinter.E, 'columnspan':2}})
        ifd.append( {'name': 'swLab',
            'widgetType':tkinter.Label,
            'text': 'Solis & Wets with uniform variances:',
            'gridcfg':{'sticky':tkinter.E}})
        ifd.append({'name': 'swRb',
            'widgetType':tkinter.Radiobutton,
            'wcfg': {'value':0},
            'variable': self.set_psw1,
            'gridcfg':{'sticky':tkinter.W, 'row':-1,'column':1}})
        ifd.append( {'name': 'pswLab',
            'widgetType':tkinter.Label,
            'text': 'pseudo-Solis & Wets with relative variances:',
            'gridcfg':{'sticky':tkinter.E}})
        ifd.append({'name': 'psw',
            'widgetType':tkinter.Radiobutton,
            'wcfg': {'value':1},
            'variable': self.set_psw1,
            'gridcfg':{'sticky':tkinter.W, 'row':-1}})
        ifd.append({'name': 'acceptB',
            'widgetType': tkinter.Button,
            'text':'Accept',
            'wcfg':{'bd':4},
            'gridcfg':{'sticky':tkinter.E+tkinter.W, 'columnspan':3},
            'command':self.Accept_cb})
        ifd.append({'widgetType': tkinter.Button,
            'text':'Close',
            'wcfg':{'bd':6},
            'gridcfg':{'sticky':tkinter.E+tkinter.W, 'row':-1, 'column':3,'columnspan':3},
            'command':self.Close_cb})
    
    
    def Accept_cb(self, event=None):
        changeVals = {}
        for item in [ 'do_local_only', 'sw_max_its', 'sw_max_succ',\
            'sw_max_fail']: 
            var = eval('self.'+item)
            val = int(var.get())
            if self.vf.dpo[item]['value']!= val:
                changeVals[item] =  val
        for item in [ 'sw_rho', 'sw_lb_rho', 'ls_search_freq']:
            var = eval('self.'+item)
            val = float(var.get())
            if self.vf.dpo[item]['value']!= val:
                changeVals[item] =  val

        #have to deal with Boolean set_psw1 specially
        if self.set_psw1.get():
            if self.vf.dpo['set_psw1']['value']==0:
                changeVals['set_psw1']=1
                changeVals['set_sw1']=0
        else:
            if self.vf.dpo['set_psw1']['value']==1:
                changeVals['set_psw1']=0
                changeVals['set_sw1']=1

        if len(list(changeVals.keys()))>0:
            changeVals['topCommand'] = 0
            self.doitWrapper(*(), **changeVals)
        self.form.withdraw()


    def Close_cb(self, event=None):
        self.form.withdraw()


    def doit(self, *args, **kw):
        self.vf.ADdpf_setDpo(*(), **kw)


LSGUI=CommandGUI()
LSGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'],'Local Search Parameters ', cascadeName = menuText['SetSearchParmsMB'])



class SetDockingRunParms(MVCommand):
    """ allows user to set these parameters for  autodock job: step sizes, energy parameters and format for the output"""

    def guiCallback(self):
        """called each time the 'Set Docking Run Parameters' button is selected"""
        if not hasattr(self, 'form'):
            self.buildForm()
            #self.form = self.vf.getUserInput(self.ifd, modal=0, blocking=0)
            #self.form.root.protocol('WM_DELETE_WINDOW',self.Close_cb)
            #self.ranNumLib.set(1)
            #self.ranNumVar1.set(1)
            #self.ranNumVar2.set(2)
        else:
            self.form.root.deiconify()
        self.seed1.set(self.vf.dpo['seed']['value'][0])
        self.seed2.set(self.vf.dpo['seed']['value'][1])
        if self.seed1.get()=='time':
            self.ranNumVar1.set('1')
        elif self.seed1.get()=='pid':
            self.ranNumVar1.set('2')
        else:
            self.ranNumVar1.set('0')
        if self.seed2.get()=='time':
            self.ranNumVar2.set('1')
        elif self.seed2.get()=='pid':
            self.ranNumVar2.set('2')
        else:
            self.ranNumVar2.set('0')
        if self.seed2.get()=='':
            self.ranNumLib.set(1)
        else:
            self.ranNumLib.set(2)
        self.set_seeds()
        for item in ['userSeedLab1','userSeedEnt1','userSeedLab2','userSeedEnt2']:
            self.ifd.entryByName[item]['widget'].grid_forget()

        #??intelec, rmsnosym
        #first set the ints:

        itemList = ['outlev', 'analysis', 'write_all_flag','intelec','rmsref_flag'] 
        varList = [self.outlev, self.analysis, self.write_all_flag, self.intelec, self.rmsref_flag]
        for i in range(len(itemList)):
            varList[i].set(self.vf.dpo[itemList[i]]['value'])
        #update write_all_flag
        #self.vf.dpo['write_all_flag']['value'] = self.vf.dpo['write_all_flag']['value']
        itemList2 = ['extnrg', 'rmstol','dstep','qstep', 'rmsref']
        varList2 = [self.extnrg, self.rmstol, self.dstep, self.qstep, self.rmsref]
        for i in range(len(itemList2)):
            varList2[i].set(str(self.vf.dpo[itemList2[i]]['value']))

        ##first set the ints:
        #for item in ['outlev']:
            ##setattr(self,item,self.vf.dpo[item]['value'])
            #exec('self.'+item+'.set(self.vf.dpo[\''+item+"\']['value'])")

        #next set the floats:
        #for item in ['extnrg', 'rmstol','dstep','qstep']:
            #exec('self.'+item+'.set(str(self.vf.dpo[\''+item+"\']['value']))")
            ##setattr(self,item,self.vf.dpo[item]['value'])

        #next set the strings:
        #for item in ['rmsref']:
            #exec('self.'+item+'.set(str(self.vf.dpo[\''+item+"\']['value']))")
            ##setattr(self,item,self.vf.dpo[item]['value'])

        #next set the booleans:
        #for item in ['analysis', 'write_all_flag']:
            #exec('self.'+item+'.set(self.vf.dpo[\''+item+"\']['value'])")
            ##setattr(self,item,self.vf.dpo[item]['value'])

        #last fudge the lists:
        #oldval = self.vf.dpo['tstep']['value']
        #newval = ''
        #for item in oldval:
        #    newval = newval + str(item) + ','
        #newval = newval[:-1]
        newval = self.vf.dpo['tstep']['value']
        self.tstep.set(newval)
        oldval = self.vf.dpo['e0max']['value']
        self.e0max.set(str(oldval[0]))
        self.emaxRetries.set(str(oldval[1]))


    def buildForm(self):
        ifd = self.ifd = InputFormDescr(title = "Set Docking Run Options")
        #ranNumLib: 2 is platform independent one from UTexasBiomedicalSchool
        #ranNumLib: 1 is system's own implementation
        self.showLibOpts = tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.showEnergyOpts = tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.showStepSizeOpts = tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.showOutputOpts = tkinter.IntVar(master=self.vf.GUI.ROOT)
        for v in [self.showLibOpts, self.showEnergyOpts, self.showStepSizeOpts, 
                        self.showOutputOpts]:
            v.set(0)

        self.ranNumLib=tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.ranNumVar1=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.ranNumVar2=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.seed1=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.seed2=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.intelec=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.extnrg=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.e0max=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.emaxRetries=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.tstep = tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.qstep = tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.dstep = tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.outlev=tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.rmstol=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.rmsref=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.rmsref_flag=tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.rmsref_flag.set(0)
        self.rmsnosym=tkinter.StringVar(master=self.vf.GUI.ROOT)
        self.analysis=tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.write_all_flag=tkinter.IntVar(master=self.vf.GUI.ROOT)
        self.write_all_flag.set(0)
        ifd.append({'widgetType': tkinter.Label,
            'text':'for random number generator:',
            'wcfg':{'bd':1},
            'gridcfg':{'sticky':tkinter.W, 'columnspan':1}})
        ifd.append( {'name': 'LibOpts1',
            'widgetType':tkinter.Radiobutton,
            'variable': self.showLibOpts,
            'wcfg': {
                    'text': 'Use defaults',
                    'value':0,
                    'command': self.hideLibOpts,
                    },
            'gridcfg':{'row':-1, 'column':1, 'columnspan':2}})
            #'gridcfg':{'sticky':Tkinter.W, 'columnspan':2}})
        ifd.append( {'name': 'LibOpts2',
            'widgetType':tkinter.Radiobutton,
            'variable': self.showLibOpts,
            'wcfg': {
                    'text': 'Select library + set seeds',
                    'value':1,
                    'command': self.hideLibOpts,
                    },
            'gridcfg':{'row': -1, 'sticky':tkinter.W, 'column':3, 'columnspan':3}})
            
        ifd.append({'name': 'ranLibLabel',
            'widgetType': tkinter.Label,
            'text':'For RANDOM NUMBER GENERATOR LIBRARY:',
            'wcfg':{'bd':6},
            'gridcfg':{'sticky':tkinter.W, 'columnspan':6}})
        ifd.append( {'name': 'sysRanNumLibRB1',
            'widgetType':tkinter.Radiobutton,
            'variable': self.ranNumLib,
            'text': 'Built-In Library',
            'wcfg': {'value':1},
            'command': self.set_seeds,
            'gridcfg':{'sticky':tkinter.W}})
        ifd.append( {'name': 'indRanNumLibRB1',
            'widgetType':tkinter.Radiobutton,
            'variable': self.ranNumLib,
            'text': 'Platform-Independent Library\n(from UTexas Biomedical School):',
            'wcfg': {'value':2},
            'command': self.set_seeds,
            'gridcfg':{'sticky':tkinter.W,'row':-1,'column':1, 'columnspan':5}})
        ifd.append( {'name': 'ranNumChoiceLab',
            'widgetType':tkinter.Label,
            'text': 'SELECT ONE RANDOM NUMBER GENERATOR SEED:',
            'gridcfg':{'sticky':tkinter.W + tkinter.E, 'columnspan':6}})
        ifd.append({'name': 'time1',
            'widgetType':tkinter.Radiobutton,
            'wcfg': {'value':'1'},
            'text': 'time',
            'variable': self.ranNumVar1,
            'gridcfg':{'sticky':tkinter.W},
            'command': self.getUserSeed1 })
        ifd.append({'name': 'pid1',
            'widgetType':tkinter.Radiobutton,
            'wcfg': {'value':'2'},
            'text': 'pid',
            'variable': self.ranNumVar1,
            'gridcfg':{'sticky':tkinter.W, 'row':-1,'column':1, 'columnspan':2},
            'command': self.getUserSeed1 })
        ifd.append({'name': 'userSeedRb1',
            'widgetType':tkinter.Radiobutton,
            'wcfg': {'value':'0'},
            'text': 'user defined',
            'variable': self.ranNumVar1,
            'gridcfg':{'sticky':tkinter.W, 'row':-1,'columnspan':2, 'column':4},
            'command': self.getUserSeed1 })
        ifd.append({'name': 'userSeedLab1',
            'widgetType':tkinter.Label,
            'text': 'Enter Seed 1:',
            'gridcfg':{'sticky':tkinter.W+tkinter.E, 'columnspan':4}})
        ifd.append( {'name': 'userSeedEnt1',
            'widgetType':tkinter.Entry,
            'wcfg':{ 'textvariable': self.seed1},
            'gridcfg':{'sticky':tkinter.W, 'row':-1, 'column':4}})
        ifd.append({'name': 'time2',
            'widgetType':tkinter.Radiobutton,
            'wcfg': {'value':'1'},
            'text': 'time',
            'variable': self.ranNumVar2,
            'gridcfg':{'sticky':tkinter.W},
            'command': self.getUserSeed2 })
        ifd.append({'name': 'pid2',
            'widgetType':tkinter.Radiobutton,
            'text': 'pid',
            'wcfg': {'value':'2'},
            'variable': self.ranNumVar2,
            'gridcfg':{'sticky':tkinter.W, 'row':-1,'column':1,'columnspan':2},
            'command': self.getUserSeed2 })
        ifd.append({'name': 'userSeedRb2',
            'widgetType':tkinter.Radiobutton,
            'wcfg': {'value':'0'},
            'text': 'user defined',
            'variable': self.ranNumVar2,
            'gridcfg':{'sticky':tkinter.W, 'row':-1, 'column':4,'columnspan':2},
            'command': self.getUserSeed2 })
        ifd.append({'name': 'userSeedLab2',
            'widgetType':tkinter.Label,
            'text': 'Enter Seed 2:',
            'gridcfg':{'sticky':tkinter.W+tkinter.E, 'columnspan':4}})
        ifd.append( {'name': 'userSeedEnt2',
            'widgetType':tkinter.Entry,
            'wcfg':{'textvariable': self.seed2},
            'gridcfg':{'sticky':tkinter.W, 'row':-1, 'column':4}})
        ifd.append({'widgetType': tkinter.Label,
            'text':'_______________________________________',
            'wcfg':{'bd':1},
            'gridcfg':{'sticky':tkinter.E+tkinter.W, 'columnspan':6}})
        ifd.append({'widgetType': tkinter.Label,
            'text':'for energy parameters:',
            'wcfg':{'bd':1},
            'gridcfg':{'sticky':tkinter.W, 'columnspan':1}})
        ifd.append( {'name': 'EnergyOpts1',
            'widgetType':tkinter.Radiobutton,
            'variable': self.showEnergyOpts,
            'wcfg': {
                    'text': 'Use defaults',
                    'value':0,
                    'command': self.hideEnergyOpts,
                    },
            'gridcfg':{'row':-1, 'column':1, 'columnspan':2}})
            #'gridcfg':{'sticky':Tkinter.W, 'columnspan':2}})
        ifd.append( {'name': 'EnergyOpts2',
            'widgetType':tkinter.Radiobutton,
            'variable': self.showEnergyOpts,
            'wcfg': {
                    'text': 'Customize energy parameters',
                    'value':1,
                    'command': self.hideEnergyOpts,
                    },
            #'gridcfg':{'row': -1, 'column':3, 'columnspan':2}})
            'gridcfg':{'row': -1, 'sticky':tkinter.W, 'column':3, 'columnspan':3}})
            
        ifd.append({'name': 'energyLabel',
            'widgetType': tkinter.Label,
            'text':'ENERGY PARAMETERS:',
            'wcfg':{'bd':1},
            'gridcfg':{'sticky':tkinter.W+ tkinter.E, 'columnspan':6}})
        ifd.append( {'name': 'extnrgLab',
            'widgetType':tkinter.Label,
            'wcfg':{
                'text': 'External Grid Energy',
             },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':4}})
        ifd.append( {'name': 'extnrgEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'textvariable': self.extnrg
            },
            'gridcfg':{'sticky':tkinter.W,'row':-1, 'column':5, 'columnspan':2}})
        ifd.append( {'name': 'e0maxLab',
            'widgetType':tkinter.Label,
            'wcfg':{
                'text': 'Maximum allowable initial energy:',
             },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':4}})
        ifd.append( {'name': 'e0maxEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'textvariable': self.e0max
            },
            'gridcfg':{'sticky':tkinter.W,'row':-1, 'column':5, 'columnspan':2}})
            #'gridcfg':{'sticky':Tkinter.E, 'columnspan':6}})
        ifd.append( {'name': 'emaxRetriesLab',
            'widgetType':tkinter.Label,
            'wcfg':{
                'text': 'Maximum Number of Retries:',
             },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':4}})
        ifd.append( {'name': 'emaxRetriesEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'textvariable': self.emaxRetries
            },
            'gridcfg':{'sticky':tkinter.W,'row':-1, 'column':5, 'columnspan':2}})
            #'gridcfg':{'sticky':Tkinter.E, 'columnspan':6}})
        ifd.append( {'name': 'intelecLab1',
            'widgetType':tkinter.Label,
            'text': 'Calculate internal electrostatic energy:',
            'gridcfg':{'sticky':tkinter.E, 'columnspan':2}})
        ifd.append( {'name': 'intelecRB1',
            'widgetType':tkinter.Radiobutton,
            'variable': self.intelec,
            'text': 'Yes',
            'wcfg': {'value':'1'},
            'gridcfg':{'sticky':tkinter.E,'row':-1, 'columnspan':2,'column':2}})
        ifd.append( {'name': 'intelecRB0',
            'widgetType':tkinter.Radiobutton,
            'variable': self.intelec,
            'text': 'No',
            'wcfg': {'value':'0'},
            'gridcfg':{'sticky':tkinter.E + tkinter.W,'row':-1,'columnspan':2, 'column':4}})
        ifd.append({'widgetType': tkinter.Label,
            'text':'_______________________________________',
            'wcfg':{'bd':1},
            'gridcfg':{'sticky':tkinter.E+tkinter.W, 'columnspan':6}})
        ifd.append({'widgetType': tkinter.Label,
            'text':'for step size parameters:',
            'wcfg':{'bd':1},
            'gridcfg':{'sticky':tkinter.W, 'columnspan':1}})
        ifd.append( {'name': 'StepSizeOpts1',
            'widgetType':tkinter.Radiobutton,
            'variable': self.showStepSizeOpts,
            'wcfg': {
                    'text': 'Use defaults',
                    'value':0,
                    'command': self.hideStepSizeOpts,
                    },
            'gridcfg':{'row':-1, 'column':1, 'columnspan':2}})
            #'gridcfg':{'sticky':Tkinter.W, 'columnspan':2}})
        ifd.append( {'name': 'StepSizeOpts2',
            'widgetType':tkinter.Radiobutton,
            'variable': self.showStepSizeOpts,
            'wcfg': {
                    'text': 'Customize step size parameters',
                    'value':1,
                    'command': self.hideStepSizeOpts,
                    },
            #'gridcfg':{'row': -1, 'column':3, 'columnspan':2}})
            'gridcfg':{'row': -1, 'sticky':tkinter.W, 'column':3, 'columnspan':3}})
            
        ifd.append({'name':'stepSizeLab',
            'widgetType': tkinter.Label,
            'text':'STEP SIZE PARAMETERS:',
            'wcfg':{'bd':1},
            'gridcfg':{'sticky':tkinter.W+tkinter.E, 'columnspan':6}})
        ifd.append( {'name': 'tstepLab',
            'widgetType':tkinter.Label,
            'wcfg':{
                'text': 'Translation (Angstrom/step):\nEnter values for 1st , last cycles to have AutoDock calculate trnrf',
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':4}})
        ifd.append( {'name': 'tstepEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                #'label': 'Translation (Angstrom/step):\nEnter values for 1st , last cycles to have AutoDock calculate trnrf',
                'textvariable': self.tstep
            },
            'gridcfg':{'sticky':tkinter.E, 'row':-1, 'column':4}})
        ifd.append( {'name': 'qstepLab',
            'widgetType':tkinter.Label,
            'wcfg':{ 'text': 'Quaternion (Degree/step):',
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':3}})
        ifd.append( {'name': 'qstepEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'textvariable': self.qstep
            },
            'gridcfg':{'sticky':tkinter.E, 'row':-1, 'column':4}})
        ifd.append( {'name': 'dstepLab',
            'widgetType':tkinter.Label,
            'wcfg':{
                'text': 'Torsion (Degree/step):',
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':3}})
        ifd.append( {'name': 'dstepEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'textvariable': self.dstep,
            },
            'gridcfg':{'sticky':tkinter.E, 'row':-1, 'column':3, 'columnspan':3}})
        ifd.append({'widgetType': tkinter.Label,
            'text':'_______________________________________',
            'wcfg':{'bd':1},
            'gridcfg':{'sticky':tkinter.E+tkinter.W, 'columnspan':6}})
        ifd.append({'widgetType': tkinter.Label,
            'text':'for output format parameters:',
            'wcfg':{'bd':1},
            'gridcfg':{'sticky':tkinter.W, 'columnspan':1}})
        ifd.append( {'name': 'OutputOpt1',
            'widgetType':tkinter.Radiobutton,
            'variable': self.showOutputOpts,
            'wcfg': {
                    'text': 'Use defaults',
                    'value':0,
                    'command': self.hideOutputOpts,
                    },
            'gridcfg':{'row':-1, 'column':1, 'columnspan':2}})
            #'gridcfg':{'sticky':Tkinter.W, 'columnspan':2}})
        ifd.append( {'name': 'OutputOpt1',
            'widgetType':tkinter.Radiobutton,
            'variable': self.showOutputOpts,
            'wcfg': {
                    'text': 'Customize output format parameters',
                    'value':1,
                    'command': self.hideOutputOpts,
                    },
            #'gridcfg':{'row': -1, 'column':3, 'columnspan':2}})
            'gridcfg':{'row': -1, 'sticky':tkinter.W, 'column':3, 'columnspan':3}})
            
        ifd.append({'name':'outputOptsLab',
            'widgetType': tkinter.Label,
            'text':'OUTPUT FORMAT PARAMETERS:',
            'wcfg':{'bd':6},
            'gridcfg':{'sticky':tkinter.W+ tkinter.E, 'columnspan':6}})
            #'gridcfg':{'sticky':Tkinter.E+Tkinter.W, 'columnspan':6}})
        ifd.append( {'name': 'outlevLab',
            'widgetType':tkinter.Label,
            'text': 'Level of Detail for output:\n (use no output for GA and minimal for SA)',
            'gridcfg':{'sticky':tkinter.W+tkinter.E, 'columnspan':6}})
        ifd.append( {'name': 'outlevLab0',
            'widgetType':tkinter.Label,
            'text': 'no output:',
            'gridcfg':{'sticky':tkinter.E, 'columnspan':2}})
        ifd.append( {'name': 'outlevRB0',
            'widgetType':tkinter.Radiobutton,
            'variable': self.outlev,
            'wcfg': {'value':0},
            'gridcfg':{'sticky':tkinter.E,'row':-1, 'column':2}})
        ifd.append( {'name': 'outlevLab1',
            'widgetType':tkinter.Label,
            'text': 'minimal output:',
            'gridcfg':{'sticky':tkinter.E, 'row':-1, 'columnspan':2,'column':3}})
        ifd.append( {'name': 'outlevRB1',
            'widgetType':tkinter.Radiobutton,
            'variable': self.outlev,
            'wcfg': {'value':1},
            'gridcfg':{'sticky':tkinter.E,'row':-1, 'column':5}})
        ifd.append( {'name': 'outlevLab2',
            'widgetType':tkinter.Label,
            'text': 'full state output\nat end of each cycle:',
            'gridcfg':{'sticky':tkinter.E, 'columnspan':2}})
        ifd.append( {'name': 'outlevRB2',
            'widgetType':tkinter.Radiobutton,
            'variable': self.outlev,
            'wcfg': {'value':2},
            'gridcfg':{'sticky':tkinter.E,'row':-1, 'column':2}})
        ifd.append( {'name': 'outlevLab3',
            'widgetType':tkinter.Label,
            'text': 'detailed output \nfor each step:',
            'gridcfg':{'sticky':tkinter.E, 'row':-1,'columnspan':2,'column':3}})
        ifd.append( {'name': 'outlevRB3',
            'widgetType':tkinter.Radiobutton,
            'variable': self.outlev,
            'wcfg': {'value':3},
            'gridcfg':{'sticky':tkinter.E,'row':-1, 'column':5}})
        ifd.append( {'name': 'rmstolLab',
            'widgetType':tkinter.Label,
            'wcfg':{
                'text': 'RMS Cluster Tolerance (Angstrom):',
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':2}})
        ifd.append( {'name': 'rmstolEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'textvariable': self.rmstol,
            },
            'gridcfg':{'sticky':tkinter.E,'row':-1, 'column':2, 'columnspan':6}})
        ifd.append( {'name': 'rmsrefLab',
            'widgetType':tkinter.Label,
            'wcfg':{
                'text': 'Reference structure file for RMS calc:',
            },
            'gridcfg':{'sticky':tkinter.E, 'columnspan':2}})
        ifd.append( {'name': 'rmsrefEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{
                'textvariable': self.rmsref,
            },
            'gridcfg':{'sticky':tkinter.E,'row':-1, 'column':2, 'columnspan':6}})
        ifd.append( {'name': 'analysisChoiceLab',
            'widgetType':tkinter.Label,
            'text': 'for results of a docking: ',
            'gridcfg':{'sticky':tkinter.W +tkinter.E, 'columnspan':6}})
        ifd.append( {'name': 'analysisLab',
            'widgetType':tkinter.Label,
            'text': 'perform a cluster analysis:',
            'gridcfg':{'sticky':tkinter.E}})
        ifd.append({'name': 'analysisRB',
            'widgetType':tkinter.Radiobutton,
            'wcfg': {'value':1},
            'variable': self.analysis,
            'gridcfg':{'sticky':tkinter.W, 'row':-1,'column':1}})
        ifd.append( {'name': 'noanalysisLab',
            'widgetType':tkinter.Label,
            'text': 'do no analysis:',
            'gridcfg':{'sticky':tkinter.E,'row':-1,'columnspan':3, 'column':2}})
        ifd.append({'name': 'noanalysisRB',
            'widgetType':tkinter.Radiobutton,
            'wcfg': {'value':0},
            'variable': self.analysis,
            'gridcfg':{'sticky':tkinter.W, 'row':-1, 'column':5}})
        ifd.append( {'name': 'write_allLab',
            'widgetType':tkinter.Label,
            'text': '(for clustering multi-job output only:\nWrite all conformations in a cluster:',
            'gridcfg':{'sticky':tkinter.E}})
        ifd.append({'name': 'write_allRB',
            'widgetType':tkinter.Radiobutton,
            'wcfg': {'value':1},
            'variable': self.write_all_flag,
            'gridcfg':{'sticky':tkinter.W, 'row':-1,'column':1}})
        ifd.append( {'name': 'nowrite_allLab',
            'widgetType':tkinter.Label,
            'text': 'do not write all conformations:',
            'gridcfg':{'sticky':tkinter.E,'row':-1,'columnspan':3, 'column':2}})
        ifd.append({'name': 'nowrite_allRB',
            'widgetType':tkinter.Radiobutton,
            'wcfg': {'value':0},
            'variable': self.write_all_flag,
            'gridcfg':{'sticky':tkinter.W, 'row':-1, 'column':5}})
        ifd.append({'name': 'acceptB',
            'widgetType': tkinter.Button,
            'text':'Accept',
            'wcfg':{'bd':4},
            'gridcfg':{'sticky':tkinter.E+tkinter.W, 'columnspan':3},
            'command':self.Accept_cb})
        ifd.append({'widgetType': tkinter.Button,
            'text':'Close',
            'wcfg':{'bd':6},
            'gridcfg':{'sticky':tkinter.E+tkinter.W, 'row':-1, 'column':3,\
                'columnspan':3},
            'command':self.Close_cb})
        self.form = self.vf.getUserInput(self.ifd, modal=0, blocking=0)
        self.form.root.protocol('WM_DELETE_WINDOW', self.Close_cb)
        self.showLibOpts.set(0)
        self.hideLibOpts()
        self.hideEnergyOpts()
        self.hideStepSizeOpts()
        self.hideOutputOpts()
    

    def hideLibOpts(self, event=None):
        entry_names = [ 'ranLibLabel','sysRanNumLibRB1', 'indRanNumLibRB1', 
                        'ranNumChoiceLab', 'time1', 'pid1', 'userSeedRb1', 
                        'userSeedLab1','userSeedEnt1',
                        'userSeedLab2','userSeedEnt2', 
                        'time2', 'pid2', 'userSeedRb2']
        for n in entry_names:
            e = self.ifd.entryByName[n] 
            if not self.showLibOpts.get():
                e['widget'].grid_forget()
                self.seed1.set('time')
                self.seed2.set('pid')
                #restore the defaults random library seeds here
                self.vf.dpo['seed']['value'] = self.vf.dpo['seed']['default']

            else:
                if n not in ['userSeedLab1','userSeedEnt1','userSeedLab2','userSeedEnt2']:
                    e['widget'].grid(e['gridcfg'])
                if self.ranNumVar1.get()=='0':
                    w=self.ifd.entryByName['userSeedEnt1']
                    w['widget'].grid(w['gridcfg'])
                if self.ranNumVar2.get()=='0':
                    w=self.ifd.entryByName['userSeedEnt2']
                    w['widget'].grid(w['gridcfg'])
        self.form.autoSize()
    

    def hideEnergyOpts(self, event=None):
        entry_names = [ 'energyLabel', 'extnrgLab','extnrgEnt', 
                        'e0maxLab','e0maxEnt','emaxRetriesLab',
                        'emaxRetriesEnt', 'intelecLab1', 'intelecRB1',
                        'intelecRB0' ]
        for n in entry_names:
            e = self.ifd.entryByName[n] 
            if not self.showEnergyOpts.get():
                e['widget'].grid_forget()
                #restore defaults for intelec(0), extnrg(1000),e0max(0,10000)
                for p in ['intelec','extnrg','e0max']:
                    self.vf.dpo[p]['value'] = self.vf.dpo[p]['default']
                self.extnrg.set(str(self.vf.dpo['extnrg']['value']))
                self.intelec.set(str(self.vf.dpo['intelec']['value']))
                self.e0max.set(str(self.vf.dpo['e0max']['value'][0]))
                self.emaxRetries.set(str(int(self.vf.dpo['e0max']['value'][1])))
            else:
                e['widget'].grid(e['gridcfg'])
        self.form.autoSize()


    def hideStepSizeOpts(self, event=None):
        entry_names = [ 'stepSizeLab', 'tstepLab', 'tstepEnt',
                        'qstepLab', 'qstepEnt', 'dstepLab', 'dstepEnt', ]
        for n in entry_names:
            e = self.ifd.entryByName[n] 
            if not self.showStepSizeOpts.get():
                e['widget'].grid_forget()
                for p in ['dstep','qstep', 'tstep']:
                    self.vf.dpo[p]['value'] = self.vf.dpo[p]['default']
                self.dstep.set(str(self.vf.dpo['dstep']['value']))
                self.qstep.set(str(self.vf.dpo['qstep']['value']))
                self.tstep.set(str(self.vf.dpo['tstep']['value']))
                #tstep can have 3 values..(i've never seen it)
                #here val is probably [2.0]
                #val = self.tstep.get()
                #if val.find(',')<0:
                #    self.tstep.set(")
            else:
                e['widget'].grid(e['gridcfg'])
        self.form.autoSize()
    


    def hideOutputOpts(self, event=None):
        entry_names = [ 'outputOptsLab', 'outlevLab', 'outlevLab0',
                        'outlevRB0', 'outlevLab1', 'outlevRB1',
                        'outlevLab2', 'outlevRB2', 'outlevLab3',
                        'outlevRB3',  'rmstolLab', 'rmstolEnt', 
                        'rmsrefLab', 'rmsrefEnt', 'analysisChoiceLab', 
                        'analysisLab', 'analysisRB', 'noanalysisLab',
                        'noanalysisRB', 'write_allLab', 'write_allRB',
                        'nowrite_allLab', 'nowrite_allRB']
        for n in entry_names:
            e = self.ifd.entryByName[n] 
            if not self.showOutputOpts.get():
                e['widget'].grid_forget()
                for p in ['outlev','rmstol', 'rmsref_flag', 'analysis','write_all_flag']:
                    self.vf.dpo[p]['value'] = self.vf.dpo[p]['default']
                self.outlev.set(str(self.vf.dpo['outlev']['value']))
                self.rmstol.set(str(self.vf.dpo['rmstol']['value']))
                self.rmsref.set(str(self.vf.dpo['rmsref']['default']))
                self.analysis.set(str(self.vf.dpo['analysis']['default']))
                self.write_all_flag.set(self.vf.dpo['write_all_flag']['default'])
                self.rmsref_flag.set(self.vf.dpo['rmsref_flag']['default'])
            else:
                e['widget'].grid(e['gridcfg'])
        self.form.autoSize()
    
    

    def Accept_cb(self, event=None):
        changeVals = {}
        #first check the ints and booleans:
        for item in [ 'outlev', 'analysis', 'intelec', 'write_all_flag']:
            var = eval('self.'+item)
            val = int(var.get())
            if self.vf.dpo[item]['value']!= val:
                changeVals[item] =  val
        #next check the floats:
        for item in ['extnrg', 'rmstol','dstep','qstep']:
            var = eval('self.'+item)
            val = float(var.get())
            if self.vf.dpo[item]['value']!= val:
                changeVals[item] =  val
        #next check the only string->rmsref
        val = self.rmsref.get()
        #self.vf.dpo.ligand.name
        ligand_filename=None
        if hasattr(self.vf.dpo, 'ligand'):
            ligand_filename = os.path.basename(self.vf.dpo.ligand.parser.filename)
        if self.vf.dpo['rmsref']['value']!= val:
            if val!=ligand_filename:
                changeVals['rmsref'] =  val
                if val=="": 
                    changeVals['rmsref_flag'] = 0
                elif ligand_filename is not None and val==ligand_filename:
                    changeVals['rmsref'] =  ""
                    changeVals['rmsref_flag'] = 0
                else:
                    changeVals['rmsref_flag'] = 1
            else:
                changeVals['rmsref'] =  ""
                changeVals['rmsref_flag'] = 0
        #last special treatment for  the lists:
        vals = [self.seed1.get(), self.seed2.get()]
        if self.vf.dpo['seed']['value']!= vals:
            changeVals['seed'] =  vals
        #for item in ['e0max', 'tstep']:
        val = [float(self.e0max.get()), float(self.emaxRetries.get())]
        if self.vf.dpo['e0max']['value']!= val:
            changeVals['e0max'] =  val
        val = []
        if self.tstep.get().find(',')>-1:
            for item in string.split(self.tstep.get(), ','):
                val.append(float(item))
        elif self.tstep.get()[0]=='[':
            val = [float(self.tstep.get()[1:-1])]
        else:
            val = float(self.tstep.get())
            
        oldval = self.vf.dpo['tstep']['value']
        if val != oldval:
            changeVals['tstep'] = val

        if len(list(changeVals.keys()))>0:
            changeVals['topCommand'] = 0
            self.doitWrapper(*(), **changeVals)
        self.form.withdraw()


    def Close_cb(self, event=None):
        self.form.withdraw()


    def doit(self, *args, **kw):
        #print 'in docking run doit'
        self.vf.ADdpf_setDpo(*(), **kw)


    def set_seeds(self):
        if self.ranNumLib.get()==1:
            self.ifd.entryByName['ranNumChoiceLab']['widget'].config(text= 'SELECT ONE RANDOM NUMBER GENERATOR SEED:')
            for item in ['time2','pid2','userSeedRb2','userSeedLab2','userSeedEnt2']:
                self.ifd.entryByName[item]['widget'].grid_forget()
        elif self.showLibOpts.get(): 
            self.ifd.entryByName['ranNumChoiceLab']['widget'].config(text= 'SELECT TWO RANDOM NUMBER GENERATOR SEEDS:')
            for item in ['time2','pid2','userSeedRb2']:
                w=self.ifd.entryByName[item]
                w['widget'].grid(w['gridcfg'])
            if self.ranNumVar2.get()=='0':
                w=self.ifd.entryByName['userSeedEnt2']
                w['widget'].grid(w['gridcfg'])
        if self.seed1.get()==self.seed2.get():
            if self.ranNumVar1.get()=='2':
                self.ranNumVar2.set('1')
                self.seed2.set('time')
            elif self.ranNumVar1.get()=='1':
                self.ranNumVar2.set('2')
                self.seed2.set('pid')


    def getUserSeed1(self):
        w=self.ifd.entryByName['userSeedEnt1']
        w2=self.ifd.entryByName['userSeedLab1']
        val1= self.ranNumVar1.get()
        if val1 !='0':
            w2['widget'].grid_forget()
            w['widget'].grid_forget()
            #possibly do something else here
            if val1 =='1': 
                self.seed1.set('time')
                if self.seed1.get()==self.seed2.get():
                    self.ranNumVar2.set('2')
                    self.seed2.set('pid')
            else:
                self.seed1.set('pid')
                if self.seed1.get()==self.seed2.get():
                    self.ranNumVar2.set('1')
                    self.seed2.set('time')
        else:
            w2['widget'].grid(w2['gridcfg'])
            w['widget'].grid(w['gridcfg'])
            self.seed1.set('')

    def getUserSeed2(self):
        w=self.ifd.entryByName['userSeedEnt2']
        w2=self.ifd.entryByName['userSeedLab2']
        val2=self.ranNumVar2.get()
        if val2!='0':
            w2['widget'].grid_forget()
            w['widget'].grid_forget()
            #possibly do something else here
            if val2=='1':
                self.seed2.set('time')
                if self.ranNumVar2.get()==self.ranNumVar1.get():
                    self.ranNumVar1.set('2')
                    self.seed1.set('pid')
            else:
                self.seed2.set('pid')
                if self.ranNumVar2.get()==self.ranNumVar1.get():
                    self.ranNumVar1.set('1')
                    self.seed1.set('time')
        else:
            w2['widget'].grid(w2['gridcfg'])
            w['widget'].grid(w['gridcfg'])
            self.seed2.set('')


SetDockingRunParmsGUI=CommandGUI()
SetDockingRunParmsGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'],
    menuText['SetDockingRunParmsMB'])



class SetAutoDock4Parameters(MVCommand):
    """ allows user to set these parameters for autodock4 job: (1)parameter_file, (2)include_1_4_interactions (+float),  (3)unbound float, (4)epdb """


    def __init__(self):
        MVCommand.__init__(self)
        self.ifd_title = 'Set AutoDock4 Options'
        self.param_file = 'AD4_parameters.dat'

    def guiCallback(self):
        """called each time the 'Set AutoDock4 Parameters' button is selected"""
        if not hasattr(self, 'form'):
            self.buildForm()
            self.form = self.vf.getUserInput(self.ifd, modal=0, blocking=0)
            self.form.root.protocol('WM_DELETE_WINDOW',self.Close_cb)
            self.showParamFile_cb()
        else:
            self.form.root.deiconify()


    def buildForm(self):
        ifd = self.ifd = InputFormDescr(title = self.ifd_title)
        #ifd = self.ifd = InputFormDescr(title = "Set AutoDock4 Options")
        # for the new AD4 parameters
        self.pfile_flag = tkinter.IntVar(master=self.vf.GUI.ROOT)   #whether to include parameter_file keyword in dpf
        self.pfile_flag.set(0)
        self.parameter_file=tkinter.StringVar(master=self.vf.GUI.ROOT)                 #AD4 specific
        self.parameter_file.set("")
        #self.parameter_file.set(self.param_file)
        #self.parameter_file.set("AD4_parameters.dat")
        #self.include_1_4_interactions_flag = Tkinter.IntVar(master=self.vf.GUI.ROOT)   #AD4 specific
        #self.include_1_4_interactions = Tkinter.StringVar(master=self.vf.GUI.ROOT)     #AD4 specific
        #self.include_1_4_interactions.set("1.0")
        self.unbound = tkinter.StringVar(master=self.vf.GUI.ROOT)                      #AD4 specific
        self.unbound.set("0.0")
        self.epdb_flag = tkinter.IntVar(master=self.vf.GUI.ROOT)                       #AD4 specific
        self.epdb_flag.set(0)
        #self.epdb = Tkinter.StringVar(master=self.vf.GUI.ROOT)                         #AD4 specific
        self.compute_unbound_extended_flag = tkinter.IntVar(master=self.vf.GUI.ROOT)   #AD4 specific
        self.compute_unbound_extended_flag.set(1)
        self.rmsatoms = tkinter.StringVar(master=self.vf.GUI.ROOT)                     #AD4 specific
        self.rmsatoms.set('ligand_only')                        #AD4 specific
        #start AD4 specific:
        ifd.append({'name': 'pfileLab',
            'widgetType':tkinter.Label,
            'text': 'Include parameter_file in dpf:',
            'gridcfg':{'sticky':tkinter.W+tkinter.E, 'columnspan':2}})
        ifd.append( {'name': 'include_pfile_flag',
            'widgetType':tkinter.Radiobutton,
            'variable': self.pfile_flag,
            'text': 'Yes',
            'wcfg': {'value':'1'},
            'command': self.showParamFile_cb,
            'gridcfg':{'sticky':tkinter.W,'row':-1, 'columnspan':2,'column':2}})
        ifd.append( {'name': 'dont_include_pfile_flag',
            'widgetType':tkinter.Radiobutton,
            'variable': self.pfile_flag,
            'text': 'No',
            'command': self.showParamFile_cb,
            'wcfg': {'value':'0'},
            'gridcfg':{'sticky':tkinter.W,'row':-1,'columnspan':2, 'column':4}})
        ifd.append({'name': 'parameter_fileLab',
            'widgetType':tkinter.Label,
            'text': 'Enter parameter filename:',
            'gridcfg':{'sticky':tkinter.W+tkinter.E, 'columnspan':4}})
        ifd.append( {'name': 'parameter_fileEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{ 'textvariable': self.parameter_file},
            'gridcfg':{'sticky':tkinter.W, 'row':-1, 'column':4}})
        ifd.append( {'name': 'compute_unbound_extended_flagLab',
            'widgetType':tkinter.Label,
            'text': 'compute_unbound_extended?:',
            'gridcfg':{'sticky':tkinter.E, 'columnspan':2}})
        ifd.append( {'name': 'compute_unbound_extendedRB1',
            'widgetType':tkinter.Radiobutton,
            'variable': self.compute_unbound_extended_flag,
            'text': 'Yes',
            'wcfg': {'value':1},
            'gridcfg':{'sticky':tkinter.W,'row':-1, 'columnspan':2,'column':2}})
        ifd.append( {'name': 'compute_unbound_extendedRB0',
            'widgetType':tkinter.Radiobutton,
            'variable': self.compute_unbound_extended_flag,
            'text': 'No',
            'wcfg': {'value':0},
            'gridcfg':{'sticky':tkinter.W,'row':-1,'columnspan':2, 'column':4}}) 
        ifd.append({'name': 'unboundLab',
            'widgetType':tkinter.Label,
            'text': 'Enter unbound ligand energy:',
            'gridcfg':{'sticky':tkinter.W+tkinter.E, 'columnspan':4}})
        ifd.append( {'name': 'unboundEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{ 'textvariable': self.unbound},
            'gridcfg':{'sticky':tkinter.W, 'row':-1, 'column':4}})
        ifd.append( {'name': 'rmsatoms_lab',
            'widgetType':tkinter.Label,
            'text': 'for rms calculation use:',
            'gridcfg':{'sticky':tkinter.E, 'columnspan':2}})
        ifd.append( {'name': 'rmsatomsRB1',
            'widgetType':tkinter.Radiobutton,
            'variable': self.rmsatoms,
            'text': 'ligand_only',
            'wcfg': {'value':'ligand_only'},
            'gridcfg':{'sticky':tkinter.W,'row':-1, 'columnspan':2,'column':2}})
        ifd.append( {'name': 'rmsatomsRB0',
            'widgetType':tkinter.Radiobutton,
            'variable': self.rmsatoms,
            'text': 'ligand and flexible residues',
            'wcfg': {'value':'all'},
            'gridcfg':{'sticky':tkinter.W,'row':-1,'columnspan':2, 'column':4}}) 
        ifd.append( {'name': 'epdb_flagLab',
            'widgetType':tkinter.Label,
            'text': 'epdb?:',
            'gridcfg':{'sticky':tkinter.E, 'columnspan':2}})
        ifd.append( {'name': 'epdbRB1',
            'widgetType':tkinter.Radiobutton,
            'variable': self.epdb_flag,
            'text': 'Yes',
            'wcfg': {'value':'1'},
            'gridcfg':{'sticky':tkinter.W,'row':-1, 'columnspan':2,'column':2}})
        ifd.append( {'name': 'epdbRB0',
            'widgetType':tkinter.Radiobutton,
            'variable': self.epdb_flag,
            'text': 'No',
            'wcfg': {'value':'0'},
            'gridcfg':{'sticky':tkinter.W,'row':-1,'columnspan':2, 'column':4}})
        #ifd.append({'name': 'epdbLab',
        #    'widgetType':Tkinter.Label,
        #    'text': 'Enter filename for epdb calculation:',
        #    'gridcfg':{'sticky':Tkinter.W+Tkinter.E, 'columnspan':4}})
        #ifd.append( {'name': 'epdbFilename',
        #    'widgetType':Tkinter.Entry,
        #    'wcfg':{ 'textvariable': self.epdb},
        #    'gridcfg':{'sticky':Tkinter.W, 'row':-1, 'column':4}})
        #end AD4 specific:
        ifd.append({'name': 'acceptB',
            'widgetType': tkinter.Button,
            'text':'Accept',
            'wcfg':{'bd':4},
            'gridcfg':{'sticky':tkinter.E+tkinter.W, 'columnspan':3},
            'command':self.Accept_cb})
        ifd.append({'widgetType': tkinter.Button,
            'text':'Close',
            'wcfg':{'bd':6},
            'gridcfg':{'sticky':tkinter.E+tkinter.W, 'row':-1, 'column':3,\
                'columnspan':3},
            'command':self.Close_cb})
    
    
    def showParamFile_cb(self, event=None):
        #w=self.ifd.entryByName['barrierLab']
        pf_lab = self.ifd.entryByName['parameter_fileLab']
        pf_ent = self.ifd.entryByName['parameter_fileEnt']
        if self.pfile_flag.get()==0:
            self.parameter_file.set("")
            self.vf.dpo['custom_parameter_file']['value'] = 0
            self.vf.dpo['parameter_file']['value'] = ""
            pf_lab['widget'].config(state='disabled')
            pf_ent['widget'].config(state='disabled')
        else:
            pf_lab['widget'].config(state='normal')
            pf_ent['widget'].config(state='normal')
            #pf_lab['widget'].grid(pf_lab['gridcfg'])
            #pf_ent['widget'].grid(pf_ent['gridcfg'])
            self.vf.dpo['custom_parameter_file']['value'] = 1
            self.vf.dpo['parameter_file']['value'] = self.parameter_file.get()
        
    
    def Accept_cb(self, event=None):
        changeVals = {}
        #first check the ints and booleans:
        for item in [ 'epdb_flag', 'compute_unbound_extended_flag']:
            var = eval('self.'+item)
            val = var.get()
            if self.vf.dpo[item]['value']!= val:
                changeVals[item] =  val
        #next check the floats:
        for item in [ 'unbound']:
            var = eval('self.'+item)
            val = float(var.get())
            if self.vf.dpo[item]['value']!= val:
                changeVals[item] =  val
        #check the 2 new strings: parameter_file and epdb
        val = self.parameter_file.get()
        if self.vf.dpo['parameter_file']['value']!= val:
            if val=="":
                changeVals['custom_parameter_file'] =  0
            else:
                changeVals['custom_parameter_file'] =  1
            changeVals['parameter_file'] =  val
            
        #val = self.epdb.get()
        #if self.vf.dpo['epdb']['value']!= val:
        #    changeVals['epdb'] =  val
        val = self.rmsatoms.get()
        if self.vf.dpo['rmsatoms']['value']!= val:
            changeVals['rmsatoms_flag'] =  1
            changeVals['rmsatoms'] =  val
        if len(list(changeVals.keys()))>0:
            changeVals['topCommand'] = 0
            self.doitWrapper(*(), **changeVals)
        self.form.withdraw()

    def Close_cb(self, event=None):
        self.form.withdraw()


    def doit(self, *args, **kw):
        #print 'in set AutoDock4 parameters doit'
        self.vf.ADdpf_setDpo(*(), **kw)


SetAutoDock4ParametersGUI=CommandGUI()
SetAutoDock4ParametersGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'],
    menuText['SetAutoDock4Parameters'], cascadeName = menuText['OtherOptionsMB'])



class SetAutoDock41Parameters(MVCommand):
    """ allows user to set these parameters for autodock4 job: (1)parameter_file, (2)include_1_4_interactions (+float),  (3)unbound float, (4)epdb filename"""
    
    def __init__(self):
        MVCommand.__init__(self)
        self.ifd_title = 'Set AutoDock4.2 Options'
        self.param_file = 'AD4.1_bound.dat'


    def guiCallback(self):
        """called each time the 'Set AutoDock4.2 Parameters' button is selected"""
        if not hasattr(self, 'form'):
            self.buildForm()
            self.form = self.vf.getUserInput(self.ifd, modal=0, blocking=0)
            self.form.root.protocol('WM_DELETE_WINDOW',self.Close_cb)
            self.showParamFile_cb()
        else:
            self.form.root.deiconify()


    def buildForm(self):
        ifd = self.ifd = InputFormDescr(title = self.ifd_title)
        #ifd = self.ifd = InputFormDescr(title = "Set AutoDock4 Options")
        # for the new AD4 parameters
        self.pfile_flag = tkinter.IntVar(master=self.vf.GUI.ROOT)   #whether to include parameter_file keyword in dpf
        self.pfile_flag.set(0)
        self.parameter_file=tkinter.StringVar(master=self.vf.GUI.ROOT)                 #AD4 specific
        self.parameter_file.set("")
        #self.parameter_file.set(self.param_file)
        #self.parameter_file.set("AD4_parameters.dat")
        #self.include_1_4_interactions_flag = Tkinter.IntVar(master=self.vf.GUI.ROOT)   #AD4 specific
        #self.include_1_4_interactions = Tkinter.StringVar(master=self.vf.GUI.ROOT)     #AD4 specific
        #self.include_1_4_interactions.set("1.0")
        self.unbound = tkinter.StringVar(master=self.vf.GUI.ROOT)                      #AD4 specific
        self.unbound.set("0.0")
        self.epdb_flag = tkinter.IntVar(master=self.vf.GUI.ROOT)                       #AD4 specific
        self.epdb_flag.set(0)
        #self.epdb = Tkinter.StringVar(master=self.vf.GUI.ROOT)                         #AD4 specific
        #self.compute_unbound_extended_flag = Tkinter.IntVar(master=self.vf.GUI.ROOT)   #AD4 specific
        #self.compute_unbound_extended_flag.set(1)
        self.rmsatoms = tkinter.StringVar(master=self.vf.GUI.ROOT)                     #AD4 specific
        self.rmsatoms.set('ligand_only')                        #AD4 specific
        #start AD4.1 specific:
        ifd.append({'name': 'pfileLab',
            'widgetType':tkinter.Label,
            'text': 'Include parameter_file in dpf:',
            'gridcfg':{'sticky':tkinter.W+tkinter.E, 'columnspan':2}})
        ifd.append( {'name': 'include_pfile_flag',
            'widgetType':tkinter.Radiobutton,
            'variable': self.pfile_flag,
            'text': 'Yes',
            'wcfg': {'value':'1'},
            'command': self.showParamFile_cb,
            'gridcfg':{'sticky':tkinter.W,'row':-1, 'columnspan':2,'column':2}})
        ifd.append( {'name': 'dont_include_pfile_flag',
            'widgetType':tkinter.Radiobutton,
            'variable': self.pfile_flag,
            'text': 'No',
            'command': self.showParamFile_cb,
            'wcfg': {'value':'0'},
            'gridcfg':{'sticky':tkinter.W,'row':-1,'columnspan':2, 'column':4}})
        ifd.append({'name': 'parameter_fileLab',
            'widgetType':tkinter.Label,
            'text': 'Enter parameter_file:',
            'gridcfg':{'sticky':tkinter.W+tkinter.E, 'columnspan':4}})
        ifd.append( {'name': 'parameter_fileEnt',
            'widgetType':tkinter.Entry,
            'wcfg':{ 'textvariable': self.parameter_file},
            'gridcfg':{'sticky':tkinter.W, 'row':-1, 'column':4}})
        #ifd.append( {'name': 'include_1_4_Lab1',
        #    'widgetType':Tkinter.Label,
        #    'text': 'Include 1_4 interactions:',
        #    'gridcfg':{'sticky':Tkinter.E, 'columnspan':2}})
        #ifd.append( {'name': 'include_1_4_RB1',
        #    'widgetType':Tkinter.Radiobutton,
        #    'variable': self.include_1_4_interactions_flag,
        #    'text': 'Yes',
        #    'wcfg': {'value':'1'},
        #    'gridcfg':{'sticky':Tkinter.W,'row':-1, 'columnspan':2,'column':2}})
        #ifd.append( {'name': 'include_1_4_RB0',
        #    'widgetType':Tkinter.Radiobutton,
        #    'variable': self.include_1_4_interactions_flag,
        #    'text': 'No',
        #    'wcfg': {'value':'0'},
        #    'gridcfg':{'sticky':Tkinter.W,'row':-1,'columnspan':2, 'column':4}})
        #ifd.append({'name': 'include_1_4_Lab2',
        #    'widgetType':Tkinter.Label,
        #    'text': 'Enter scaling factor for 1_4_interactions:',
        #    'gridcfg':{'sticky':Tkinter.W+Tkinter.E, 'columnspan':4}})
        #ifd.append( {'name': 'include_1_4_value',
        #    'widgetType':Tkinter.Entry,
        #    'wcfg':{ 'textvariable': self.include_1_4_interactions},
        #    'gridcfg':{'sticky':Tkinter.W, 'row':-1, 'column':4}})
        ifd.append( {'name': 'rmsatoms_lab',
            'widgetType':tkinter.Label,
            'text': 'for rms calculation use:',
            'gridcfg':{'sticky':tkinter.E, 'columnspan':2}})
        ifd.append( {'name': 'rmsatomsRB1',
            'widgetType':tkinter.Radiobutton,
            'variable': self.rmsatoms,
            'text': 'ligand_only',
            'wcfg': {'value':'ligand_only'},
            'gridcfg':{'sticky':tkinter.W,'row':-1, 'columnspan':2,'column':2}})
        ifd.append( {'name': 'rmsatomsRB0',
            'widgetType':tkinter.Radiobutton,
            'variable': self.rmsatoms,
            'text': 'ligand and flexible residues',
            'wcfg': {'value':'all'},
            'gridcfg':{'sticky':tkinter.W,'row':-1,'columnspan':2, 'column':4}}) 
        ifd.append( {'name': 'epdb_flagLab',
            'widgetType':tkinter.Label,
            'text': 'epdb?:',
            'gridcfg':{'sticky':tkinter.E, 'columnspan':2}})
        ifd.append( {'name': 'epdbRB1',
            'widgetType':tkinter.Radiobutton,
            'variable': self.epdb_flag,
            'text': 'Yes',
            'wcfg': {'value':'1'},
            'gridcfg':{'sticky':tkinter.W,'row':-1, 'columnspan':2,'column':2}})
        ifd.append( {'name': 'epdbRB0',
            'widgetType':tkinter.Radiobutton,
            'variable': self.epdb_flag,
            'text': 'No',
            'wcfg': {'value':'0'},
            'gridcfg':{'sticky':tkinter.W,'row':-1,'columnspan':2, 'column':4}})
        #ifd.append({'name': 'epdbLab',
        #    'widgetType':Tkinter.Label,
        #    'text': 'Enter filename for epdb calculation:',
        #    'gridcfg':{'sticky':Tkinter.W+Tkinter.E, 'columnspan':4}})
        #ifd.append( {'name': 'epdbFilename',
        #    'widgetType':Tkinter.Entry,
        #    'wcfg':{ 'textvariable': self.epdb},
        #    'gridcfg':{'sticky':Tkinter.W, 'row':-1, 'column':4}})
        #end AD4 specific:
        ifd.append({'name': 'acceptB',
            'widgetType': tkinter.Button,
            'text':'Accept',
            'wcfg':{'bd':4},
            'gridcfg':{'sticky':tkinter.E+tkinter.W, 'columnspan':3},
            'command':self.Accept_cb})
        ifd.append({'widgetType': tkinter.Button,
            'text':'Close',
            'wcfg':{'bd':6},
            'gridcfg':{'sticky':tkinter.E+tkinter.W, 'row':-1, 'column':3,\
                'columnspan':3},
            'command':self.Close_cb})
    

    def showParamFile_cb(self, event=None):
        pf_lab = self.ifd.entryByName['parameter_fileLab']
        pf_ent = self.ifd.entryByName['parameter_fileEnt']
        if self.pfile_flag.get()==0:
            self.parameter_file.set("")
            self.vf.dpo['custom_parameter_file']['value'] = 0
            self.vf.dpo['parameter_file']['value'] = ""
            pf_lab['widget'].config(state='disabled')
            pf_ent['widget'].config(state='disabled')
        else:
            pf_lab['widget'].config(state='normal')
            pf_ent['widget'].config(state='normal')
            self.vf.dpo['custom_parameter_file']['value'] = 1
            self.vf.dpo['parameter_file']['value'] = self.parameter_file.get()

    
    def Accept_cb(self, event=None):
        changeVals = {}
        #first check the ints and booleans:
        for item in [ 'epdb_flag' ]:
            var = eval('self.'+item)
            val = var.get()
            if self.vf.dpo[item]['value']!= val:
                changeVals[item] =  val
        #check the 2 new strings: parameter_file and epdb
        val = self.parameter_file.get()
        if self.vf.dpo['parameter_file']['value']!= val:
            if val=="":
                changeVals['custom_parameter_file'] =  0
            else:
                changeVals['custom_parameter_file'] =  1
            changeVals['parameter_file'] =  val
            
        #val = self.epdb.get()
        #if self.vf.dpo['epdb']['value']!= val:
        #    changeVals['epdb'] =  val
        val = self.rmsatoms.get()
        if self.vf.dpo['rmsatoms']['value']!= val:
            changeVals['rmsatoms_flag'] =  1
            changeVals['rmsatoms'] =  val
        if len(list(changeVals.keys()))>0:
            changeVals['topCommand'] = 0
            self.doitWrapper(*(), **changeVals)
        self.form.withdraw()


    def Close_cb(self, event=None):
        self.form.withdraw()


    def doit(self, *args, **kw):
        self.vf.ADdpf_setDpo(*(), **kw)


SetAutoDock41ParametersGUI=CommandGUI()
SetAutoDock41ParametersGUI.addMenuCommand('AutoToolsBar', menuText['AutoDpfMB'],
    menuText['SetAutoDock41Parameters'], cascadeName = menuText['OtherOptionsMB'])



class StopAutoDpf(MVCommand):
    """hides the autotools menubar and """

    def __call__(self, **kw):
        self.doitWrapper(*(), **kw)

    def doit(self):
        if hasattr(self.vf,'ADtors_init'):
            self.vf.ADtors_init.rootSph.Set(visible=0)
        self.vf.GUI.VIEWER.Redraw()
        self.vf.GUI.menuBars['AutoToolsBar'].forget()
        
    def guiCallback(self):
        self.doitWrapper(log=1,redraw=0)



commandList = [
    {'name':'ADdpf_read','cmd':DpfLoadDefaults(),'gui':DpfLoadDefaultsGUI},
    {'name':'ADdpf_setDpo','cmd':DpfSetDpo(),'gui':None},
    #{'name':'ADdpf4_chooseMacromolecule','cmd':Dpf4MacroChooser(),'gui':Dpf4MacroChooserGUI},
    #macromolecule
    {'name':'ADdpf4_readMacromolecule','cmd':Dpf4MacroSelector(),'gui':Dpf4MacroSelectorGUI},
    {'name':'ADdpf4_readFlexRes','cmd':Dpf4FlexResSelector(),'gui':Dpf4FlexResSelectorGUI},
    {'name':'ADdpf_chooseMacromolecule','cmd':DpfMacroChooser(),'gui':DpfMacroChooserGUI},
    {'name':'ADdpf_readMacromolecule','cmd':DpfMacroSelector(),'gui':DpfMacroSelectorGUI},
    #ligand
    {'name':'ADdpf4_chooseFormattedLigand','cmd':Dpf4LigandChooser(),'gui':Dpf4LigandChooserGUI},
    {'name':'ADdpf4_readFormattedLigand','cmd':Dpf4LigPDBQReader(),'gui':Dpf4LigPDBQReaderGUI},
    {'name':'ADdpf4_initLigand','cmd':Dpf4InitLigand(),'gui':Dpf4InitLigandGUI},
    {'name':'ADdpf_chooseFormattedLigand','cmd':DpfLigandChooser(),'gui':DpfLigandChooserGUI},
    {'name':'ADdpf_readFormattedLigand','cmd':DpfLigPDBQReader(),'gui':DpfLigPDBQReaderGUI},
    {'name':'ADdpf_initLigand','cmd':DpfInitLigand(),'gui':DpfInitLigandGUI},

    {'name':'ADdpf_setGAparameters','cmd':GA(),'gui':GAGUI},
    {'name':'ADdpf_setSAparameters','cmd':SimAnneal(),'gui':SimAnnealGUI},
    {'name':'ADdpf_setLSparameters','cmd':LS(),'gui':LSGUI},
    {'name':'ADdpf_setDockingParameters','cmd':SetDockingRunParms(),'gui':SetDockingRunParmsGUI},
    {'name':'ADdpf_setAutoDock4Parameters','cmd':SetAutoDock4Parameters(),'gui':SetAutoDock4ParametersGUI},
    {'name':'ADdpf_setAutoDock41Parameters','cmd':SetAutoDock41Parameters(),'gui':SetAutoDock41ParametersGUI},
    {'name':'ADdpf4_writeGALS','cmd':Dpf4GALSWriter(),'gui':Dpf4GALSWriterGUI},
    {'name':'ADdpf41_writeGALS','cmd':Dpf41GALSWriter(),'gui':Dpf41GALSWriterGUI},
    {'name':'ADdpf4_writeGA','cmd':Dpf4GAWriter(),'gui':Dpf4GAWriterGUI},
    {'name':'ADdpf41_writeGA','cmd':Dpf41GAWriter(),'gui':Dpf41GAWriterGUI},
    {'name':'ADdpf4_writeSA','cmd':Dpf4SAWriter(),'gui':Dpf4SAWriterGUI},
    {'name':'ADdpf41_writeSA','cmd':Dpf41SAWriter(),'gui':Dpf41SAWriterGUI},
    {'name':'ADdpf4_writeLS','cmd':Dpf4LSWriter(),'gui':Dpf4LSWriterGUI},
    {'name':'ADdpf41_writeLS','cmd':Dpf41LSWriter(),'gui':Dpf41LSWriterGUI},
    {'name':'ADdpf41_writeConfig','cmd':Dpf41CONFIGWriter(),'gui':Dpf41CONFIGWriterGUI},
    {'name':'ADdpf41_writeEPDB','cmd':Dpf41EPDBWriter(),'gui':Dpf41EPDBWriterGUI},
    #{'name':'ADdpf4_writeCluster','cmd':Dpf4ClusterWriter(),'gui':Dpf4ClusterWriterGUI},
    {'name':'ADdpf_writeGALS','cmd':DpfGALSWriter(),'gui':DpfGALSWriterGUI},
    {'name':'ADdpf_writeGA','cmd':DpfGAWriter(),'gui':DpfGAWriterGUI},
    {'name':'ADdpf_writeSA','cmd':DpfSAWriter(),'gui':DpfSAWriterGUI},
    {'name':'ADdpf_writeLS','cmd':DpfLSWriter(),'gui':DpfLSWriterGUI},
    #{'name':'ADdpf_writeCluster','cmd':DpfClusterWriter(),'gui':DpfClusterWriterGUI},
    {'name':'ADdpf_edit','cmd':DpfEditor(),'gui':DpfEditorGUI},
]



def initModule(vf):

    for dict in commandList:
        vf.addCommand(dict['cmd'], dict['name'], dict['gui'])

    if vf.hasGui:
        for item in list(vf.GUI.menuBars['AutoToolsBar'].menubuttons.values()):
            item.configure(background = 'tan')
        if not hasattr(vf.GUI, 'adtBar'):
            vf.GUI.adtBar = vf.GUI.menuBars['AutoToolsBar']
            vf.GUI.adtFrame = list(vf.GUI.adtBar.menubuttons.values())[0].master




        
