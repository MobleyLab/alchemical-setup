#!/bin/env python
# Written by David L. Mobley and Pavel V. Klimovich at UC Irvine, November 2014.

#===================================================================================================
# IMPORTS
#===================================================================================================

import os
import numpy
from StringIO import StringIO
from optparse import OptionParser
from collections import defaultdict

#===================================================================================================
# INPUT OPTIONS
#===================================================================================================

parser = OptionParser()
parser.add_option('-m', dest = 'maptxt',  help = 'The atom map .txt file.',      default = 'input_files/map.txt')
parser.add_option('-a', dest = 'top_A',   help = 'The .top file of molecule A.', default = 'input_files/sampl4_1.top')
parser.add_option('-b', dest = 'top_B',   help = 'The .top file of molecule B.', default = 'input_files/sampl4_42.top')
parser.add_option('-s', dest = 'mcss',    help = 'The substructure .mol2 file.', default = '')
parser.add_option('-o', dest = 'out_dir', help = 'The output directory name.',   default = 'output_files')
parser.add_option('-v', dest = 'verbose', help = 'Verbosity.',                   default = 0, type=int)

#===================================================================================================
# The FEP_Molecule object to be built on the .top and .gro files.
#===================================================================================================

class FEP_Molecule:

   def __init__(self, mol_ID, top):

      def grep_b_Py(f):
         """ grep ^\\[ f -b """
         if not type(f) is file:
            with open(f, 'r') as f:
               return grep_b_Py(f)
         ordered_sn = [('uppermost', 0)]
         for line in iter(f.readline, ''):
            if line.startswith('['):
               ordered_sn += [(line.split()[1], f.tell()-len(line))]
         f.seek(0,2)
         ordered_sn += [('eof', f.tell())]
         return ordered_sn

      def dictionarize(l):
         """ [('uppermost', 0), (), ... ] ==> {'moleculetype': (583, 64), 'pairs': (2869, 974), 'angles': (3843, 1777), ... } """
         print 200*'lllllllllllllllll', l
	 global SN
         #sn, d = zip(*l)[0], {}
         sn, d = [], {}
         for i in range(len(l)-1):
            (k, v1), v2 =l[i], l[i+1][1]
            if k == 'dihedrals':
               k += '_%s' % '8'
            d[k] = (v1, v2-v1)
            sn.append(k)
	 if O.verbose:
            print "Molecule %s topology file's sections:\n%s\n%s" % (mol_ID, sn, d)
	 SN = sn[:-1]
	 print SN
	 return sn[:-1], d

      # Grep the directive titles from the topology file along with their byte offset.
      self.sn, self.d = dictionarize(grep_b_Py(top))

      self.ID               = mol_ID           # The molecule ID: either 'A' (the initial state) or 'B' (the final state).
      self.top              = top              # The topology file name.
      self.gro              = top[:-3] + 'gro' # The coorinates file name.
      self.mol2             = top[:-3] + 'mol2'# The .mol2 file name.
      self.old_new          = {}               # Initialize a dictionary to store the atom numbers mapping.
      self.type_charge_mass = {}

   def extractSection(self, section):
      """Extract the topology file section (directive)."""
      starting_position, up_to = self.d[section]
      with open(self.top, 'r') as in_file:
         in_file.seek(starting_position)
         section_string = in_file.read(up_to)
         if section.split('_')[0] not in ['atoms', 'atomtypes', 'pairs', 'bonds', 'angles', 'dihedrals']:
            if section=='uppermost':
               section_string = '#include "amber99sb.ff/forcefield.itp"\n' + section_string
            return section_string
         else:
	    if O.verbose:
               print "%s\n%s\n%s" % (('Molecule '+self.ID).center(36,'='), section_string, '='*36)
            print section_string
            return numpy.loadtxt(StringIO(section_string), dtype=str, skiprows=1, comments=';')

def initializeChimericFrom(mol):
   """Initialize molecule X that would contain the output data."""
   import copy
   molX = copy.copy(mol)              # Create a (shallow) copy of the mol object.
   molX.oo_bonded = defaultdict(list) # A dictionary to store the lines for the output [ bonds ], [ angles ], and [ dihedrals ] sections.
   molX.bonded_at = defaultdict(dict) # A dictionary to store the atom type sequences (will be of service in Part IV).
   molX.first = 'Integer'             # The index of the first atom in the chimeric molecule.
   molX.dummy = []                    # A list to keep track of the types of dummy atoms.
   molX.nrs_to_atom_types = {}        # A dictionary to store a pair (typeA, typeB) correponding to the new atom number.
   return molX

#===================================================================================================
# Part I. Figure out the matched atoms.
#===================================================================================================

def getRetainedIndices():

   def readMolecule(filename):
      """Read in a molecule from a .mol2 file as the OE object."""
      istream = oechem.oemolistream()
      istream.open(filename)
      molecule = oechem.OEMol()
      oechem.OEReadMolecule(istream, molecule)
      istream.close()
      return molecule

   def getMatchedAtoms(mol2_molecule, mol2_substructure, atomexpr, bondexpr, debug=False):
      """Map the mol2_molecule atoms onto those of mol2_substructure and return their indices as a list."""
 
      molecule = readMolecule(mol2_molecule)
      substructure = readMolecule(mol2_substructure)
 
      oechem.OEPerceiveChiral( molecule )
      oechem.OEPerceiveChiral( substructure )
      initial_substr_size = substructure.NumAtoms()
 
      mcss = oechem.OEMCSSearch( substructure, atomexpr, bondexpr, oechem.OEMCSType_Approximate )
      mcss.SetMinAtoms( initial_substr_size - 1 )
      mcss.SetMaxMatches( 5 )
      mcss.SetMCSFunc( oechem.OEMCSMaxAtomsCompleteCycles() )
 
      for match in mcss.Match(molecule, True):
         retained_indices = [matchpair.target.GetIdx() for matchpair in match.GetAtoms()]
 
         if match.NumAtoms() == initial_substr_size:
            return retained_indices
      parser.error("\nOE can not perform the mapping.")

   maptxt = O.maptxt
   if maptxt:
      try:
         molA.retained_indices, molB.retained_indices = numpy.genfromtxt(maptxt, dtype=int).transpose() - 1
      except:
         parser.error("\nCheck the format of the file %s.\nEach line of the file should represent a pair of the matched atom indices of the molecules 'A' (left) and 'B' (right)." % maptxt)
   else:
      substructure_file = O.mcss
      if not os.path.isfile(substructure_file):
         parser.error("\nFile '%s' not found and the OpenEye mapping search cannot be performed.\nTo get mapping from the .txt file rerun with the option '--maptxt'." % substructure_file)
      from openeye import oechem ## this is not a built-in module ##
      ae = oechem.OEExprOpts_EqONS
      be = (oechem.OEExprOpts_BondOrder | oechem.OEExprOpts_EqSingleDouble | oechem.OEExprOpts_EqAromatic)
      molA.retained_indices = getMatchedAtoms(molA.mol2, substructure_file, ae, be)
      molB.retained_indices = getMatchedAtoms(molB.mol2, substructure_file, ae, be)
   return

#===================================================================================================
# Part II: Build up the [ atoms ] section of the output topology file.
#===================================================================================================

def buildAtoms():

   def preprocessGro(f, residue='TMP', overlay=False):
      """Extract the coordinates of the 'residue' atoms and, if needed, modify them by overlaying onto the 'xyz_A'."""

      def grepPy(f, s):
         """ grep s f """
         if not type(f) is file:
            with open(f, 'r') as f:
               return grepPy(f, s)
         gl = ''
         for line in f:
            if s in line:
               gl += line
         return gl

      def optimalOverlayWithSVD(B_full, A_full, atoms_map={'A': [1, 0, 3, 8, 5], 'B': [2, 7, 1, 4, 3]}):
         """Return the coordinates of B_full optimally overlayed onto A_full (both are numpy arrays of shape Num_Atoms x 3) based on the atoms mapping."""
  
         # Define the substructures to be superimposed; move the origin to the geometrical center.
         B = B_full[atoms_map['B']]
         A = A_full[atoms_map['A']]
         centerB = B.mean(axis=0)
         centerA = A.mean(axis=0)
         B_full -= centerB
         B      -= centerB
         A      -= centerA
  
         # Perform the singular value decomposition of the cov(BA), find the optimal rotation matrix UV, and rototranslate B_full.
         U, S, V = numpy.linalg.svd(numpy.dot(numpy.transpose(B), A))
         B_full = numpy.dot(B_full, numpy.dot(U, V))
         return B_full+centerA

      xyz = grepPy(f, residue)
      if O.verbose:
         print "The excerpt of the residue %s from the .gro file %s:\n%s" % (residue, f, xyz)
      uc = (4,5,6) if len(xyz[:10].split())==2 else (3,4,5) # First two columns are sometimes glued up (as a result of "%5s%-5s"). Reading columns backwards would require a check for the presence of the velocities.
      xyz = numpy.genfromtxt(StringIO(xyz), dtype=numpy.float32, usecols=uc)
      if overlay:
         atoms_map = {'A':molA.retained_indices, 'B':molB.retained_indices}
         xyz = optimalOverlayWithSVD(xyz, xyz_A, atoms_map=atoms_map)
         if O.verbose:
            print "Coordinates after overlay:\n%s" % xyz
      return xyz

   def preprocessAtoms(mol, name_alpha_counter, hh, oo):
      """Read in the [ atoms ] section of the topology file of the FEP_Molecule 'mol' representing
      each data entry in a form of dictionary."""
      atoms = mol.extractSection('atoms')
      mol.ligand_acronym = atoms[0][3]
      nrs = atoms[:,0]
      molX.first = int(nrs[0])
      mol.retained_nrs = nrs[mol.retained_indices]
      xyz_A = preprocessGro(mol.gro, mol.ligand_acronym, overlay=(mol.ID!='A'))
      for (i, els) in enumerate(atoms):
         # Convert the data entry in a dictionary format.
         l = dict(zip(['nr', 'type', 'resnr', 'residue', 'atom', 'cgnr_old', 'charge', 'mass'], els))
         # Append atom's coordinates.
         l['xyz'] = "%8.3f%8.3f%8.3f" % tuple(xyz_A[i])
         # Memorize of which molecule the atoms are.
         l['mol_ID'] = mol.ID
         # Store the parameters of the properties that may potentially be perturbed.
         mol.type_charge_mass[l['nr']] = (l['type'], l['charge'], l['mass'])
         # To avoid duplicates, skip the common substructure atoms of the B molecule.
         if mol.ID == 'B':
            if l['nr'] in mol.retained_nrs: continue
         # Assign new atom name.
         name_alpha = filter(type(l['atom']).isalpha, l['atom'])
         name_alpha_counter[name_alpha] += 1
         l['name_new'] = '%s%d' % (name_alpha, name_alpha_counter[name_alpha])
         # This will help to sort out the entries.
         l['alpha_new'] = name_alpha
         l['digit_new'] = name_alpha_counter[name_alpha]
         # Store the hydrogens and heavy atoms separately.
         if name_alpha == 'H':
            hh.append(l)
         else:
            oo.append(l)
      if mol.ID=='A':
         return xyz_A, name_alpha_counter, hh, oo
      # Sort the atoms by new name.
      molX.oo_atoms = sorted(oo, key=lambda i: (i['alpha_new'], i['digit_new'])) + hh
      return

   def postprocessAtoms():
      dA_dB = dict(zip(molA.retained_nrs, molB.retained_nrs)) # Map the common substructure atom numbers.
      ind_new = molX.first-1
      molX.newtop_atoms = []
      for l in molX.oo_atoms:
         ind_old = l['nr']
         ind_new += 1
         l['ind_new'] = ind_new
         l['cgnr'] = ind_new
         if l['mol_ID'] == 'A':
            if ind_old in molA.retained_nrs:
               molB.old_new[ dA_dB[ind_old] ] = ind_new
               l['comments'] = '; MCSS'
               l['typeB'], l['chargeB'], l['massB'] = molB.type_charge_mass[ dA_dB[ind_old] ]
               if l['typeB'] != l['type']:
                  l['comments'] += '; to be perturbed'
            else:
               l['comments'] = '; to be annihilated'
               l['typeB'], l['chargeB'], l['massB'] = l['type']+'_dummy', '0.00000', l['mass']
               molX.dummy.append(l['type'])
            molA.old_new[ ind_old ] = ind_new
         else:
            molB.old_new[ ind_old ] = ind_new
            l['typeB'], l['chargeB'], l['massB'] = l['type'], l['charge'], l['mass']
            l['type'] = l['typeB'] + '_dummy'
            l['charge'] = '0.00000'
            l['comments'] = '; to be appeared'
            molX.dummy.append(l['typeB'])
         molX.nrs_to_atom_types[ind_new] = (l['type'][:2], l['typeB'][:2])
         molX.newtop_atoms.append(getString['atoms'] % l)
      return

   xyz_A, name_alpha_counter, hh, oo = preprocessAtoms(molA, defaultdict(int), [], [])
   preprocessAtoms(molB, name_alpha_counter, hh, oo)
   postprocessAtoms()
   return

#===================================================================================================
# Part III: Build up the [ atomtypes ] and [ pairs ] section.
#===================================================================================================

def buildAtomTypes():
   
   def preprocessAtomTypes(mol):
      """Preprocess atomtypes."""
      oo_atomtypes = []
      oo_atomtypes2 = []
      for els in mol.extractSection('atomtypes'):
         l = dict(zip(['name', 'bond_type', 'mass', 'charge', 'ptype', 'sigma', 'epsilon'], els))
         oo_atomtypes.append(getString['atomtypes'] % l)
         if l['name'] in molX.dummy:
            molX.dummy.remove(l['name'])
            l['name'] += '_dummy'; l['epsilon'] = '0.00000e+00'
            oo_atomtypes2.append(getString['atomtypes'] % l)
      return oo_atomtypes, oo_atomtypes2

   def zigZag((l1, l2), (l3, l4)):
      return l1+l3+l2+l4

   molX.newtop_atomtypes = zigZag(preprocessAtomTypes(molA), preprocessAtomTypes(molB))
   return

def buildPairs():
   
   def preprocessPairs(mol):
      """Preprocess pairs."""
      oo_pairs = []
      map_old_new = mol.old_new
      mol_ID = mol.ID
      for i, j, funct in mol.extractSection('pairs'):
         i_new = map_old_new[i]
         j_new = map_old_new[j]
         ijf = tuple(sorted((i_new, j_new))) + (funct,)
 
         if mol_ID=='B':
            if all(ind in mol.retained_nrs for ind in [i, j]): continue
         oo_pairs.append(getString['pairs'] % ijf)
      return oo_pairs

   molX.newtop_pairs = sorted(preprocessPairs(molA) + preprocessPairs(molB))
   return

#===================================================================================================
# Part IV: Build up the [ bonds ], [ angles ], and [ dihedrals ] sections.
#===================================================================================================

def buildBonded():

   def preprocessBonded(mol, section):
      """Preprocess bonds, angles, dihedrals."""
      nn = {'bonds':2, 'angles':3, 'dihedrals':4}[section.split('_')[0]]
      map_old_new = mol.old_new
      map_type = mol.type_charge_mass
      mol_ID = mol.ID
      print "molX.nrs_to_atom_types\n", molX.nrs_to_atom_types
      print "maptype\n", map_type
      print "section", section
      for els in mol.extractSection(section):
         print mol_ID, els
         ijkl, funct, params_list = els[:nn], els[nn], els[nn+1:]
         #mult = params_list[-1][0] # Does not make sense.
         mult = 'X' #params_list[-1][0] # Does not make sense.
         parameters = ''
         for p in params_list:
            parameters += "%14s" % p
         ijkl_new, k, k2 = '', '-', '-'
         for ind in ijkl:
            ind_new = map_old_new[ind]
            k += '%s-' % map_type[ind][0] # Old type.
            k2 += '%s-' % molX.nrs_to_atom_types[ind_new][1] # New type.
            ijkl_new += "%6s" % ind_new
 
         k2 = mult + k2 + mult # New.
         k = mult + k + mult   # Old.
         if O.verbose:
            print (mol_ID, k, k2)
         if k not in molX.bonded_at:
            molX.bonded_at[section][k] = parameters
         if mol_ID=='B':
            if all(ind in mol.retained_nrs for ind in ijkl): continue
         molX.oo_bonded[section].append(dict(zip(['type_sequence', 'ijkl', 'funct', 'parameters', 'k2'], [k, ijkl_new, funct, parameters, k2])))
      for k,v in molX.bonded_at[section].items():
         molX.bonded_at[section]['-'.join(k.split('-')[::-1])] = v
      print "\n\n\nmolX.bonded_at[section]\n", molX.bonded_at[section]
      return

   def postprocessDihedrals(): 
      d_funct = {'1':0, '2':0, '3':0, '4':1, '5':0, '8':8, '9':1}
      for (i, l) in enumerate(molX.oo_bonded[bt]):
         #print l
	 try:
            l['parametersB'] = ' '*8 + molX.bonded_at[bt][ l['k2'] ]
	 except KeyError:
            l['parametersB'] = l['parameters']

         if d_funct[ l['funct'] ]:
            # The to-be-appeared set of parameters.
            if l['parameters'][-1] != parametersB[-1]:
               parametersB_left = parametersB.split()
               parametersB_left[1] = '0.00000'
               l['parametersB_left'] = "%14s%14s%14s" % tuple(parametersB_left)

               parameters_right = l['parameters'].split()
               parameters_right[1] = '0.00000'
               l['parameters_right'] = ' '*8 + "%14s%14s%14s" % tuple(parameters_right)
               molX.oo_bonded[bt][i] = "%(ijkl)s%(funct)8s%(parameters)s%(parameters_right)s\n%(ijkl)s%(funct)8s%(parametersB_left)s%(parametersB)s\n" % l
            else:
               if O.verbose:
                  print "Dih else"
               molX.oo_bonded[bt][i] = "%(ijkl)s%(funct)8s%(parameters)s%(parametersB)s\n" % l

         else:
            if O.verbose:
               print "    Dih else22"
            molX.oo_bonded[bt][i] = "%(ijkl)s%(funct)8s%(parameters)s%(parametersB)s\n" % l

      molX.offset_dih = (len(l['parameters']), len(l['parametersB']))
      return

   def postprocessAnglesBonds():
      for (i, l) in enumerate(molX.oo_bonded[bt]):
         #print l
	 try:
            l['parametersB'] = ' '*8 + molX.bonded_at[bt][ l['k2'] ]
	 except KeyError:
            l['parametersB'] = l['parameters']
         if O.verbose:
            print ('aa', l['type_sequence'], l['k2'])
            print "zz%(type_sequence)s %(ijkl)s %(funct)s %(parameters)s %(parametersB)s %(k2)s" % l
         molX.oo_bonded[bt][i] = "%(ijkl)s%(funct)8s%(parameters)s%(parametersB)s\n" % l
      return

   lll = [z for z in SN if (numpy.array(['bo', 'an', 'di']) == z[:2]).any()]
   #for bt in ['bonds', 'angles', 'dihedrals']:
   for bt in lll:
      preprocessBonded(molA, bt)
      preprocessBonded(molB, bt)
      print "=== %s ===" % bt
      print lll
      print SN
      if bt.split('_')[0] == 'dihedrals':
         postprocessDihedrals()
      else:
         postprocessAnglesBonds()
   return

#===================================================================================================
# Part V: Write out the .gro file.
#===================================================================================================

def prepareGroFile(fi, fo):
   """Write out the .gro file linewise."""
   # Line number counter.
   ind = 0
   # The molecule starts at this line number.
   beg_mol = molX.first
   # The number of new atom lines to be inserted.
   extra_atoms = len([l for l in molX.oo_atoms if l['mol_ID']=='B'])
   # Generate the file.
   with open(fi, 'r') as in_file:
      with open(fo, 'w') as out_file:
         # The first line of the file: generate a new title.
         out_file.write(in_file.next()[:-1]+' -copied from- '+fi+'\n')
         # The second line of the file: update the number of the atom entries.
         n_lines = int(in_file.next())
         out_file.write("%5d\n" % (n_lines + extra_atoms))
         # Modify the rest of the file.
         for line in in_file:
            ind += 1
            # The lines prior to the molecule, as well as the last line (the box vector) are to be written unchanged,
            if (ind < beg_mol or ind > n_lines+2):
               out_file.write(line)
            # while those following the molecule need their indices to be updated.
            elif ind > beg_mol:
               out_file.write("%s%5d%s" % (line[:15], ind, line[20:45]))
            # As for the molecule lines themselves, they are to be replaced with those prepared in the 'molX.oo_atoms'.
            else:
               # Eliminate the old molecule lines.
               len_new = len(molX.oo_atoms)-1
               ind += len_new
               for i in range(len_new-extra_atoms):
                  in_file.next()
               # Write out new molecule lines.
               for l in molX.oo_atoms:
                  line = "%(resnr)5s%(residue)-5s%(name_new)5s%(ind_new)5s%(xyz)s\n" % l
                  out_file.write(line)
   return
       
#===================================================================================================
# Part VI: Write out the .top file.
#===================================================================================================

def writeOutTop(mol, fo):
   """Write out the topology file of the chimeric molecule."""
   print "Writing out the topology file %s" % fo
   rename = {'ind_new':'nr', 'residue':'res', 'name_new':'atom'}

   def formatHeadString(s, sec_name, names):
      if O.verbose:
         print (sec_name, names)
      i1, i2 = s.find('('), s.find(')')
      if i1<0:
         if not names:
            import re
            fi, se, th = [int(i) for i in re.sub('[sd%\n]', ' ', s).split()]
            return "[ %s ]\n; %*s%*s%*s\n" % (sec_name, fi-2, 'ai', se, 'aj', th, 'funct')
         n, s = s.split('s', 1)
         n = "[ %s ]\n; %s" % (sec_name, '%'+str(int(n.strip('%s'))-2)+'s')
         return (n+s) % tuple(names)
      else:
         chunk = s[i1:i2+1]
         name = chunk[1:-1]
         if name in rename.keys():
            name = rename[name]
         names.append(name)
         s = s.replace(chunk, '')
         return formatHeadString(s, sec_name, names)

   with open(fo, 'w') as out_file:
      for sec_name in SN:
         print "Building up the [ %s ] section ..." % sec_name
         if sec_name in ['atoms', 'atomtypes', 'pairs']:
            out_file.write(formatHeadString(getString[sec_name], sec_name, []))
            out_file.writelines(getattr(mol, 'newtop_'+sec_name) + ["\n"])
         #elif sec_name in ['bonds', 'angles', 'dihedrals']:
         elif sec_name[:2]=='bo' or sec_name[:2]=='an' or sec_name[:2]=='di':
            out_file.write("[ %s ]\n" % sec_name.split('_')[0])
            offset_ndx = {'bonds':2, 'angles':3, 'dihedrals':4}[sec_name.split('_')[0]]*6
            offset_prm = molX.offset_dih if sec_name=='dihedrals' else (28,28+8)
            #out_file.write("; %*s %7s %*s %*s\n" % (offset_ndx-2, 'atomnrs', 'funct', offset_prm[0]-1, 'parametersA', offset_prm[1]-1, 'parametersB'))
            out_file.write("; %*s %7s %s %s\n" % (offset_ndx-2, 'atomnrs', 'funct', 'parametersA'.rjust(offset_prm[0]-1, '<'), 'parametersB'.rjust(offset_prm[1]-1, '<')))
            out_file.writelines(mol.oo_bonded[sec_name] + ["\n"])
         else:
            out_file.write(mol.extractSection(sec_name))
   return

#===================================================================================================
# MAIN
#===================================================================================================

if __name__ == "__main__":

   getString = {   'atoms': "%(ind_new)6s %(type)10s %(resnr)6s %(residue)6s %(name_new)6s %(cgnr)6s %(charge)10s %(mass)10s %(typeB)10s %(chargeB)10s %(massB)10s %(comments)s\n",
               'atomtypes': "%(name)17s %(bond_type)10s %(mass)12s %(charge)12s %(ptype)12s %(sigma)12s %(epsilon)12s\n",
                   'pairs': "%6d%6d%8s\n",
                  'bonded': "%(ijkl)s%(funct)8s%(parameters)s%(parametersB)s\n"} # TODO: 'bonded' (offset_dih).
  
   O = parser.parse_args()[0]
   outputdir = O.out_dir
   if not os.path.isdir(outputdir):
      os.mkdir(outputdir)
  
   # Initialize molecules A and B, as well as X which is to contain the output data.
   molA = FEP_Molecule('A', O.top_A)
   molB = FEP_Molecule('B', O.top_B)
   molX = initializeChimericFrom(molA)
  
   getRetainedIndices()                                         # Part   I.
   buildAtoms()                                                 # Part  II.
   buildAtomTypes()                                             # Part III.
   buildPairs()                                                 # Part III.
   buildBonded()                                                # Part  IV.
   prepareGroFile(molA.gro, os.path.join(outputdir, 'out.gro')) # Part   V.
   writeOutTop(molX, os.path.join(outputdir, 'out.top'))        # Part  VI.

#===================================================================================================
#                                   End of the script 
#===================================================================================================
