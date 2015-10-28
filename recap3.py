# recap3: from 
# Copyright 2009 TJ O'Donnell

import sys
import openbabel as ob


# each smarts must contain only two atoms representing the bond to be broken.
# of course, each atom may be a complex atom smarts, ala [$(whatever)]

BondBreakSmarts = {
' 1.amide':'[$([C;!$(C([#7])[#7])](=!@[O]))]!@[$([#7;+0;!D1])]',
' 2.ester':'[$(C=!@O)]!@[$([O;+0])]',
' 3.amine':'[$([N;!D1;+0;!$(N-C=[#7,#8,#15,#16])](-!@[*]))]-!@[$([*])]',
' 4.urea':'[$(C(=!@O)([#7;+0;D2,D3])!@[#7;+0;D2,D3])]!@[$([#7;+0;D2,D3])]',
' 5.ether':'[$([O;+0](-!@[#6!$(C=O)])-!@[#6!$(C=O)])]-!@[$([#6!$(C=O)])]',
' 6.olefin':'C=!@C',
' 7.quaternaryN':'[N;+1;D4]!@[#6]',
' 8.aromaticN-aliphaticC':'[$([n;+0])]-!@C',
' 9.lactamN-aromaticC':'[$([O]=[C]-@[N;+0])]-!@[$([C])]',
'10.aromaticC-aromaticC':'c-!@c',
'11.sulphonamide':'[$([#7;+0;D2,D3])]-!@[$([S](=[O])=[O])]'
}

class Recap:
  # RECAP-Retrosynthetic Combinatorial Analysis Procedure
  # J. Chem. Inf. Comput. Sci. 1998, 38, 511-522
  def __init__(self, mol, minsize=5,newBondSmarts=None):
	self.mol = mol;
	# minimum allowed size (atom count) of fragment
	self.minsize = minsize;
	# bonded atom pairs populated by the apply method,
	# subsequently used by split and add_star
	self.atom_pairs = list()

	if newBondSmarts != None:
		self.smarts = newBondSmarts
		 
	else:
		self.smarts = BondBreakSmarts
 
	# NB: order is important to Recap.apply()!
	self.bondNames = self.smarts.keys()
	self.bondNames.sort() 

  def apply(self, pat, patnum):
	if pat.Match(self.mol):
	  # find all atom pairs that match
	  for p in pat.GetUMapList():
		i = 0
		atoms = list()
		for a in ob.OBMolAtomIter(self.mol):
		  i += 1
		  if i in p:
			atoms.append(a)
		if self.small_fragment(atoms[0], atoms[1]):
		  #print True
		  pass
		else:
		  atoms.append(patnum)
		  self.atom_pairs.append(atoms)
		  #print False
	  return True
	else:
	  return False

  def split(self, label=None):
	for a in self.atom_pairs:
	  if label:
		a[0].SetIsotope(a[2])
		a[1].SetIsotope(a[2])
	  bond = a[0].GetBond(a[1])
	  # bond could be null if already deleted when smarts matched multiple times
	  if bond: self.mol.DeleteBond(bond)

  def add_star(self):
	for pair in self.atom_pairs:
	  self.mol.AddBond(pair[0].GetIdx(),self.mol.NewAtom().GetIdx(),1)
	  self.mol.AddBond(pair[1].GetIdx(),self.mol.NewAtom().GetIdx(),1)

  def decide_multiples(self):
	# some smarts (e.g. ether, amine) allow multiple bonds to the
	#  central atom to be broken.  Yet it appears the central atom
	#  needs to be retained in one of the multiple fragments.
	#  If multiple fragments, let it stay with the smallest fragment.
	#  If tied, pick the first fragment.
	multiples = [0]*self.mol.NumAtoms()
	for pair in self.atom_pairs:
	  multiples[pair[0].GetIdx()] += 1
	  multiples[pair[1].GetIdx()] += 1
	#print multiples

	currsize = -1
	currpair = None
	for pair in self.atom_pairs:
	  a = pair[0]
	  b = pair[1]
	  if multiples[a.GetIdx()] > 1 or multiples[b.GetIdx()] > 1:
		# remove larger fragment(s) if a-b were broken
		#print a.GetIdx(),b.GetIdx(),
		fsize = self.fragment_size(a,b)
		if currpair == None:
		  currpair = pair
		  currsize = fsize
		else:
		  if fsize < currsize:
			self.atom_pairs.remove(pair)
		  else:
			self.atom_pairs.remove(currpair)
			currpair = pair
			currsize = fsize

  def fragment_size(self, a, b):
	# size of fragment b if a-b were broken
	c1 = ob.vectorInt()
	self.mol.FindChildren(c1,a.GetIdx(),b.GetIdx())
	#for atom in c1:
	#  if self.mol.GetAtom(atom).GetValence() == 1:
	return 1+len(c1)
 
  def small_fragment(self, a, b):
	# if we were to break the bond between a and b,
	#  would either fragment be too small?
	#print a.GetIdx(), b.GetIdx(),
	if self.fragment_size(a,b) < self.minsize: return True
	if self.fragment_size(b,a) < self.minsize: return True

	return False

if __name__ == '__main__':

	import glob 
	
	pat = ob.OBSmartsPattern()
	obcM2C = ob.OBConversion()
	obcM2C.SetOutFormat("can")
	obcM2C.SetOptions("-i", obcM2C.OUTOPTIONS)
	obcP2M = ob.OBConversion()
	obcP2M.SetInAndOutFormats('pdbqt','mol')

	# hancock_VB
	# dockDir = '/media/sf_sharedData/coevol-HIV/WCG/anal/recapTest/'
	# monk
	dockDir = '/media/rik/e856a0e4-e02a-4d58-b68a-8f1aac37f2c1/coevol-HIV/WCG/anal/recapTest/'
	
	testFiles = glob.glob(dockDir+'*.pdbqt')
	print 'recap3 test: %d pdbqt found' % (len(testFiles))

	for pdbqf in testFiles:		
		print pdbqf
		ligMol = ob.OBMol()
		obcP2M.ReadFile(ligMol,pdbqf)
	
		cansmile0 = obcM2C.WriteString(ligMol,1)
		# this returns both the canonSmiles string, but also PDBQT file name?!
		# str: O=C1/C(=C/c2cccc(c2)N(=O)=O)/C[C@]2(C/C/1=C\c1cccc(c1)N(=O)=O)C(=O)Nc1c2cccc1    /media/sf_sharedData/coevol-HIV/WCG/anal/PrAS_115-120_150122/Dock_lowE_W_R/116/0406688_fahv.x3KFN_prASw0c0_ZINC40146065_1288553872_out_Vina_VS.pdbqt
    
		lsbits = cansmile0.split()
		canon2 = lsbits[0]
		
		print 'Canon:',canon2
	
		currRecap = Recap(ligMol,4)
		for si,bondName in enumerate(currRecap.bondNames):
			pat.Init(currRecap.smarts[bondName])
			currRecap.apply(pat, si)
		currRecap.decide_multiples()
		currRecap.split()
		# Recap.add_star()
		recapStr = obcM2C.WriteString(ligMol,1)
		lsbits = recapStr.split()
		recapStr2 = lsbits[0]
		print 'RECAP:',recapStr2
