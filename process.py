# process: process docking files for RACCOON2 interactions
#		 ALSO CRAWL before tgz-ing processed files
# v0.2
# 160108
# rbelew@ucsd.edu
#
#  Using Raccoon2	 
#  Copyright 2013, Stefano Forli
#  Molecular Graphics Lab
#  v.0.4 (modified version to include vdw interactions)
# 
#############################################################################

import datetime
import time
import os, sys
from glob import glob  
import re
from sys import argv, exc_info
import getopt
import shutil
import socket
import tempfile
import tarfile

import config
import crawl_ADV

from numpy import array, zeros, sqrt  
from bhtree import * # bhtreelib

from mglutil.math.rmsd import RMSDCalculator

from CADD.Raccoon2.piStackingAndRingDetection import *  #(pi-stacking interatcion)

from CADD.Raccoon2.HelperFunctionsN3P import getCoords, pmvAtomStrip, getLines, quickdist, dist, percent, getAtype, atomCoord # TODO
from CADD.Raccoon2.HelperFunctionsN3P import pathToList, writeList
from CADD.Raccoon2.WaterProcessing import processHydroDocking

## Utilities

FloatRegEx = r'[-+]?\d*\.\d+|\d+'

class VSDockingResultsGenerator:
	"""VS result class for extracting the results (LE,LC) from a series of DLG files
	   and calculate properties (L.eff, cluster statistics, receptor contacts)
	   - flexible residue interactions are calculated
	   - simple planar ring detection on ligand (and target co-factors) is implemented for pi-stackings
	   - Pi-stackings (parallel, t-stack) interactions are detected

	   if auto is True, the list of dlgs will be used immediately to run recluster results,
	   otherwise methods recluster() and calcInteractions() must be used.

	   Receptor structure can be specified as filename or as {text|array|atype} dictionary.

	   if a filename is specified, the file is processed by self.getCoords() and the self.recname is set.

	   if a pre-processed receptor dictionary is provided, self.recname must be specified explicitly.

	   receptor = { 'text' : atoms, 'coord' : array( coord, 'f'), 'atype': atype }

	   A simple usage could be this:

		result = VSDockingResultsGenerator(dlgList = dlg_list, rmsTol = 2.0, receptor = receptor_file, recname = None, auto = True)
		print result
		pdbqt = result.generatePDBQTplus()

	# print "ITERTOOLZ!" <- improve speed?
	# http://docs.python.org/library/itertools.html
	# http://docs.python.org/release/2.5.2/tut/node7.html
	# http://www.doughellmann.com/PyMOTW/itertools/

	"""
	#
	# TODO Note on the self.receptor
	# initialize the object with self.recPiVectors:
	# i.e. result = VSDockingResultsGenerator(dlgList = list)
	#	  result.recPiVectors(vectorsCached)
	# to save time from re-calculating every time the PiVectors
	# TODO ...or not necessary? Pi vectors should be attached to the receptor parsed
	#		  coords?  {text : ...., pi_vectors : [ ] }

	def __init__(self, input_files = None,		   # list of DLG[AD]/*_out.pdbqt[VINA] contaning poses
					   #rmsTol = None,		  # clustering RMSD tolerance
					   recname = None,		 # receptor name (if not provided it will be taken from 'self.receptor'
					   receptor = None,		# receptor coordinates: PDBQT or getCoord() output (CACHE)
					   doInteractions = True,  # calculate ligand-target interactions
					   debug = False,		  # debug mode
					   hbtol = 0.0,			# hb distance extra tolerance
					   water_map = None,	   # water map used for rescoring W atoms
					   fullentropy=0,		  # include entropy penalty for conserved waters # TODO REMOVE THIS?
					   ):
		
		# default values
		self.DEBUG = debug

		if self.DEBUG:
			print "VSDockingResultsGenerator> init:", input_files, recname , receptor, doInteractions, debug, hbtol, water_map
		# SF-DBG 
#		 if not hbtol == 0 :
#			 print "DEBUG: extra hbtolerance:", hbtol

		#print "VSDockingResultsGenerator> init:", input_files, recname , receptor, doInteractions, debug, hbtol, water_map
		
		HB_CUTOFF=3.21 + hbtol # X-Y hydrogen-bond distance cutoff -X-H...Y- # from AutoGrid 4.2 values
		VDW_SMOOTHING = .5 # vdw tolerance (account for AutoGrid map smoothing)
		# TODO TODO TODO TODO TODO
		# TODO these values must be adapted *differently* for VINA AND AUTODOCK? TODO TODO TODO 

		# /default values

		self.VDW_TOL = VDW_SMOOTHING
		self.HB_CUTOFF = HB_CUTOFF**2 # save time

		# from AD4_parameter.dat Rii/2 values
		self.vdw_radii = { 'H': 1.00, 'HD': 1.00, 'HS': 1.00, 'C': 2.00,
			'A': 2.00, 'N': 1.75, 'NA': 1.75, 'NS': 1.75, 'OA': 1.60,
			'OS': 1.60, 'F': 1.54, 'Mg': 0.65, 'MG': 0.65, 'P': 2.10,
			'SA': 2.00, 'S': 2.00, 'Cl': 2.04, 'CL': 2.04, 'Ca': 0.99,
			'CA': 0.99, 'Mn': 0.65, 'MN': 0.65, 'Fe': 0.65, 'FE': 0.65,
			'Zn': 0.74, 'ZN': 0.74, 'Br': 2.165, 'BR':2.165, 'I':2.36,
			'Z' : 2.00, 'G' : 2.00, 'GA': 2.00, 'J' :2.00, 'Q' :2.00,
			'W' : 3.0,  # TODO evaluate this value...
			'X': 2 } # default vdW for unknown atom

		# list of atoms that are ignored when calculating RMSD
		self.ignore_at = [ 'HD', 'H', 'W' ]
		#self.ignore_at = [ 'HD', 'H' ]
		#print " \n\n***\n\n Change the management of W atom types! contacts and so on...\notherwise they will not be correctly removed\n\n***\n\n"
		self.water_types = [ 'W' ]
		self.water_map = water_map
		self.isHydrated = False
		self.fullentropy = fullentropy
		self.recname = recname
		self.receptor = receptor
		self.doInteractions = doInteractions
		self.hbtol = hbtol

		# HB?
		#		self.hb_radii = {
		#'N'	  3.50  0.160  22.4493  -0.00162  0.0  0.0  0  -1  -1  1	# Non H-bonding Nitrogen
		#'NA'	 3.50  0.160  22.4493  -0.00162  1.9  5.0  4  -1  -1  1	# Acceptor 1 H-bond Nitrogen
		#'NS'	 3.50  0.160  22.4493  -0.00162  1.9  5.0  3  -1  -1  1	# Acceptor S Spherical Nitrogen
		#'OA'	 3.20  0.200  17.1573  -0.00251  1.9  5.0  5  -1  -1  2	# Acceptor 2 H-bonds Oxygen
		#'OS'	 3.20  0.200  17.1573  -0.00251  1.9  5.0  3  -1  -1  2	# Acceptor S Spherical Oxygen
		#'SA'	 4.00  0.200  33.5103  -0.00214  2.5  1.0  5  -1  -1  6	# Acceptor 2 H-bonds Sulphur
		#'S'	  4.00  0.200  33.5103  -0.00214  0.0  0.0  0  -1  -1  6	# Non H-bonding Sulphur
		#}

		self.metals = [ 'Mg', 'MG', 'Mn', 'MN', 'Fe', 'FE', 'Zn', 'ZN', 'CA', 'Ca'] # metals
		self.met_coord = [ 'OA', 'OS', 'SA', 'S', 'N', 'NA', 'NS' ]				 # metal-binders

		if input_files:
			self.setLigands(input_files = input_files)

		if receptor:
			self.setReceptor(receptor=receptor, recname=recname)

		
	def setLigands(self, input_files=[]): 
		""" set the input files to be processed
		"""
		self.input_files = input_files
		self.poses = []
		self.problematic = []
		self.ligName = None
		self.atomTypes = []
		self.results = []
		self.flexres = False
		self.interactions = False
		self.histogram = []
			

	def setReceptor(self, receptor, recname=None):
		# pre-process receptor if filename
		if isinstance(receptor, str):
			if recname == None:
				r = os.path.basename(receptor)
				recname = os.path.splitext(r)[0]
			try:
				receptor = getCoords(getLines(receptor))
			except:
				print "Problem in reading the receptor [%s]: (%s) " % (receptor, sys.exc_info()[1])
				return False

		self.receptor = receptor
		self.recname = recname

		# initialize the bh_tree for the 
		# interactions calculated later
		try:
			freeBHtree(self.rec_bht)
		except:
			pass
		self.rec_bht		   = bhtreelib.BHtree( self.receptor['coord'], None, 10)
		self.rec_bht_indices   = zeros( (len(self.receptor['coord']),) ).astype('i')
		self.rec_bht_distances = zeros( (len(self.receptor['coord']),) ).astype('f') 


	def getPath(self):
		""" return the path where the input data is stored
			by reading it from the first ligand file
		"""
		if not self.input_files:
			return
		if isinstance(self.input_files,list):
			f = self.input_files[0]
		elif isinstance(self.input_files, str):
			f = self.input_files
		path = os.path.dirname(f)
		return path

	
	def findHbAccepDon(self, atom_list):
		"""identifies HB donors and acceptors in a list of PDBQT atoms
		   returns : acceptors[] and donors[] lists
		"""
		# sf-DBG print "\n\nCALLEX"
		H_COV_BOND = 1.19 # value adapted for some distorted structures (ZINC21002974)
		H_COV_BOND  = H_COV_BOND ** 2  
		acceptor_types = ['OA', 'NA', 'SA']
		donor_types = ['N', 'O', 'OA', 'NA', 'SA']
		acceptors = []
		donors = []
		h = []
		dcandidate = []
		for l in atom_list:
			# sf-DBG print "LATOM", l
			if l.startswith("ATOM") or l.startswith("HETATM"):
				l = l.strip()
				atype=l.split()[-1]
				if atype in acceptor_types:
					if not l in acceptors:
						acceptors.append(l)
				if atype in donor_types:
					if not l in dcandidate:
						dcandidate.append(l)
				if atype == 'HD':
					if not l in h:
						h.append(l)
		for a in dcandidate:
			for x in h:
				if dist(a, x) <= H_COV_BOND:
					donors.append(a)
					break

		return acceptors, donors 

	def calcInteractions(self):
		if not self.results:
			if self.DEBUG: print "Warning: no results."
			return False
		self.getContactAtoms()
		self.getHbInteractions()

		# XXX think about the optimal place for water processing
		# XXX
		# water processing happens after 
		# close contacts and HB acc/don are known...
		if self.isHydrated: 
			if self.water_map:
				self.processWaters(poses = self.results)
			else:
				print " # WARNING: inaccurate energy! [%s] hydrated ligand found, but no water map specified" %  self.ligName

		self.getPiStackInteraction()
		self.getMetalCoordInteraction()
		self.interactions = True

	# Phase 2 functions
	def getContactAtoms(self):
		""" populate self.result[]['vdw_contacts'] list 
			saves receptor metals for potential usage later (metal-coordination)
		"""
		cutoff	= 5. # initial cutoff for bhtree
		#print len(self.receptor['coord'])


		# XXX THIS PART COULD BE MOVED IN THE INIT/set receptor part
		#bht	   = bhtreelib.BHtree( self.receptor['coord'], None, 10)
		#indices   = zeros( (len(self.receptor['coord']),) ).astype('i')
		#distances = zeros( (len(self.receptor['coord']),) ).astype('f')
		for r in xrange(len(self.results)):
			self.results[r]['metal'] = []
			# A. find ligand atoms close c.ontacts:
			# sf-DBGDBG print "\nPROCESSING", self.results[r]
			self.xxx = self.results
			for i in xrange( len(self.results[r]["coord"])):
				try:
					l_vdw = self.vdw_radii[ self.atomTypes[i] ] + self.VDW_TOL
				except:
					l_vdw = self.vdw_radii['X'] + self.VDW_TOL
					print "Warning! Unrecognized ligand atom type!", self.atomTypes[i]
				#print self.atomTypes[i], l_vdw,

				l_coord = self.results[r]["coord"][i]
				nb = self.rec_bht.closePointsDist( tuple(l_coord), 
						cutoff, self.rec_bht_indices, self.rec_bht_distances)
				# 1. rigid receptor atoms
				for j in xrange(nb):
					rec_index = self.rec_bht_indices[j]
					rec_atom = self.receptor['text'][rec_index]
					if not rec_atom in self.results[r]['vdw_contacts']: # slight speed up?
						d = self.rec_bht_distances[j]
						try:
							r_vdw = self.vdw_radii[self.receptor['atype'][rec_index]]
						except:
							r_vdw = self.vdw_radii['X']
							print "Warning! Unrecognized receptor atom type!"
						if d <= ( r_vdw + l_vdw ): 
							# vdW contact
							self.results[r]['vdw_contacts'].append((self.results[r]['true_ligand'][i], rec_atom))
							# metal detector
							if self.receptor['atype'][rec_index] in self.metals:
								self.results[r]['metal'].append(rec_atom)	  
				# sf-DBG print "VVVV", self.results[r]['vdw_contacts'] 
				self.results[r]['vdw_contacts'] = list(set(  self.results[r]['vdw_contacts'] ))
				del nb

				# 2. flexible residues atoms
				if self.flexres: 
					if self.DEBUG: print "get_ContactAtoms> calc_distances for FlexRes too!"
					for fres_idx in range(len(self.results[r]['flex_res'])):
						#print "FRESINDDX", fres_idx, self.results[r]['flex_res']['text'][fres_idx]
						r_coord = self.results[r]['flex_res']['coord'][fres_idx]
						r_vdw = self.vdw_radii[self.results[r]['flex_res']['atype'][fres_idx]]
						r_atom = self.results[r]['flex_res']['text'][fres_idx]
						#print "RATOM IS", fres_idx, r_atom
						if not r_atom in self.results[r]['vdw_contacts']:
							if quickdist(l_coord, r_coord) <= ( r_vdw + l_vdw ):
								# vdW contact (flex res)
								# XXX MODIFY FOR VDW PAIRS
								self.results[r]['vdw_contacts'].append(r_atom) 
								# metal detector (flex res; it shouldn't happen...?)
								if self.receptor['atype'][rec_index] in self.metals:
									self.results[r]['metal'].append(rec_atom)  
								if self.DEBUG:
									print "contact with flex_res atom:", r_atom
									print quickdist(l_coord, r_coord), r_vdw , l_vdw , r_vdw+l_vdw
							#else:
							#	print "MIOSSING"
			# B. find water close contacts
			if self.isHydrated:
				if self.DEBUG: print "getContactAtoms> calc.distances for W atoms too"
				self.results[r]['water_bridge_contacts'] = []
				for w in self.results[r]['water_bridge']:
					contacts = []

					w_coord = atomCoord(w)
					# w_coord = hf.atomCoord(w) # sf_150818 ??

					w_vdw = 2.0 + self.VDW_TOL # XXX water radius
					nb = self.rec_bht.closePointsDist(tuple(w_coord), cutoff, 
							self.rec_bht_indices, self.rec_bht_distances)
					for j in xrange(nb):
						rec_index = self.rec_bht_indices[j]
						rec_atom = self.receptor['text'][rec_index]
						if not rec_atom in contacts:
							d = self.rec_bht_distances[j]
							try:
								r_vdw = self.vdw_radii[self.receptor['atype'][rec_index]]
							except:
								r_vdw = self.vdw_radii['X']
								print "Warning! Unrecognized receptor atom type!"
							if d <= ( r_vdw + w_vdw ): 
								contacts.append(rec_atom)   
					del nb
					self.results[r]['water_bridge_contacts'].append(contacts)

			self.results[r]['rec_hb_candidates'] = {}
			acc, don = self.findHbAccepDon([x[1] for x in self.results[r]['vdw_contacts']])
			if self.DEBUG:
# 				writeList("rec_acceptors.pdb", acc, addNewLine=True)
# 				writeList("rec_donors.pdb", don, addNewLine=True)
				print "============================"
				print "ACC", acc
				print "============================"
				print "DON", don
				print "============================"
			self.results[r]['rec_hb_candidates']['acc'] = acc
			self.results[r]['rec_hb_candidates']['don'] = don
		# XXX this part should move to the initialization and be called every time 
		# a new receptor is called...
		# freeBHtree(bht)


	def getHbInteractions(self):
		"""
		input   : pdb lines (ligand, receptor)
		output  : list of acceptor/donor pairs (lig:rec) <= PDBQT+ format
		notes   : distance only is used to characterize HB (consistant with AG maps)
				  flex res atoms are handled in a transparent manner being added to:
					 self.results[r]['rec_hb_candidates']['acc'] = acc
					 self.results[r]['rec_hb_candidates']['don'] = don.
		"""
		for i in range(len(self.results)):
			# use the close contact atoms to calculate the HB (faster)
			hb_acc = []
			hb_don = []
			# HB acceptors
			for a in self.results[i]['hba_atoms']: # ligand hb acceptor atoms
				ltype = a.rsplit()[-1] 
				stol=0.
				if ltype=='SA': stol += 3.5 # TODO CHECK THIS VALUE!
					#print "A>",a, "[%2.2f]" % stol
				for r in self.results[i]['rec_hb_candidates']['don']:
					#rtype = r.rsplit()[-1] # TODO check this value
					rtype = getAtype(r)
					if rtype=='SA': stol += 3.5
					if dist(a,r,sq=False)-stol <= self.HB_CUTOFF:
						#print  dist(a,r)-stol, self.HB_CUTOFF
						#print  sqrt( dist(a,r)-stol), sqrt(self.HB_CUTOFF)
						if not [a,r] in hb_acc:
							hb_acc.append([a,r])
			# HB donors
			for a in self.results[i]['hbd_atoms']: # ligand hb donor atoms
				if a.rsplit()[-1]=='SA': stol=2
				else: stol=0
				for r in self.results[i]['rec_hb_candidates']['acc']:
					if dist(a,r)-stol <= self.HB_CUTOFF:
						if not [a,r] in hb_don:
							hb_don.append([a,r])
			self.results[i]['hb'] = { 'acceptors' : hb_acc, 'donors': hb_don}
			# TODO WATERS GO HERE?
			# Water bridges
			if self.isHydrated:
				pass


	def getPiStackInteraction(self):
		# TODO process the receptor in advance, cache the vectors!
		receptor = self.receptor['text']
		for i in range(len(self.results)):
			if self.flexres:
				#pstack, tstack = findLigRecPiStack(self.results[i]['text'], receptor+self.results[i]['flex_res']['text'])
				# TODO fix the flex_res bug! This godedoesn'twork
				true_lig_atoms = list(set(self.results[i]['text']) - set(self.results[i]['flex_res']['text']))
				# DEBUG XXX
				#writeList(str(i)+"_flex.pdb", self.results[i]['flex_res']['text'], addNewLine=True)
				#writeList(str(i)+"_lig.pdb", true_lig_atoms)
				pstack, tstack = findLigRecPiStack(true_lig_atoms[i], receptor+self.results[i]['flex_res']['text'])
			else:
				pstack, tstack = findLigRecPiStack(self.results[i]['text'], receptor)
			self.results[i]['pi'] = { 'p-stack' : pstack , 't-stack': tstack}

	def getMetalCoordInteraction(self):
		# it could be incorporated with getHbInteractions leading
		# to getPolarInteractions, but in this way it's easier to 
		# maintain... probably.
		metal_coordination = []
		for i in range(len(self.results)):
			for a in self.results[i]['hba_atoms']: # ligand hb acceptor atoms
				for r in self.results[i]['metal']: # receptor metals
					if dist(a,r) <= self.HB_CUTOFF:
						if not [a,r] in metal_coordination:
							metal_coordination.append([a,r])
			self.results[i]['metal_coord'] = metal_coordination

	def calcReferenceRms(self, reference, keepH=False, debug=False):
		ref = getCoords(getLines(reference), include_hydrogens = keepH)
		if not len(ref['coord']) == len(self.results[0]['coord']):
			print "[atoms mismatch] Warning! The reference stricture doesn't match the docked ligand!"
			print "Reference [ %d ] | Docked ligand [ %d ]" % (  len(ref['coord']), len(self.results[0]['coord']) )
			return 
		rmsdcalc = RMSDCalculator(ref['coord'])
		rmsd = [ rmsdcalc.computeRMSD(self.results[0]['coord']) ]
		if len(self.results) > 1:
			rmsd.append(rmsdcalc.computeRMSD(self.results[1]['coord']))
		if debug:
			dist_pool = []
			for px in self.poses:
				p = px['coord']
				d = []
				for i in range(len(p)):
					x = quickdist(p[i], ref['coord'][i], sq=True)
					d.append(x)
					d.sort()
				dist_pool.append(d)
				#print "================="
				#print d
			return rmsd, dist_pool,self.poses
		else:
			return rmsd

	def getInteractSummary(self, pose=0, recOnly=False):
		""" collect and format the interaction information for the pose
			if recOnly, the summary will inlcude only one entry
			for receptor residues involved in interactions
		"""
		interactions = {}
		pose = self.results[pose]
		if recOnly:
			if pose['hb']['acceptors']:
				interactions['hba'] = [ pmvAtomStrip(p[0]) for p in pose['hb']['acceptors'] ]
			if pose['hb']['donors']:
				interactions['hbd'] = [ pmvAtomStrip(p[0]) for p in pose['hb']['donors'] ]
			if pose['vdw_contacts']:
				interactions['vdw'] = [ pmvAtomStrip(p) for p in pose['vdw_contacts'] ]
			if pose['metal_coord']:
				interactions['metal'] = [ pmvAtomStrip(p[0]) for p in pose['metal_coord'] ]
			if pose['pi']['p-stack']:
				interactions['ppi'] = [ p[0] for p in pose['pi']['p-stack'] ]
			if pose['pi']['t-stack']:
				interactions['tpi'] = [ p[0] for p in pose['pi']['t-stack'] ]
		else:
			if pose['hb']['acceptors']:
				interactions['hba'] = [ (pmvAtomStrip(p[0]), pmvAtomStrip(p[1])) for p in pose['hb']['acceptors'] ]
			if pose['hb']['donors']:
				interactions['hbd'] = [ (pmvAtomStrip(p[0]), pmvAtomStrip(p[1])) for p in pose['hb']['donors'] ]
			if pose['vdw_contacts']:
				interactions['vdw'] = [ pmvAtomStrip(p) for p in pose['vdw_contacts'] ]
			if pose['metal_coord']:
				interactions['metal'] = [ (pmvAtomStrip(p[0]), pmvAtomStrip(p[1])) for p in pose['metal_coord'] ]
			if pose['pi']['p-stack']:
				interactions['ppi'] = [ pose['pi']['p-stack'] ]
			if pose['pi']['t-stack']:
				interactions['tpi'] = [ pose['pi']['t-stack'] ]
		return interactions



################################################## 
################################################## AutoDock class 
################################################## 

class AutoDockVsResult(VSDockingResultsGenerator):

	def __init__(self, input_files = None, rmsTol = None, 
				recname = None, receptor = None, auto = True, 
				doInteractions = True, hbtol = 0.0, 
				ignoreWaters = False, # XXX
				water_map = None,
				fullentropy=0,
				debug = False):

		VSDockingResultsGenerator.__init__(self,
			input_files = input_files,		  # list of DLG[AD] from wich poses will be extracted
			recname = recname,				  # receptor name (if not provided it will be taken from 'self.receptor'
			receptor = receptor,				# receptor coordinates: PDBQT or getCoord() output (CACHE)
			doInteractions = doInteractions,	# calculate ligand-target interactions
			debug = debug,					  # debug mode
			hbtol = hbtol,
			water_map = water_map,
			fullentropy=fullentropy,
			)
		# defaults
		DEFAULT_RMSTOL = 2.0
		
		# AutoDoc specific properties
		self.tag = ["LE", "LC"]
		self.clusters = None
		self.histogram = []
		self.rmsTol = DEFAULT_RMSTOL # angstro		self.ignoreWaters = ignoreWaters
		self.isHydrated = False
		self.doInteractions = doInteractions
		if rmsTol:
			self.rmsTol = rmsTol
		if water_map:
			self.water_map = water_map
		#
		if auto and input_files:
			self.process()

	def __str__(self):
		text = "\n===================================================\n"
		text +=" This is my VSDockingResults object. [AutoDock]\n"
		text +=" There are many like it, but only this one\n"
		text +=" has been made from => %s\n" % self.ligName
		text +="  - receptor name			: %s\n" % self.recname
		text +="  - total runs			   : %s\n" % self.totRuns
		text +="  - results poses			: %s\n" % len(self.results)
		text +="  - total clusters		   : %s\n" % len(self.histogram)
		text +="  - flex res				 : %s" % self.flexres
		text +="\n Result(s) ------------------------------------------\n"
		if self.interactions:
			text += "	   Energy	L.eff  Csize	Cpc   vdW	  HB	Pi-Pi\n"
			for i in range(len(self.results)):
				text += " %s   %2.3f\t%2.3f\t %d\t%3.2f%%\t%d\t%d\t%d\n" % \
					( self.tag[i], self.results[i]['energy'], self.results[i]['leff'],
					  self.results[i]['csize'], self.results[i]['cpercent'],\
					  len( self.results[i]['vdw_contacts']),
					  (len(self.results[i]['hb']['acceptors'])+len(self.results[i]['hb']['donors'])),
					  len(self.results[i]['pi']['p-stack'])+len(self.results[i]['pi']['t-stack'])
					)
		else:
			text += "	   Energy	L.eff  Csize	Cpc\n"
			for i in range(len(self.results)):
				text += " %s   %2.3f\t%2.3f\t %d\t%3.2f%%\n" % \
					( self.tag[i], self.results[i]['energy'], self.results[i]['leff'],
					self.results[i]['csize'], self.results[i]['cpercent'])

		text +=" Histogram ------------------------------------------\n"
		c = 1
		for i in self.histogram:
			text+= "   %d.  %2.3f | %d\t%s" % (c, i[0],i[1] ,"#"*i[1])
			try:
				text += " ( %s )\n" % i[2]
			except:
				text += "\n"
			c+=1
			if c == 6:
				text += "	...   ...   ...\n"
				break
		text +="=====================================================\n"
		return text

	def process(self): # AutoDock version
		self.getPoses()
		# check here we have poses
		if len(self.poses):
			if self.DEBUG: print "found poses:", len(self.poses)
			self.FastReclustering()
			self.extractResultsPoses()
			if self.doInteractions:
				self.calcInteractions()		
		else:
			if self.DEBUG: print "[ WARNING: No poses found in DLG ! ]" 


	def appendDlg(self, file_list): # allows to add a dlg to the list TODO potentially useless
		if self.dlgList:
			for d in file_list:
				if not d in self.dlgList:
					self.dlgList.append(d)
		else:
			self.dlgList = file_list


	def getPoses(self, include_hydrogens=False, ignoreWaters = False):
		"""
		- parser of docked poses from the dlg
		- populate self.atomTypes[]
		- define self.ligName (str)
		- populate poses dictionary: { "text", "coord", "energy", 'source'}
					  ( for every pose the 'source' specify dlg containing it and the position)
		"""
		mismatchMsg = ("\n\n ****************** WARNING: POSE REJECTED **************\t\t [VsResultsGenerator.py]\n" 
					   "	WARNING: %s atom count mismatch.\n"
					   "	Ligand	: %s\n"
					   "	DLG file  : %s\n"
					   " *********************************************************\n")
		#accepted_kw = [ "ATOM", "HETATM", "ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF", "REMARK"  ]
		accepted_kw = [ "ATOM", "HETATM", "TORSDOF"  ]
		atype_list_complete = False
		#water_found = False
		atomCount = None
		atomCountFlex = None
		for d in self.input_files: # we expect dlgs
			inside = False
			c = 0
			## for l in getLines(d)[50:]:
			## XXX N3P 2012.9.20 removed to deal with new 
			##		 condenset FA@H results
			for l in getLines(d):
				if l[0:7] == "DOCKED:":
					inside = True # used to stop parsing as soon as another keyword (i.e. "DPF> ") is found
					l = l[8:]
					if "MODEL" in l: # initialize pose
						text_pose = []
						true_ligand = []
						coord = []
						flex_res = []
						water = []
						c += 1
						in_res = False
					elif "ENDMDL" in l:
						coord = array( coord, 'f')
						if len(coord) == 0:
							if DEBUG: print "[ Warning! no coordinates found in this dlg: %s ]" % (d)
							break
						if flex_res:
							flex_res = getCoords(flex_res) # transform coordinates and process them
																# as we do with receptor atoms   
						if atomCount == None:
							atomCount = len(coord)
						if atomCountFlex == None:
							atomCountFlex = len(flex_res)
							
						if len(coord) == atomCount:
							if len(flex_res) == atomCountFlex:
								# flex res atoms are added to the pose here
								self.poses.append( { "text"		 : text_pose, 
													 "coord"		: coord, 
													 "energy"	   : e ,
													 'source'	   : (os.path.basename(d)+":"+str(c)), 
													 "flex_res"	 : flex_res,
													 "water_bridge" : water,
													 'true_ligand'   : true_ligand,
												   } ) 
								atype_list_complete = True
							else:
								print mismatchMsg % ( 'Flex residue', self.ligName, d) #, getLines(d).index(l) )
						else:
							print mismatchMsg % ( 'Ligand', self.ligName, d) #, getLines(d).index(l) )

					elif "BEGIN_RES" in l:  
						in_res = True
						self.flexres = True # flex_res trigger
						flex_res.append(l)
						text_pose.append(l)
					elif "END_RES" in l: 
						in_res = False
						flex_res.append(l)
						text_pose.append(l)
					elif l.startswith("ATOM") or l.startswith("HETATM"): # if the line is atom but *not* flexres
						text_pose.append(l)
						if not in_res:
							true_ligand.append(l)
						if not in_res:
							atype = l.rsplit()[-1]
							if (atype == 'HD' and include_hydrogens) or not ( atype in self.ignore_at):
								try:
									# NOTE the flex res are not considered in the clustering
									coord.append([float(l[30:38]),float(l[38:46]),float(l[46:54])])
									if not atype_list_complete: # FUTURE: to be used for the David's reclustering method
										self.atomTypes.append( atype ) 
								except:
									print ">WARNING! error in parsing coords in file :",d
									self.problematic.append(d)
									# the entire dlg is ignored if something bad happens
									# rude, but efficient on the large scale...
									break 
							elif atype in self.water_types and not ignoreWaters:
								water.append(l)
								self.isHydrated = True
						else:
							flex_res.append(l)
					elif "USER    Estimated Free Energy of Binding    =" in l:
						e = l.split("=")[1]
						e = float(e.split("k")[0])
						if e == 0.:
							if DEBUG: print "[ Warning! FEB is zero... and it shouldn't ]"
							
					elif l.split(None, 1)[0] in accepted_kw:
						text_pose.append(l)
				elif l.startswith("DPF> move") and not self.ligName:
					self.ligName = l.split("DPF> move ")[1].split(".pdbqt")[0]
				if l.startswith("DPF >") and inside: # stop reading dlg file (average lines skipped: ~30%)
					break
			self.totRuns = len(self.poses) # TODO check this indentation?
			# Debugging printing of poses
			if self.DEBUG:
				for p in self.poses:
					for f in p:
						print f
						print p[f]
						print "==================================="

			# THERE
			#if self.isHydrated and water_auto:
			#	self.processWaters()

	def checkPoses(self):
		# sf-DBG?
#		 return
#		 """ check that all poses are consistent """
#		 print "CALLED"
		
		discard = []
		correct = len(self.poses[0]['coord'])
		correctFlex = len(self.poses[0]['flex_res'])

		if self.DEBUG:
			print "> reference structure coord len [ %d ] " % correct
		for i in range(len(self.poses)-1):
			reject = True
			pose = self.poses[i+1]
			if len(pose['coord']) == correct:
				if len(pose['flex_res']) == correctFlex:
					reject = False
				else:
					if self.DEBUG or 1: print "\tPose discarded [ flexRes atoms mismatch ]"
			else:
				if self.DEBUG or 1:
					print "\tPose discarded [ ligand atoms mismatch ]", len(pose['coord']), correct
			if reject:
				discard.append(i)
				if self.DEBUG or 1: print "Pose discarded [ vect len mismatch ]"
		if self.DEBUG or 1: print "%d poses to remove" % len(discard)
		discard.sort(reverse=True)
		for d in discard:
			self.poses.pop(d)

			


	# WATER
	def processWaters(self, poses = None):
		#print "...acshuly, I don't know what to do, nao..."

		self.hydro_docking_processor = processHydroDocking(gridmap=self.water_map, 
					   ignore_types = self.ignore_at, mapdistrange = 1.0,
					   conservedwaterentropy=self.fullentropy)

		#print "DEBUG>", self.ligName
		if poses == None:
			poses = self.poses

		

		for i in range(len(poses)):
			self.hydro_docking_processor.pose = poses[i]
			poses[i] = self.hydro_docking_processor.process()
			
			#poses[i] = processHydroDocking(poses[i], gridmap = self.water_map,
			#		ignore_types = self.ignore_at, mapdistrange = 1.0)
			#poses[i] = processHydroDocking(poses[i], receptor=self.receptor, gridmap = self.water_map,
			 #	   ignore_types = self.ignore_at, mapdistrange = 1.0)
			if self.DEBUG: print "[ CALC WATER_INTERACTIONS HERE WHEN WET WILL BE READY ]"
		#for pose in poses:
		#	for f in pose:
		#		#if "water" in f :
		#		#   for w in pose[f]:
		#		#		print w.strip()
		#		print f
		#		print pose[f]
		#print "===================="
		return


	def FastReclustering(self):
		""" self.cluster results : poses [x,y,x,z,y,z,x, ...] =>  [ [x,x,x], [y,y], [ z,z,z,z,z,z,z,z], ... ] """
		self.clusters = []
		self.poses.sort(key=lambda x: x['energy'])
		poses = self.poses[:]
		while len(poses) > 1:
			rmsdcalc = RMSDCalculator(poses[0]["coord"])
			def func(x): return rmsdcalc.computeRMSD(x["coord"]) <= self.rmsTol
			cluster = filter(func, poses[1:]) + [poses[0]]
			cluster.sort(key=lambda x: x['energy'])
			# XXX	N3P the next line has weird behaviors basing on the Numpy used...
			#		<type 'exceptions.ValueError'>: 
			#		The truth value of an array with more than one element is ambiguous. Use a.any() or a.all()
			# poses = [ p for p in poses if not p in cluster] XXX REMOVED XXX
			reminders = []
			for p in poses: # TODO change this to another function as above to use 'filter'
				used=False
				for c in cluster:
					if p['text'] == c['text']: # XXX INEFFICIENT HERE!!
						used = True
						break
				if not used:
					reminders.append(p)
			poses = reminders[:]
			self.clusters.append(cluster)
		if len(poses):
			self.clusters.append(poses)

	def extractResultsPoses(self):
		"""find LE and LC (if available) and the self.histogram"""
		# TODO speed here can be improved
		LCpose = ""
		LCsize = -1
		# LCenergy = 10000. Changed for catastrophyc results with energy bigger than that
		# (i.e. too-small grid box)
		LCenergy = 1E99
		LEpose = ""
		LEsize = -1
		LEenergy = 1E99
		c = 0


		for pop in self.clusters:
			curr_csize = len(pop)
			curr_cenergy = pop[0]['energy']
			self.histogram.append([curr_cenergy, curr_csize])
			if (curr_csize > LCsize) or (curr_csize == LCsize and curr_cenergy < LCenergy ):
				LCsize = curr_csize
				LCpose = pop[0]
				LCenergy = curr_cenergy
				LCindex = c
			if curr_cenergy < LEenergy:
				LEenergy = curr_cenergy
				LEpose = pop[0]
				LEsize = curr_csize
				LEindex = c
			c += 1
		self.results.append(LEpose)
		self.results[0]['csize'] = LEsize
		self.results[0]['cpercent'] = (float(LEsize) / float(self.totRuns))*100
		self.results[0]['leff'] = self.results[0]['energy'] / float( len( self.atomTypes))
		self.results[0]['vdw_contacts'] = []
		self.results[0]['metal_coord'] = [] 
		acc, don = self.findHbAccepDon( self.results[0]['text'] )
		self.results[0]['hba_atoms'] = acc
		self.results[0]['hbd_atoms'] = don
		if self.DEBUG:
# 		   writeList(self.ligName+"_lig_acceptors.pdb", acc)
# 		   writeList(self.ligName+"_lig_donors.pdb", don)
			print self.ligName+"_lig_acceptors.pdb", acc
			print self.ligName+"_lig_donors.pdb", don
		self.histogram[LEindex].append("**") # <-marker for LE
		if not LEpose == LCpose:
			self.results.append( LCpose )
			self.results[1]['csize'] = LCsize
			self.results[1]['cpercent'] =  (float(LCsize) / float(self.totRuns))*100
			self.results[1]['leff'] = self.results[1]['energy'] / float( len( self.atomTypes))
			self.results[1]['vdw_contacts'] = []
			self.results[1]['metal_coord'] = []
			acc, don = self.findHbAccepDon(self.results[1]['text'])
			self.results[1]['hba_atoms'] = acc
			self.results[1]['hbd_atoms'] = don
			self.histogram[LCindex].append("*") # <-marker for LC
		self.histogram.sort(key = lambda x: x[0])

	def generatePDBQTplus(self):
		# AutoDock version
		""" creates a formatted AutoDock PDBQT plus file 
			with all the data extracted from the re-clustering
		"""
		time_info = time.localtime()[0:6]
		# pack the header:
		buff  = "USER    ADVS_result> %d-%d-%d %d:%d:%d\n"  % (time_info)
		buff += "USER    AD_rec> %s\n" % ( self.recname )
		buff += "USER    AD_runs,rmstol,tot_clusters> %d,%1.2f,%d\n" % ( self.totRuns,\
					self.rmsTol, len(self.histogram))
		buff += "USER    AD_dlg_list> " 
		for d in self.input_files: # we expect dlgs here
			buff+="%s," % (os.path.basename(d))
		buff = buff[:-1]+"\n"
		buff += "USER    AD_results> %d\n" % ( len(self.results) ) # add the results count  (1,2)
		buff += "USER    AD_histogram> " # add the histogram
		for h in self.histogram:
			buff += "%2.3f:%d" % (h[0], h[1] )
			try:
				buff += ":"+h[2]+","
			except:
				buff += ","
		buff = buff[:-1]+"\n" # cut out the last ","
		for p in range(len(self.results)): # separate pose MODELS are handled here
			pose = self.results[p]
			buff += "MODEL   %d\n" % ( p+1)
			buff += "USER    #	 energy,\tleff,\tc_size,\tc_pc\n"
			buff += "USER    AD_%s> %2.3f,\t%2.3f,\t%d,\t%3.2f\n" % (self.tag[p],\
									pose['energy'], pose['leff'], pose['csize'], pose['cpercent']) 
			#
			# XXX change the ligand atom numbering with the code Michel sent 3.25.2012
			#

			if self.interactions:
				if pose['hb']['acceptors']:		 # add hb acceptors info
					buff += "USER    AD_%s_hba> " % (self.tag[p])
					for pair in pose['hb']['acceptors']:
						buff += "%s~~%s," % (pmvAtomStrip(pair[0]), pmvAtomStrip(pair[1]))
					buff = buff[:-1]+"\n"
				if pose['hb']['donors']:			# add hb donors info
					buff += "USER    AD_%s_hbd> " % (self.tag[p])
					for pair in pose['hb']['donors']:
						buff += "%s~~%s," % (pmvAtomStrip(pair[0]), pmvAtomStrip(pair[1]))
					buff = buff[:-1]+"\n"
				if pose['vdw_contacts']:			# add close contacts info
					buff += "USER    AD_%s_vdw> " % (self.tag[p])
					for pair in pose['vdw_contacts']:
						#buff += "%s," % (pmvAtomStrip(a))
						buff += "%s~~%s," % (pmvAtomStrip(pair[0]), pmvAtomStrip(pair[1]))
					buff = buff[:-1]+"\n"
				if pose['metal_coord']:			 # add metal coordination info
					buff += "USER    AD_%s_mtl> " % (self.tag[p])
					for pair in pose['metal_coord']:
						buff += "%s~~%s," % (pmvAtomStrip(pair[0]), pmvAtomStrip(pair[1]))
					buff = buff[:-1]+"\n"
				if pose['pi']['p-stack']:
					buff += "USER    AD_%s_ppi> " % (self.tag[p])
					for a in pose['pi']['p-stack']: # add pi interaction info (t-stack)
						buff += "%s~~(%2.3f,%2.3f,%2.3f:%2.3f,%2.3f,%2.3f);" %\
						(a[0], a[1][0], a[1][1], a[1][2], a[2][0], a[2][1], a[2][2])
					buff = buff[:-1]+"\n"
				if pose['pi']['t-stack']:
					buff += "USER    AD_%s_tpi> " % (self.tag[p])
					for a in pose['pi']['t-stack']: # add pi interaction info (p-stack)
						buff += "%s~~(%2.3f,%2.3f,%2.3f:%2.3f,%2.3f,%2.3f);" %\
						(a[0], a[1][0], a[1][1], a[1][2], a[2][0], a[2][1], a[2][2])
					buff = buff[:-1]+"\n"
			buff += "USER    AD_%s_source> %s\n" % (self.tag[p], pose['source'])
			for i in pose['text']:
				if not i.startswith("MODEL") and not i.startswith("ENDMDL"):
					buff += i
			buff +="ENDMDL  %d\n" % ( p+1 )
		return buff
#
################################################## AutoDock class [end]

class AutoDockVinaVsResult(VSDockingResultsGenerator):
	def __init__(self, input_files = None, mode = 1, water_map = None, # XXX useless, but kept for compatibility
				recname = None, receptor = None, auto = True, 
				doInteractions = True, hbtol = 0.0, debug = False):


		VSDockingResultsGenerator.__init__(self,
			input_files = input_files,		   # list of vina ligand_out.pdbqt from wich poses will be extracted
			recname = recname,		 # receptor name (if not provided it will be taken from 'self.receptor'
			receptor = receptor,		# receptor coordinates: PDBQT or getCoord() output (CACHE)
			doInteractions = doInteractions,  # calculate ligand-target interactions
			hbtol = 0.0,
			debug = debug,		  # debug mode
			)
		self.totRuns = 0 # Vina poses in the file
		self.mode = 1	# poses saved in the result file
		self.histogram = [] # energy | lb RMSD | ub RMSD
		if not input_files == None:
			name = os.path.basename(input_files) # the name of the ligand is guessed from the filename
			name = os.path.splitext(name)[0]  # name can be eventually changed by
											  # overriding the self.ligName value
			self.ligName = name
			if auto:
				self.process()

	def __str__(self):
		text = "\n===================================================\n"
		text +=" This is my VSDockingResults object. [ Vina ]\n"
		text +=" There are many like it, but only this one\n"
		text +=" has been made from => %s\n" % self.ligName
		text +="  - receptor name			: %s\n" % self.recname
		text +="  - total binding modes	  : %s\n" % self.totRuns
		text +="  - results saved			: %s\n" % len(self.results)
		text +="  - flex res				 : %s" % self.flexres
		text +="\n Result(s) ------------------------------------------\n"
		if self.interactions:
			text += "	   Energy	L.eff  vdW	  HB	Pi-Pi\n"
			for i in range(len(self.results)):
				text += "  %s   %2.3f\t%2.3f\t%d\t%d\t%d\n" % \
					( i+1, self.results[i]['energy'], self.results[i]['leff'],
					  len( self.results[i]['vdw_contacts']),
					  (len(self.results[i]['hb']['acceptors'])+len(self.results[i]['hb']['donors'])),
					  len(self.results[i]['pi']['p-stack'])+len(self.results[i]['pi']['t-stack'])
					)
		else:
			text += "	   Energy	L.eff\n"
			for i in range(len(self.results)):
				text += " %s   %2.3f\t%2.3f" % \
					( i+1, self.results[i]['energy'], self.results[i]['leff'])

		text +=" Histogram ------------------------------------------\n"
		c = 1
		text +=" Pose	E	 lbRMS	 ubRMS\n"
		for i in self.histogram:
			# it should be # energy | lb RMSD | ub RMSD
			text+= " %d   %2.3f |  %2.3f | %2.3f\n" % (c, i[0],i[1],i[2])
			c+=1
			if c == 9:
				text += "		...   ...   ...\n"
				break
		text +="=====================================================\n"
		return text

	def process(self): # AutoDock Vina version
		self.getPoses()
		self.extractResultsPoses()
		if self.doInteractions:
			self.calcInteractions()

	def getPoses(self, include_hydrogens=False):
		"""
		- parser of docked poses from a Vina output
		- populate self.atomTypes[]
		- populate poses dictionary: { "text", "coord", "energy", 'source'}
					  ( for every pose the 'source' specify dlg containing it and the position)
		"""
		mismatchMsg = ("\n\n ****************** WARNING: POSE REJECTED **************\t\t [VsResultsGenerator.py]\n" 
					   "	WARNING: %s atom count mismatch.\n"
					   "	Ligand	: %s\n"
					   "	DLG file  : %s\n"
					   " *********************************************************\n")
		atomCount = None
		atomCountFlex = None

		accepted_kw = [ "ATOM", "HETATM", "ROOT", "ENDROOT", "BRANCH", "ENDBRANCH", "TORSDOF", "REMARK"  ]
		atype_list_complete = False
		inside = False
		c = 0
		lines = getLines(self.input_files)
			
		if ((config.Reprocess and "ADVS_Vina_result" not in lines[0]) or \
			(not config.Reprocess and "REMARK VINA RESULT" not in lines[1])):
			print "getPoses: WARNING! neither Vina nor ADVS output PDBQT file :",self.input_files
			self.poses = []
			return
	
		if DEBUG: print "[found Vina result pdbqt]"
		nskip = 0
		nempty = 0
		defaultEnergy = 1000. 
		
		for il, l in enumerate(lines):
			# print il,l
			try:
				if len(l.strip()) == 0:
					nempty += 1
					continue
					
				if "MODEL" in l: # initialize pose
					text_pose = []
					true_ligand = []
					coord = []
					flex_res = []
					c += 1
					in_res = False
					e = defaultEnergy
				elif "ENDMDL" in l:
					coord = array( coord, 'f')
					if flex_res:
						flex_res = getCoords(flex_res) # transform coordinates and process them
													   # as we do with receptor atoms   
					if atomCount == None:
						atomCount = len(coord)
					if atomCountFlex == None:
						atomCountFlex = len(flex_res)
					if len(coord) == atomCount:
						if len(flex_res) == atomCountFlex:
							# flex res atoms are added to the pose here
							if e == defaultEnergy:
								print "getPoses: WARNING! e not set line=%d?!" % (il)
							self.poses.append( { "text" : text_pose, "coord" : coord, "energy" : e ,\
								'source' : (os.path.basename(self.input_files)+":"+str(c)),
								'true_ligand': true_ligand, "flex_res" : flex_res } ) 
							atype_list_complete = True
						else:
							# 150925: what is `d`?
							print mismatchMsg % ( 'Flex residue', self.ligName, 'd') #, getLines(d).index(l) )
					else:
						print mismatchMsg % ( 'Ligand', self.ligName, 'd') #, getLines(d).index(l) )
						  
				elif l.startswith("BEGIN_RES"):
					in_res = True
					self.flexres = True # flex_res trigger
					flex_res.append(l)
					text_pose.append(l)
				elif "END_RES" in l: 
				# elif l.startswith("END_RES"): # sf_150818 ??
					in_res = False
					flex_res.append(l)
					text_pose.append(l)
				elif l.startswith("ATOM") or l.startswith("HETATM"):
					text_pose.append(l)
					if not in_res:
						true_ligand.append(l)
					if not in_res:
						atype = l.rsplit()[-1]
						if (not atype == "HD") or include_hydrogens : # HD should be excluded from array/RMSD calculation
							try:
								coord.append([float(l[30:38]),float(l[38:46]),float(l[46:54])])
								if not atype_list_complete: # FUTURE: to be used for the David's reclustering method
									self.atomTypes.append( atype ) 
							except:
								print "getPoses: WARNING! error in parsing coords in file :",self.input_files
								self.problematic.append(self.input_files)
								break 
					else:
						flex_res.append(l)
				elif (not config.Reprocess and "REMARK VINA RESULT:" in l):
					values = l.split(":")[1]
					values = values.split()
					e = float(values[0])
					lbrms = float(values[1])
					ubrms = float(values[2])
					self.histogram.append([e, lbrms, ubrms])
				elif (config.Reprocess and "ADVina_pose1>" in l):
					# USER    #	 energy,	leff
					# USER    ADVina_pose1> -8.700,	-0.395
					gpos = l.find('>')
					postPrefix = l[gpos+1:]
					bits = re.findall(FloatRegEx, postPrefix)
					# bits = [b for b in l.split(' \t') if len(b) > 0]
					# bits2 = bits[2].split('\t')
					efield = bits[0]
					e = float(efield)
					# 2do: LB, UB from processed file?
					lbrms = 0. # was float(values[1])
					ubrms = 0. # was float(values[2])
					self.histogram.append([e, lbrms, ubrms])
					
				elif l.split(None, 1)[0] in accepted_kw:
					text_pose.append(l)
				else:
					nskip += 1
			except Exception, e:
				print 'getPoses: exception',e
				
		self.totRuns = len(self.poses)
		if DEBUG: print "[NPose=%d NEmpty=%d NSkip=%d]" % (len(self.poses), nempty, nskip)
	
	def extractResultsPoses(self):
		"""extract the best result (the first?) accordingly to Vina"""
		# TODO figure out if we want some more from the result?
		# so far it is good in this way...
		if self.mode > self.totRuns:
			self.mode = self.totRuns
		for i in range(0, self.mode):
			self.results.append(self.poses[i])
			self.results[i]['leff'] = self.results[i]['energy'] / float( len( self.atomTypes))
			self.results[i]['vdw_contacts'] = []
			self.results[i]['metal_coord'] = [] # TODO 
			acc, don = self.findHbAccepDon( self.results[i]['true_ligand'] )
			self.results[i]['hba_atoms'] = acc
			self.results[i]['hbd_atoms'] = don
		if self.DEBUG:
# 		   writeList(self.ligName+"_lig_acceptors.pdb", acc)
# 		   writeList(self.ligName+"_lig_donors.pdb", don)
			print self.ligName+"_lig_acceptors.pdb", acc
			print self.ligName+"_lig_donors.pdb", don
#			writeList(self.ligName+"_BABABA_lig_acceptors.pdb", acc, addNewLine=1)
#			writeList(self.ligName+"_BABABA_lig_donors.pdb", don, addNewLine=1)

	def generatePDBQTplus(self):
		# Vina version
		""" creates a formatted AutoDock PDBQT plus file 
			with all the data extracted from the re-clustering
		"""
		time_info = time.localtime()[0:6]
		# pack the header:
		buff  = "USER    ADVS_Vina_result> %d-%d-%d %d:%d:%d\n"  % (time_info)
		buff += "USER    ADVina_rec> %s\n" % ( self.recname )
		buff += "USER    ADVina_poses> %d\n" % ( self.totRuns)
		buff += "USER    ADVina_input_file> %s\n" % (os.path.basename(self.input_files))
		buff += "USER    ADVina_results> %d\n" % ( len(self.results) ) # add the results count  (1,2,...)
		buff += "USER    ADVina_histogram> " # add the histogram
		for h in self.histogram:
			buff += "%2.3f:%2.3f:%2.3f," % (h[0], h[1], h[2] )
		buff = buff[:-1]+"\n" # cut out the last ","
		for p in range(len(self.results)): # separate pose MODELS are handled here
			pose = self.results[p]
			buff += "MODEL   %d\n" % ( p+1)
			buff += "USER    #	 energy,\tleff\n"
			buff += "USER    ADVina_pose%d> %2.3f,\t%2.3f\n" % (p+1,\
									pose['energy'], pose['leff'])
			if self.interactions:
				if pose['hb']['acceptors']:	 # add hb acceptors info
					buff += "USER    ADVina_pose%d_hba> " % (p+1)
					for pair in pose['hb']['acceptors']:
						buff += "%s~~%s," % (pmvAtomStrip(pair[0]), pmvAtomStrip(pair[1]))
					buff = buff[:-1]+"\n"
				if pose['hb']['donors']:		# add hb donors info
					buff += "USER    ADVina_pose%d_hbd> " % (p+1)
					for pair in pose['hb']['donors']:
						buff += "%s~~%s," % (pmvAtomStrip(pair[0]), pmvAtomStrip(pair[1]))
					buff = buff[:-1]+"\n"
				if pose['vdw_contacts']:		# add close contacts info
					buff += "USER    ADVina_pose%s_vdw> " % (p+1)
					for pair in pose['vdw_contacts']:
						buff += "%s~~%s," % (pmvAtomStrip(pair[0]), pmvAtomStrip(pair[1]))
					buff = buff[:-1]+"\n"
				if pose['metal_coord']:		 # add metal coordination info
					buff += "USER    ADVina_pose%s_mtl> " % (p+1)
					for pair in pose['metal_coord']:
						buff += "%s~~%s," % (pmvAtomStrip(pair[0]), pmvAtomStrip(pair[1]))
					buff = buff[:-1]+"\n"   
				if pose['pi']['p-stack']:	   # add pi interaction info
					buff += "USER    ADVina_pose%s_ppi> " % (p+1)
					for a in pose['pi']['p-stack']:
						buff += "%s~~(%2.3f,%2.3f,%2.3f:%2.3f,%2.3f,%2.3f);" %\
						(a[0], a[1][0], a[1][1], a[1][2], a[2][0], a[2][1], a[2][2])
					buff = buff[:-1]+"\n"
				if pose['pi']['t-stack']:	   # add pi interaction info
					buff += "USER    ADVina_pose%s_tpi> " % (p+1)
					for a in pose['pi']['t-stack']:
						buff += "%s~~(%2.3f,%2.3f,%2.3f:%2.3f,%2.3f,%2.3f);" %\
						(a[0], a[1][0], a[1][1], a[1][2], a[2][0], a[2][1], a[2][2])
					buff = buff[:-1]+"\n"
			buff += "USER    ADVina_pose%s_source> %s\n" % (p+1, pose['source'])
			for i in pose['text']:
				if not i.startswith("MODEL") and not i.startswith("ENDMDL"):
					buff += i
			buff +="ENDMDL  %d\n" % ( p+1 )
		return buff

if __name__ == '__main__':

	HostName = socket.gethostname()
	if HostName == 'mgl0':
		print 'running on mgl0, good!'
		BaseDir = '/export/wcg/'
	
	elif HostName == 'mgl3':
		print 'running on mgl3, slow(:'
		BaseDir = '/mgl/storage/wcg/'

	elif HostName.startswith('hancock'):
		print 'running local on hancock'
		BaseDir = '/Data/sharedData/coevol-HIV/WCG/'

	elif HostName.startswith('mjq'):
		print 'running local on mjq'
		BaseDir = '/home/Data/coevol-HIV/WCG/'

	else:
		print 
		sys.exit( ('unknown host %s' % (HostName)) )
	
	# 160411
	# NNRTI_FDA_x3LAN x3LAN False True
	
	# 160405
	# DUDE-HIVRT-x3LAN x3LAN False True
	# DUDE-HIVPR-x1XL2 x1XL2 False True
	
	# 160523: Local (hancock) reprocessing of Exp96-LEDGF
	# "['Exp96-LEDGF']" True False "['x3AO1_IN_LEDGF', 'x3AO2_IN_LEDGF']"
	
	if len(sys.argv) < 4:
		sys.exit('process: missing arguments: exptName, reprocess, makeBatches ?!')		

	ExptNameList = eval(sys.argv[1])
	config.Reprocess = eval(sys.argv[2])  # rerunning process over processed (.VS) file vs. original dockings
	makeBatches = eval(sys.argv[3])

	print 'process: arguments=%s' % (sys.argv)

	if len(sys.argv) ==5:
		RestrictReceptList = eval(sys.argv[4])
	else:
		RestrictReceptList = []
		
	verbose = False
			
	batchSize = 10000

	# 160404: Reprocess DUDE
# 	DockDir = BaseDir + 'processed/DUDE_160404/'
# 	StructDir = BaseDir + 'processed/DUDE_160404/structure/'
# 	global Reprocess
# 	Reprocess = True
# 	makeBatches = True
# 	batchSize = 10000
	
	# process defaults

	config.procVersion = '0.0'  # 0.2 for crawl_ADV.visit_ADV_tgz(), visit_ADV()
	# 160216: consuming SForli's rerun of DUDE
	# ZINC63694490_out_Vina_VS.pdbqt

	# raccoon defaults
	DEBUG = False
	rmsTol = 2.0
	mode = 1 # binding modes extracted from the result 
	header = "#name\tenergy\tligand_efficiency\ttotal_poses\treceptor\tfilename\n"
	suffix = "_Vina_VS"
	hbtol = 0.0 
	
	if config.Reprocess:
		# 160404: Reprocess non-HIV DUDE to include vdwP
		dockPattern = '*_out_Vina_VS.pdbqt'
	else:
		dockPattern = '*_out.pdbqt' 

	print '<process "%s">' % (ExptNameList)
	print '\thost=%s' % (HostName)
	print '\treprocess=%s\n\tmakeBatches=%s\n\tbatchSize=%d' % (config.Reprocess,makeBatches,batchSize)
	print '\tRestrictReceptList=%s' % (RestrictReceptList)

	for ExptName in ExptNameList:
		
		if config.Reprocess:
			DockDir = BaseDir + 'processed/' + ExptName + '/'
			StructDir = BaseDir + 'structure/'
		else:
			DockDir = BaseDir + 'docking/' + ExptName + '/'
			StructDir = DockDir
			
		ProcDir = BaseDir + 'process2/' + ExptName + '/'	
	
		begTime = datetime.datetime.now()
		begTimeStr = begTime.strftime('%y%m%d_%H%M%S')

		print '<processExpt "%s" %s>' % (ExptName,begTimeStr)
	
		resultDirList = glob(DockDir+'Results_*')
		resultDirList.sort()
					
		print 'process: NResultDir=%d' % (len(resultDirList))
		for ird,rdir in enumerate(resultDirList):
			# rdir = /Data/sharedData/coevol-HIV/WCG/processed/Exp96_LEDGF/Results_x3AO1_IN_LEDGF
			spos = rdir.rfind('/')
			dirName = rdir[spos+1:]
			upos = dirName.find('_')
			receptName = dirName[upos+1:]
			
			if len(RestrictReceptList) > 0 and receptName not in RestrictReceptList:
				print 'process: Expt=%s RdirIdx=%d Recept=%s not in RestrictReceptList; skipping' % (ExptName,ird,receptName)
				continue
			
			# Allow partial processing runs; jump over previous runs at RECEPTOR level
			outReceptDir = ProcDir + '%s/' % (receptName)
			if os.path.isdir(outReceptDir):
				print 'process: Expt=%s RdirIdx=%d Recept=%s ALREADY DONE; continuing' % (ExptName,ird,receptName)
				continue
			
			# Cache all structures
			structReceptPath = StructDir + receptName + '.pdbqt'
			
			totFiles = 0
			prevReceptor = None
			generator = None
				
			fileList = glob(rdir+'/FAHV_*_processed.tgz')
			fileList.sort()
	
			print 'process: Expt=%s RdirIdx=%d Recept=%s NBatch=%d' % (ExptName,ird,receptName,len(fileList))
			for tpi,tgzPath in enumerate(fileList):
				# tgzPath = '/Data/sharedData/coevol-HIV/WCG/processed/Exp120/Results_x3KF0_prASw0c0/FAHV_x3KF0_prASw0c0_0550347_processed.tgz'
	
				searchDir = tempfile.mkdtemp()
				allTar = tarfile.open(tgzPath)
				print 'process: Extracting %s' % (tgzPath)		
				allTar.extractall(searchDir)
				tgzBits = tgzPath.split('/')
				tballName = tgzBits[-1]
				tballBits = tballName.split('_')
				
				# ASSUME batchno is last part of name prior to "_processed"
				batchName = tballBits[-2]
				batch = int(batchName)
	
				ppos = tballName.find('.')
				tballName = tballName[:ppos]
				tballPath = searchDir + '/' + tballName
				input_files = pathToList(tballPath, recursive = False, pattern = dockPattern)
				
				if makeBatches:
					print 'process: batching %d dockings into %d batches' % (len(input_files),batchSize)
					nbatchLig = 0
					batch = 0
		
				outDir = ProcDir + '%s/batch_%07d/' % (receptName,batch) 
								   
				if not os.path.isdir(outDir):
					os.makedirs(outDir)
					
				# Cache structures
				# NB: use changed structReceptPath as indication that we need a new generator
				procReceptPath = tballPath + '/' + receptName + '.pdbqt'
				if not os.path.exists(procReceptPath):
					sys.exit('process: Missing receptor?! '+procReceptPath)
	
				## copy receptor to outDir
				shutil.copy(procReceptPath, outDir)
				 
				protRoot = structReceptPath
				if protRoot != prevReceptor:
	
					if not os.path.exists(structReceptPath):
						print 'process: Caching structure',receptName
						shutil.copy(procReceptPath,structReceptPath)
											
					print 'process: Updating current receptor for %s with %s %s' % (ExptName,protRoot,structReceptPath)
							  
					prevReceptor = protRoot
					if not generator == None:
						try:
							bht = generator.rec_bht
							freeBHtree(bht)
						except Exception,e:
							print 'process: cant free bht?',e
		 
					generator =  AutoDockVinaVsResult(input_files = None, 
						mode = mode,
						receptor = structReceptPath, 
						recname = protRoot, 
						auto = False, 
						doInteractions = True,
						hbtol = hbtol)
					print 'process: bht built natoms=%d' % (len(generator.rec_bht_indices))
										
				ndock = len(input_files)
				print 'process: Expt=%s RdirIdx=%d TPI=%d NDock=%d' % (ExptName,ird,tpi,ndock)
				
				for fi,l in enumerate(input_files):
							 
					try:
						path_root = os.path.dirname(l) 
						if not path_root:
							path_root = os.getcwd()
			
						if makeBatches:
							if nbatchLig >= batchSize:
								batch += 1
								print 'process: new batch %d nlig=%d' % (batch,fi)
								nbatchLig = 1
							else:
								nbatchLig += 1
													 
						if DEBUG: print "processing", l
						# NB:    AutoDockVinaVsResult.setLigands() assumes l is single STRING
						# contra AutoDockVsResult.setLigands(), which assumes l is a list!
						generator.setLigands(l)
							 
						# NB: need to explicitly "guess" ligName as in AutoDockVinaVsResult.__init__
						# since generator isn't being reconstructed for each ligand!
						ligPath = os.path.splitext(l)[0]
						lpbits = ligPath.split('/')
						ligName = lpbits[-1]
						generator.ligName = ligName
							 
						generator.process()
						pdbqt = generator.generatePDBQTplus()
																 
						if config.Reprocess: # or generator.ligName.endswith(suffix):
							output_filename = outDir + generator.ligName+'.pdbqt' 
						else:
							output_filename = outDir + generator.ligName+suffix+'.pdbqt' 
						fp = open(output_filename, 'w')
						fp.write(pdbqt)
						fp.close()	 
				 
					except Exception,e:
						print "ERROR: problems processing input :"
						print "file :", l
						print "error:", e
							 
					if verbose and fi % 1000 == 0: print "NInputFile=%d" % (fi)
					 
					# eo input_files loop
				
				# eo tgz batch loop
				# NB: need to remove even tempdir, or risk overfilling /tmp!
				shutil.rmtree(searchDir, True)
				totFiles += ndock
				elapTime = datetime.datetime.now() - begTime
				print '<process Expt=%s RdirIdx=%d TPI=%d TotDock=%d %s sec>' % (ExptName,ird,tpi,totFiles,elapTime.total_seconds)
			
			elapTime = datetime.datetime.now() - begTime
			print '<process ResultDir=%d %s Totdock=%d %s sec>' % (ird,receptName, totFiles, elapTime.total_seconds)
	
		elapTime = datetime.datetime.now() - begTime
		print '</processExpt: Expt=%s done %s sec>' % (ExptName,elapTime.total_seconds)
		
	print 'process: done!' 