'''
Replacement file for pseuso. Now using pymatgen.
'''
import numpy
import pymatgen
import getdata
import pymatgen.analysis.structure_matcher
import pymatgen.symmetry.analyzer
import spglib as spg

def is_symmop(p1, op, tol=.01):
	p2 = op.operate(p1)
	dsquared = (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2
	if dsquared <= tol**2: return True
	else: return False

def dsquared(p1, p2):
	return (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2

#Obtain a set of points by applying a set of operations to one point
def applyops(p1, ops):
	pos = []
	pos.append(pymatgen.core.operations.SymmOp.from_xyz_string('x,y,z').operate(p1))
	for x in ops:
		new = x.operate(p1)
		isnew = True
		for old in pos:
			if dsquared(old,new) < .0001:
				isnew = False
		if isnew: pos.append(new)
	return pos

#For a set of operations, return the set of symmetry operations for a point
def stabilizers(p1, ops, tol=.01):
	stab = []
	for i in range(len(ops)):
		p2 = ops[i].operate(p1)
		dsquared = (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2
		if dsquared <= tol**2:
			stab.append(ops[i])
	return stab

#Find the set of general Wyckoff positions corresponding to a point
def wyckset(p1, ops, tol=.01):
	pos = []
	pos.append(pymatgen.core.operations.SymmOp.from_xyz_string('x,y,z').operate(p1))
	wset = []
	wset.append(pymatgen.core.operations.SymmOp.from_xyz_string('x,y,z'))
	for x in ops:
		new = x.operate(p1)
		isnew = True
		for old in pos:
			if dsquared(old,new) < .0001:
				isnew = False
		if isnew:
			pos.append(new)
			wset.append(x)
	return wset	

#def settransform(t1, t2):
	#
	#Apply transform t1 to t2, then normalize
	#Try using class pymatgen.core.sites.PeriodicSite()

#import our base structure from a cif file
#struct1 will be stored as a pymatgen.core.structure Structure class object

'''If you want to access individual points from the Structure object,
use struct1.sites[i].frac_coords
This gives a length 3 array of floating point fractional coordinates.
For absolute coordinates, use sites[i].coords instead.
'''
struct1 = getdata.fromFile("test.cif")
#Get a list of symmetry operations for the structure
symmops = pymatgen.symmetry.analyzer.SpacegroupAnalyzer(struct1).get_symmetry_operations()

#test wyckset()
point = [0.1,0.2,0.3]
for x in wyckset(point, symmops):
	print(x.as_xyz_string())

#Find the supergroups of struct1
#Extract the general positions
