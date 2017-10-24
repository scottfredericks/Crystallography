'''
Replacement file for pseuso. Now using pymatgen.
'''
import numpy
import pymatgen
import get_data
import pymatgen.analysis.structure_matcher
import pymatgen.symmetry.analyzer
import spglib as spg
import pymatgen.symmetry.settings
import numpy as np

def filter_site(v):
	w = v
	for i in range(len(w)):
		while w[i]<0: w[i] += 1
		while w[i]>=1: w[i] -= 1
	return w

def ltransform(b, a):
	for i in range(4):
		for j in range(3):
			if b[i][j] is not float: b[i][j] = float(b[i][j])
	for i in range(3):
		if a[i] is not float: a[i] = float(a[i])
	c = [float(0), float(0), float(0)]
	i = 0
	while i < 3:
		j = 0
		while j < 3:
			c[i] += b[i][j]*a[j]
			j += 1
		c[i] += b[3][i]
		i += 1
	return c

def is_symmop(p1, op, tol=.01):
	p2 = op.operate(p1)
	dsquared = (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2
	if dsquared <= tol**2: return True
	else: return False

#square of the distance between two points
def dsquared(p1, p2):
	return (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2

#Obtain a set of points by applying a set of operations to one point
def apply_ops(p1, ops):
	#convert to np.unique(array, axis=0)
	#use filter_site()
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

#Return the subset of ops which leaves the point p1 invariant
def stabilizer(p1, ops, tol=.01):
	stab = []
	for i in range(len(ops)):
		p2 = ops[i].operate(p1)
		dsquared = (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2
		if dsquared <= tol**2:
			stab.append(ops[i])
	return stab

#Find the subsset of general Wyckoff positions corresponding to a point
#Broken: not always a conjugacy class
def wyckoff_set(p1, ops, tol=.01):
	pos = []
	pos.append(pymatgen.core.operations.SymmOp.from_xyz_string('x,y,z').operate(p1))
	wset = []
	wset.append(pymatgen.core.operations.SymmOp.from_xyz_string('x,y,z'))
	for x in ops:
		new = x.operate(p1)
		isnew = True
		for old in pos:
			if dsquared(old,new) < tol**2:
				isnew = False
		if isnew:
			pos.append(new)
			wset.append(x)
	return wset

def compose_ops(op1,op2):
	#Apply transform op1 to op2, then normalize
	m = np.dot(op1.rotation_matrix,op2.rotation_matrix)
	v = op1.operate(op2.translation_vector)
	v = filter_site(v)
	new = pymatgen.core.operations.SymmOp.from_rotation_and_translation(m,v)
	return new

#map a Wyckoff set W1 to a supergroup Wyckoff set W2
def wyckoff_set_map(W1, W2):
	if len(W2)<len(W1): print("Error: cannot map to a smaller set")
	map = []
	for i in range(len(W1)):
		for j in range(len(W2)):
			if W1[i] == W2[j]:
				map.append(j)
				break
			
		if map[i] == []:
			print("Error: could find map.")
			break
	return map

'''
def checkpseudosymmetry(species, positions, transform, generator):
	for i in len(atoms):
		a
	return
'''	

#import our base structure from a cif file
#struct1 will be stored as a pymatgen.core.structure Structure class object
'''If you want to access individual points from the Structure object,
use struct1.sites[i].frac_coords
This gives a length 3 array of floating point fractional coordinates.
For absolute coordinates, use sites[i].coords instead.'''
mystruct1 = get_data.from_file("test.cif")

#Get a list of symmetry operations for the structure
sga = pymatgen.symmetry.analyzer.SpacegroupAnalyzer(mystruct1)
symmops = sga.get_symmetry_operations()

#get atomic specie info
#Find the supergroups of struct1
#Get the transformation W,w from G to H
#Find the general positions of G


#-----------------------------------------------------------------------
#test functionality
text = ['x,y,z', '-y,x-y,z', '-x+y,-x,z', '-x,-y,z+1/2',
		'y,-x+y,z+1/2', 'x-y,x,z+1/2', 'y,x,-z+1/2', 'x-y,-y,-z+1/2',
		'-x,-x+y,-z+1/2', '-y,-x,-z', '-x+y,y,-z', 'x,x-y,-z',
		'-x,-y,-z', 'y,-x+y,-z', 'x-y,x,-z', 'x,y,-z+1/2',
		'-y,x-y,-z+1/2', '-x+y,-x,-z+1/2', '-y,-x,z+1/2', '-x+y,y,z+1/2',
		'x,x-y,z+1/2', 'y,x,z', 'x-y,-y,z', '-x,-x+y,z']
group193 = []
myop = pymatgen.core.operations.SymmOp.from_xyz_string('-x, -x+y, 1/2-z')
newops = []
for i in range(len(symmops)):
	newops.append(symmops[i])
	newops.append(compose_ops(myop,symmops[i]))
for i in range(len(mystruct1.frac_coords)):
	print(str(mystruct1.species[i])+" "+str(mystruct1.frac_coords[i]))
	for y in symmops:
		if is_symmop(mystruct1.frac_coords[i], y):
			print(y.as_xyz_string())

'''This gives a list of symmetry operations which leave a given point invariant.
It should be possible to implement this into Wyckset in order to create a unique
set of general positions corresponding to each point.'''

