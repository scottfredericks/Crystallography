'''
Replacement file for pseuso. Now using pymatgen.
'''
import pymatgen
import pymatgen.analysis.structure_matcher
import pymatgen.symmetry.analyzer
import pymatgen.symmetry.settings
import numpy as np
from timeit import default_timer as timer
import make_sitesym

#Function and class definitions
#--------------------------
def filter_site(v):
	w = v
	for i in range(len(w)):
		while w[i]<0: w[i] += 1
		while w[i]>=1: w[i] -= 1
	return w

def is_symmop(p1, op, tol=.01):
	p2 = op.operate(p1)
	dsquared = (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2
	if dsquared <= tol**2: return True
	else: return False

#square of the distance between two points
def dsquared(p1, p2):
	return (p2[0]-p1[0])**2 + (p2[1]-p1[1])**2 + (p2[2]-p1[2])**2

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

def site_stabilizer(site1, ops, tol=.01):
	#Return the subset of ops which leaves the PeriodicSite site1 invariant
	stab = []
	for i in range(len(ops)):
		site2 = pymatgen.core.sites.PeriodicSite(site1.specie, ops[i].operate(site1.frac_coords), site1.lattice)
		if site1.distance(site2) <= tol:
			stab.append(ops[i])
	return stab

def op_stabilizer(op1, ops, tol=.01):
	#return the subset of ops which leaves a SymmOp object invariant
	#used for finding site symmetry of special Wyckoff position elements
	stab = []
	for i in range(len(ops)):
		op2 = compose_ops(ops[i], op1)
		if np.linalg.norm(op1.affine_matrix-op2.affine_matrix) < tol:
			stab.append(ops[i])
	return stab

#map a Wyckoff set W1 to a supergroup Wyckoff set W2
def wyckoff_set_map(W1, W2):
	if len(W2)<len(W1): print("Error: cannot map to a smaller set")
	mapping = []
	for i in range(len(W1)):
		for j in range(len(W2)):
			if W1[i] == W2[j]:
				mapingp.append(j)
				break
		if mapping[i] == []:
			print("Error: could find map.")
			break
	return mapping

#Brute-force check whether or not an operation is pseudosmmetric
#with respect to a structure hstruct. Return a displacement map from point i to
#its closest transformed neighbor j.
def naive_check(hstruct, op):
	n = len(hstruct.frac_coords)
	hatoms = []
	gatoms = []
	species = []
	mapping = []
	for i in range(n):
		hatoms.append(hstruct.sites[i])
		gatoms.append(pymatgen.core.sites.PeriodicSite(hatoms[i].specie, myop.operate(hatoms[i].frac_coords), hatoms[i].lattice))
		species.append(hstruct.species[i])
	print("Shortest distances from point i in h to point j after transformation:")
	max_min_value = 0
	for i in range(n):
		min_index = i
		min_value = hatoms[i].distance(gatoms[i])
		for j in range(n):
			#print(hatoms[i].distance(gatoms[j])/2)
			if species[i] == species[j] and (hatoms[i].distance(gatoms[j]) < min_value):
				min_index = j
				min_value = hatoms[i].distance(gatoms[j])
		mapping.append(min_index)
		if min_value > max_min_value:
			max_min_value = min_value
		print("--------i="+str(i)+" ---->---- j= "+str(min_index)+"---------")
		print(min_value/2)
	if len(np.unique(mapping)) == n:
		print("All points uniquely mapped.")
	else: print(str(n-len(np.unique(mapping)))+" points not mapped to.")
	print("Largest displacement: "+str(max_min_value))
	return mapping

def get_wyckoff_positions(hall_number):
	'''Return wyckoff positions for a group as a 2d array of SymmOp objects.
	1st index is position (Wyckoff letter, 0 = general position), 2nd index
	gives elements within a Wyckoff position. Note: face positions (F/A/B/C
	groups) ARE included.'''
	array = []
	wyckoff_positions = make_sitesym.get_wyckoff_position_operators('database/Wyckoff.csv', hall_number)
	for x in wyckoff_positions:
		temp = []
		for y in x:
			temp.append(pymatgen.core.operations.SymmOp.from_rotation_and_translation(list(y[0]), y[1]/24))
		array.append(temp)
	return array

def wyckoff_split(WG, WH, position):
	'''Function to split a Wyckoff position into positions of a subgroup.
	WG and WH are sets of Wyckoff positions for the supergroup G and subgroup H,
	respectively. position is the desired position in G to split: 0 represents
	the general position.'''
	wg = WG[position]
	gsymm = [] #H site symmetry for elements of wg
	for gel in wg:
		gsymm.append(op_stabilizer(gel, WH[0]))
	hsymm_array = [] #array of site symmetries for elements of WH
	for wh in WH:
		hsymm = []
		for hel in wh:
			hsymm.append(op_stabilizer(hel, WH[0]))
		hsymm_array.append(hsymm)
	mapping = []
	for gel in wg:
		For wh in WH:#Loop over wyckpos's wh in WH, start with largest mult.
			temp2 = temp1
			if wh.mult <= len(temp1): \
			for j in range(wh.mult): #Loop over elements in wh
				for k in range(len(wg.pos)): #loop over elements in wg

def wyckoff_split_from_hall_number(hall1, hall2, position):
	pos1 = get_wyckoff_positions(hall1)
	pos2 = get_wyckoff_positions(hall2)
	return wyckoff_split(pos1, pos2, position)
#Main program
#-------------------------------------------
#import our base structure from a cif file
#struct1 will be stored as a pymatgen.core.structure Structure class object
'''path = input("Relative CIF file path: ")
mystruct1 = pymatgen.core.structure.Structure.from_file(path, primitive=False, sort=False, merge_tol=0.01)
sga = pymatgen.symmetry.analyzer.SpacegroupAnalyzer(mystruct1)'''


#Find the supergroups of struct1
#Get the transformation W,w from G to H
#Find a/the generator g from G to G

#test functionality
#-----------------------------------------------------------------------
#print("Input structure:")
#print(mystruct1)
#Get a list of symmetry operations for the structure

text = ['x,y,z', '-y,x-y,z', '-x+y,-x,z', '-x,-y,z+1/2',
		'y,-x+y,z+1/2', 'x-y,x,z+1/2', 'y,x,-z+1/2', 'x-y,-y,-z+1/2',
		'-x,-x+y,-z+1/2', '-y,-x,-z', '-x+y,y,-z', 'x,x-y,-z',
		'-x,-y,-z', 'y,-x+y,-z', 'x-y,x,-z', 'x,y,-z+1/2',
		'-y,x-y,-z+1/2', '-x+y,-x,-z+1/2', '-y,-x,z+1/2', '-x+y,y,z+1/2',
		'x,x-y,z+1/2', 'y,x,z', 'x-y,-y,z', '-x,-x+y,z']
group193 = []
for x in text:
	group193.append(pymatgen.core.operations.SymmOp.from_xyz_string(x))

text = [
'x,y,z', 	'-x,-y,z', 	'-x,y,-z', 	'x,-y,-z',
'z,x,y', 	'z,-x,-y', 	'-z,-x,y', 	'-z,x,-y',
'y,z,x', 	'-y,z,-x', 	'y,-z,-x', 	'-y,-z,x',
'y,x,-z', 	'-y,-x,-z',	'y,-x,z', 	'-y,x,z',
'x,z,-y', 	'-x,z,y', 	'-x,-z,-y', 'x,-z,y',
'z,y,-x', 	'z,-y,x', 	'-z,y,x', 	'-z,-y,-x',
'-x,-y,-z',	'x,y,-z', 	'x,-y,z', 	'-x,y,z',
'-z,-x,-y',	'-z,x,y', 	'z,x,-y', 	'z,-x,y',
'-y,-z,-x',	'y,-z,x', 	'-y,z,x', 	'y,z,-x',
'-y,-x,z', 	'y,x,z', 	'-y,x,-z', 	'y,-x,-z',
'-x,-z,y', 	'x,-z,-y', 	'x,z,y', 	'-x,z,-y',
'-z,-y,x', 	'-z,y,-x', 	'z,-y,-x', 	'z,y,x']
group221 = []
for x in text:
	group221.append(pymatgen.core.operations.SymmOp.from_xyz_string(x))

text = [
'x,y,z', 	'-x,-y,z', 	'-x,y,-z', 	'x,-y,-z',
'-x,-y,-z', 	'x,y,-z', 	'x,-y,z', 	'-x,y,z']
group47 = []
for x in text:
	group47.append(pymatgen.core.operations.SymmOp.from_xyz_string(x))


#H=185, G=186 (Water), index = 2
#myop = pymatgen.core.operations.SymmOp.from_xyz_string('x, y, z+1/2')
#H=185, G=186 (Water), index = 3
#myop = pymatgen.core.operations.SymmOp.from_xyz_string('x+1/3, y+2/3, z')
#H=185, G=193 (Water), index = 2
#myop = pymatgen.core.operations.SymmOp.from_xyz_string('-x, -x+y, 1/2-z')
#H=38, G=65 (BaTiO3)
myop = pymatgen.core.operations.SymmOp.from_xyz_string('x, y, -z')
#myop = pymatgen.core.operations.SymmOp.from_xyz_string('x, y, -z')

#generators from P-type lattice to F-type lattice
myop1 = pymatgen.core.operations.SymmOp.from_xyz_string('x, y+1/2, z+1/2')
myop2 = pymatgen.core.operations.SymmOp.from_xyz_string('x+1/2, y, z+1/2')
myop3 = pymatgen.core.operations.SymmOp.from_xyz_string('x+1/2, y+1/2, z')

group225 = []
for x in group221:
	group225.append(x)
	group225.append(compose_ops(myop1,x))
	group225.append(compose_ops(myop2,x))
	group225.append(compose_ops(myop3,x))

start = timer()
#--------------Timer start
wyckoff_split_from_hall_number(185, 200, 3)
#--------------Timer stop
end = timer()
print("Time elapsed: "+str(end-start))



