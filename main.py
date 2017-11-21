'''
Scott Fredericks, 2017
UNLV department of Physics
Advisor: Qiang Zhu
-----------------------------------
Code for finding pseudosymmetry in crystal structures. Input files may be CIF
or POSCAR. In general, "G" represents the supergroup/potential pseudosymmetry
group, while "H" represents the subgroup/the space group of the provided
structure. The algorithm for finding pseudosymmetry is as follows:
1. Obtain the spacegroup (H) information for the provided structure (from file)
2. Determine the chain of minimal supergroups G, with transformation matrices
3. Find the splitting of the Wyckoff positions from G to H
4. If the Wyckoff splitting meets certain criteria, map each atomic species
	into a new Wyckoff position in G (multiple mappings may be possible)
5. Determine the largest displacement produced by the mapping. If the
	displacement is smaller than a given tolerance, the supergroup is flagged
	as pseudosymmetric
-----------------------------------
Dependencies:
pymatgen (http://pymatgen.org/)
numpy (http://www.numpy.org/)
pyspglib (https://atztogo.github.io/spglib/python-spglib.html#python-spglib)
'''
import sys
import pymatgen
import pymatgen.symmetry.groups
import pymatgen.symmetry.analyzer
import pymatgen.symmetry.settings
import numpy as np
from timeit import default_timer as timer
import database.make_sitesym as make_sitesym

#Constants
#--------------------------
letters = 'abcdefghijklmnopqrstuvwxyz'

#Function and class definitions
#--------------------------
def filter_site(v):
	w = v
	for i in range(len(w)):
		while w[i]<0: w[i] += 1
		while w[i]>=1: w[i] -= 1
	return w

def compose_ops(op1,op2):
	#Apply transform op1 to op2, then normalize
	R = np.dot(op1.rotation_matrix,op2.rotation_matrix)
	t = op1.operate(op2.translation_vector)
	t = filter_site(t)
	new = pymatgen.core.operations.SymmOp.from_rotation_and_translation(R,t)
	return new

def site_stabilizer(site1, ops, tol=.01):
	#Return the subset of ops which leaves the PeriodicSite site1 invariant
	stab = []
	for i in range(len(ops)):
		site2 = pymatgen.core.sites.PeriodicSite(site1.specie, ops[i].operate(site1.frac_coords), site1.lattice)
		if site1.distance(site2) <= tol:
			stab.append(ops[i])
	return stab

def op_stabilizer(op1, ops, tol=.001):
	#return the subset of ops which leaves a SymmOp object invariant
	#used for finding site symmetry of special Wyckoff position elements
	stab = []
	for i in range(len(ops)):
		op2 = compose_ops(ops[i], op1)
		if np.linalg.norm(op1.affine_matrix-op2.affine_matrix) < tol:
			stab.append(ops[i])
	return stab

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
			temp.append(pymatgen.core.operations.SymmOp.from_rotation_and_translation(list(y[0]), filter_site(y[1]/24)))
		array.append(temp)
	return array

def wyckoff_split(WG, WH, letter):
	'''Function to split a Wyckoff position into positions of a subgroup.
	WG and WH are sets of Wyckoff positions for the supergroup G and subgroup H,
	respectively. letter is the position in G you wish to split.'''
	gops = WG[0]
	hops = WH[0] #General site-symmetries of H
	index = len(gops)/len(hops)
	position = len(WG) - 1 - letters.find(letter)
	wg = WG[position]
	gsymm = [] #H site symmetry for elements of wg with respect to H
	for i in range(len(wg)):
		gsymm.append(op_stabilizer(wg[i], hops))
	end = timer()
	hsymm_array = [] #array of site symmetries for elements of WH
	for wh in WH:
		hsymm = []
		for hel in wh:
			hsymm.append(op_stabilizer(hel, hops))
		hsymm_array.append(hsymm)
	mapping = []
	known_maps = []
	for i in range(len(wg)):
		#Loop over elements in the G Wyckoff Position
		gel = pymatgen.core.operations.SymmOp.from_rotation_and_translation(wg[i].rotation_matrix, filter_site(wg[i].translation_vector))
		found = False
		for j in range(len(known_maps)):
			op =  known_maps[j]
			for l in range(len(WH)):
				wh = WH[l]
				if len(gsymm[i]) == len(hsymm_array[l][0]):
					for k in range(len(wh)):
						hel = pymatgen.core.operations.SymmOp.from_rotation_and_translation(wh[k].rotation_matrix,filter_site(wh[k].translation_vector))
						generated = compose_ops(hel, op)
						if np.linalg.norm(generated.affine_matrix-gel.affine_matrix) < .01:
							found = True
							mapping.append([l])
							mapping[i].append(k)
							mapping[i].append(op)
						if found: break
			if found: break
		if not found:
			for l in range(len(WH)):
				#Loop over Wyckoff positions in H
				wh = WH[l]
				if len(gsymm[i]) == len(hsymm_array[l][0]):
					B, b = gel.rotation_matrix, filter_site(gel.translation_vector)
					for k in range(len(wh)):
						hel = wh[k]
						A, a = hel.rotation_matrix, filter_site(hel.translation_vector)
						t = []
						T = []
						for j in range(3):
							if A[j][j] == 1:
								T.append(B[j])
								t.append(b[j]-a[j])
							elif A[j][j] == -1:
								T.append(-1*B[j])
								t.append(-(b[j]-a[j]))
							else:
								t.append(0)
								if j == 0:
									T.append([1,0,0])
								elif j == 1:
									T.append([0,1,0])
								elif j == 2:
									T.append([0,0,1])
						trans = pymatgen.core.operations.SymmOp.from_rotation_and_translation(T,t)
						generated = compose_ops(hel, trans)
						if np.linalg.norm(generated.affine_matrix-gel.affine_matrix) < .01:
							found = True
							known_maps.append(trans)
							mapping.append([l])
							mapping[i].append(k)
							mapping[i].append(trans)
							break
		if not found:
			print("Error: mapping not found for i = "+str(i))
			print(gel.as_xyz_string())
	#return letters instead of integers
	for i in range(len(mapping)):
		mapping[i][0] = letters[len(WH)-1-mapping[i][0]]
	return mapping

def wyckoff_split_from_hall_number(hall1, hall2, position):
	pos1 = get_wyckoff_positions(hall1)
	pos2 = get_wyckoff_positions(hall2)
	#TODO find and perform transformation from group 1->2
	return wyckoff_split(pos1, pos2, position)
#Main program
#-------------------------------------------
#import our base structure from a cif file
#struct1 will be stored as a pymatgen.core.structure Structure class object
'''path = input("Relative CIF file path: ")
mystruct1 = pymatgen.core.structure.Structure.from_file(path, primitive=False, sort=False, merge_tol=0.01)
sga1 = pymatgen.symmetry.analyzer.SpacegroupAnalyzer(mystruct1)
data1 = sga1.get_symmetry_dataset()
h_m_number1 = data1['number']
sg1 = pymatgen.symmetry.groups.SpaceGroup.from_int_number(h_m_number1)'''

#TODO
'''Find minimal supergroup list. Loop over this list.
For each supergroup, find transformation matrix and origin shift if necessary.
For each structure, perform Wyckoff split, then map from H to G. Display results.'''

#test functionality
#-----------------------------------------------------------------------
'''print("Input structure:")
print(mystruct1)'''

start = timer()
print("===================Timer started===================")
#--------------Timer start
#H-M groups 225(Fm-3m)=523 (supergroup), and 221(Pm-3m)=517
supergroup, group = 523, 517
print("Splitting of Wyckoff position from supergroup "+str(supergroup)+" into group "+str(group))
letter = input("letter (in supergroup): ")
if letter not in letters:
	print("Not a valid letter.")
	sys.exit()
array = wyckoff_split_from_hall_number(supergroup, group, letter)
for x in array:
	print(x[0]+", "+str(x[1])+", ("+x[2].as_xyz_string()+")")
print("(Letter, element # in subgroup Wyckoff position, mapping)")
'''Useful functions
pymatgen.symmetry.groups.SpaceGroup.from_int_number(#)
pymatgen.symmetry.groups.get_symm_data()['maximal_subgroups']
sga.get_symmetry_dataset()['hall_number']'''
#--------------Timer stop
end = timer()
print("==========Time elapsed: "+str(end-start)+"==========")



