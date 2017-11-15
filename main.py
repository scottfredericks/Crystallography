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
	print("==========Time elapsed: "+str(end-start)+"==========")
	hsymm_array = [] #array of site symmetries for elements of WH
	for wh in WH:
		hsymm = []
		for hel in wh:
			hsymm.append(op_stabilizer(hel, hops))
		hsymm_array.append(hsymm)
	mapping = []
	T, t = list([None, None, None]), list([None, None, None])
	counter = 0
	for i in range(len(wg)):
		gel = wg[i]
		X, x = gel.rotation_matrix, filter_site(gel.translation_vector)
		found = False
		for j in range(len(WH)):
			wh = WH[j]
			if len(gsymm[i]) == len(hsymm_array[j][0]):
				for k in range(len(wh)):
					if gsymm[i] == hsymm_array[j][k]:
						#Try finding a map
						#loop over x,y,z
						hel = wh[k]
						Y, y = hel.rotation_matrix, filter_site(hel.translation_vector)
						canmap = True
						for l in range(3):
							if list(X[l]) == list([0,0,0]):
								if x[l] != y[l]:
									canmap = False
									break
						if canmap:
							found = True
							mapping.append([j])
					if found: break
			if found: break
		if found == False:
			print("Error: Could not find wyckoff mapping.")
			print("i: "+str(i)+", j: "+str(j)+", k: "+str(k))
			mapping.append("Error")
	
	#Using position-position mapping, find exact transformations
	known_maps = []
	for i in range(len(wg)):
		gel = wg[i]
		wh = WH[mapping[i][0]]
		if known_maps == []:
			T = []
			t = []
			Y, y = gel.rotation_matrix, filter_site(hel.translation_vector)
			X, x = wh[0].rotation_matrix, filter_site(wh[0].translation_vector)
			for j in range(3):
				t.append(x[j]-y[j])
				if X[j][j] == 1:
					T.append(Y[j])
				elif X[j] = [0,0,0]:
					T.append([0,0,0])
				else:
					for k in range(3):
						if X[j][k] == 1:
							
					 
			trans = pymatgen.core.operations.SymmOp.from_rotation_and_translation(T,t)
			known_maps.append(trans)
			mapping[i].append(0)
			mapping[i].append(trans)
		else:
			found = False
			for j in range(len(known_maps)):
				op =  known_maps[j]
				for k in range(len(wh)):
					hel = wh[k]
					generated = compose_ops(op, hel)
					if np.linalg.norm(generated.affine_matrix-gel.affine_matrix) < .001:
						found = true
						
						known_maps.append(op)
						mapping[i].append(k)
						mapping[i].append(op)
					if found: break
				if found: break
			if not found: pass
	print(wg[0].as_xyz_string())
	print(wh[0].as_xyz_string())
	
	#return letters instead of integers
	for i in range(len(mapping)):
		mapping[i][0] = letters[len(WH)-1-mapping[i][0]]
	return mapping

def wyckoff_split_from_hall_number(hall1, hall2, position):
	pos1 = get_wyckoff_positions(hall1)
	pos2 = get_wyckoff_positions(hall2)
	#Need to perform transformation from group 1->2
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

#myop = pymatgen.core.operations.SymmOp.from_xyz_string('x, y, -z')

start = timer()
print("===================Timer started===================")
#--------------Timer start
#H-M groups 225(Fm-3m)=523 (supergroup), and 221(Pm-3m)=517
print(wyckoff_split_from_hall_number(523, 517, 'h'))

#--------------Timer stop
end = timer()
print("==========Time elapsed: "+str(end-start)+"==========")



