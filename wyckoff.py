'''Code to calculate the Wyckoff Splittings for a group-supergroup pair
Given: group H, minimal supergroup G
Given: Symmetry operations for G, H
Given: Special Wyckoff positions for G, H
	Note: we need both to be given in the H basis
	Convert from G basis to H basis by applying g^-1 to Wg
	Also apply g to x,y,z (diagonalize)
Need: Connection between symmopps for each wyckpos.'''
#------------------------------------------------------------------------------
#Class for storing Wyckoff positions
#Code incomplete
class wyck_pos:
	def from_text(text):
		wpos = wyck_pos
		#import wyckoff position from text (Wyckoff.csv)
		#Issue: +(0,1/2,1/2) type positions are implied by spacegroup letter
		#ex: 122'F 2 2 2' position k has mult 16, but only 4 pos's are listed
		#wpos.mult =  #multiplicity of position
		#wpos.letter =  #Wyckoff letter
		#wpos.elem = [] #set of specific positions in orbit, SymmOp objects
		#wpos.symm_set = [SymmOp] #stab() subset of general positions, + (x,y,z)
		return wpos
	def from_hm_symbol(symbol):
		wpos = wyck_pos
		return wpos
	def from_hm_number(number):
		wpos = wyck_pos
		return wpos
	def from_hall_symbol(symbol):
		wpos = wyck_pos
		return wpos
	def from_hall_number(number):
		wpos = wyck_pos
		return wpos
	def transformed(transformation):
		#Transform a wyckoff position in G to the setting of H
		#Consider variable name changes (see WYCKSPLIT examples)
		#Consider periodicity constraints (no negative translations)
		for x in self.pos: pass
		for x in self.symm_set: pass

'''class xyz_map:
	def __init__:
		for i in range(3):
			self[i].vector = None
			self[i].scalar = None'''

class wyck_split(wg, WH, index):
	#class to calculate and store Wyckoff splitting from wg to H
	#It is assumed wg and Wh are in the same setting. Call wg.transform() first
	#wg is a single wyck_pos in G; Wh is the array of all wyck_pos's in H
	#index is the size of the general position of G divided by that of H
	def __init__:
		indices = []
		xyzmap = [[]]
		temp1 = wg.pos
		For wh in WH:#Loop over wyckpos's wh in WH, start with largest mult.
			temp2 = temp1
			if wh.mult <= len(temp1): \
			for j in range(wh.mult): #Loop over elements in wh
				for k in range(len(wg.pos)): #loop over elements in wg
					for k in range(3): #Loop over x,y,z coordinates
						a = wh.pos[j].affine_matrix
						hvector = [a[l][0], a[l][1], a[l][2]]
						hscalar = a[l][3]
						a = wg.pos[k].affine_matrix
						gvector = [a[l][0], a[l][1], a[l][2]]
						gscalar = a[l][3]
						if vector == [0,0,0]: #if h element is a point
						
						else: #if h element is not a point
'''We can check Wh: if linked to n Wg elements,
	then we need n atoms of Wh type'''

