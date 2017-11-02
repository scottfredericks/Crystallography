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
		#wpos.pos = [] #set of specific positions in orbit, SymmOp objects
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

class wyck_split(wg, wh):
	#class to calculate and store Wyckoff splitting from G to H
	#It is assumed wg and wh are in the same setting. Call wg.transform() first
	def __init__:
		indices = [[]]
		xyzmap = [[[]]]
		index = len(wg.pos)/len(wh.pos)
		temp_pos = wg.pos
		for i in range(len(wg)): #Loop over wyckpos's in G
			remaining = wg[i].mult
			For j in range(len(wh)):#Loop over wyckpos's in H
				if wh[j].mult =< remaining:
					#check site symmetry against wg; break if mismatch
					
					xyzmap = ['','',''] #mapping needs to be unique for each j
					for k in range(len(wh[j].pos)):
						#check if xyz mapping exists for wg.pos
						#if xyz mapping is not unique for some y, abort.
					#if unique xyz map exists, flag x as possible subset

'''We can check Wh: if linked to n Wg elements,
	then we need n atoms of Wh type'''

