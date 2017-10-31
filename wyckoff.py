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
		#wpos.pos = [] #set of specific positions in orbit
		#wpos.symm_set = [SymmOp] #set of general positions not in orbit, + (x,y,z)
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
	def transform(self, transformation):
		#Transform a wyckoff position in G to the setting of H
		#Consider variable name changes (see WYCKSPLIT examples)
		#Consider periodicity constraints (no negative translations)
		for x in self.pos: pass
		for x in self.symm_set: pass
	def split(self, wh, transformation)):
	#split a Wyckoff position in G into a set in H, given a G->H transformation
		index = len(self.pos)/len(wh.pos)		
		mapping = []
		n_remaining = self.mult
		temp_pos = self.pos
		wg.transform(transformation)
		For x in wh:
			if x.mult =< n:
				#check site against wg
				#check if all pos's could be in wg (with def of x,y,z)
				xyzmap = ['','','']
				for y in wh.pos:
					check if mapping exists for wg.pos->
		return mapping

wg = [] #set of Wyckoff positions in G
wh = [] #set of Wyckoff positions in subgroup H
#0th element in wg or wh is the 'a' orbit, to make looping easier
For i in range(len(Wg)): 
	wg[i].split(wh)


'''We can check Wh: if linked to n Wg elements,
	then we need n atoms of Wh type'''

