'''Code to calculate the Wyckoff Splittings for a group-supergroup pair
Given: group H, minimal supergroup G
Given: Symmetry operations for G, H
Given: Wyckoff positions for G, H
	Note: we need both to be given in the H basis
	Convert from G basis to H basis by applying g^-1 to Wg
	Also apply g to x,y,z (diagonalize)
Need: Connection between symmopps for each wyckpos.


For a given classification of Wyckoff position type in G:
Loop over all Wyckoff positions in G (Wg)
	Set Wg in basis of Wh
		#can be simplified if we always have coset representative
		multiply by inverse of transformation from H to G
		Check for negative numbers, add 1
		rename x,y,z if in wrong positions
	Attribute one Wh to each Wg
	Make a cross-list linking Wg to Wh
Gives list splitting(Wg), Wh
We can check Wh: if linked to n Wg elements,
	then we need n atoms of Wh type


'''

