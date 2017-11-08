import pymatgen
import make_sitesym as ss

wyckoff_positions = ss.get_wyckoff_position_operators ('database/Wyckoff.csv', 423)

print("Wyckoff positions of Hall Group 423:")
for x in wyckoff_positions:
	for y in x:
		op = pymatgen.core.operations.SymmOp.from_rotation_and_translation(list(y[0]), y[1]/24)
		print(op.as_xyz_string())
