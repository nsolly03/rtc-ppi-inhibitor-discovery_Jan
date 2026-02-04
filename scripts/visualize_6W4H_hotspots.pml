# PyMOL script to visualize 6W4H NSP10-NSP16 hot spots
# Based on Trepte et al. 2024 findings

# Load structure
load data/structures/pdb/6W4H.pdb

# Color by chain
hide everything
show cartoon

# Identify chains (usually A=NSP10, B=NSP16, but verify)
# Adjust chain IDs based on analysis script output

# Color chains
color marine, chain A     # NSP10
color salmon, chain B     # NSP16

# Show hot spot residues as spheres
# NSP10 Lys93
select nsp10_lys93, chain A and resi 93
show spheres, nsp10_lys93
color red, nsp10_lys93
label nsp10_lys93 and name CA, "NSP10 Lys93"

# NSP16 Asp106  
select nsp16_asp106, chain B and resi 106
show spheres, nsp16_asp106
color blue, nsp16_asp106
label nsp16_asp106 and name CA, "NSP16 Asp106"

# Show interface residues within 10 Å of Lys93
select interface_10A, (chain A or chain B) and (byres (chain A and resi 93) around 10)
show sticks, interface_10A
color yellow, interface_10A and chain A
color cyan, interface_10A and chain B

# Draw distance between hot spots
distance hotspot_distance, nsp10_lys93 and name CA, nsp16_asp106 and name CA
color magenta, hotspot_distance

# Center view on interface
center nsp10_lys93
zoom interface_10A

# Improve visualization
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
set sphere_scale, 0.5
set label_size, 14

# Save image
png 6W4H_hotspots.png, width=1600, height=1200, dpi=300, ray=1

# Print info
print "="*60
print "6W4H NSP10-NSP16 Visualization"
print "="*60
print "Hot spots highlighted:"
print "  - NSP10 Lys93 (RED spheres)"
print "  - NSP16 Asp106 (BLUE spheres)"
print "Interface residues (10 Å from Lys93):"
print "  - NSP10: YELLOW sticks"
print "  - NSP16: CYAN sticks"
print "Distance shown in MAGENTA"
print "="*60
