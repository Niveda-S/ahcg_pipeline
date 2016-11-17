# argv[1] is the pathogenic annotation of variants call
# argv[2] is the variants read depth file
import sys

def print_header():
	print("CHROM" + "\t" + "POSITION" + "\t" + "REF" + "\t" + "ALT" + "\t" + "PATHOGENIC")

pathogenics = []
depths = []
with open(sys.argv[1]) as fin_pathogenic:
	pathogenics = fin_pathogenic.readlines()
with open(sys.argv[2]) as fin_depth:
	depths = fin_depth.readlines()

print_header()

for pathogenic in pathogenics:
	find = 0
	for depth in depths:
		if pathogenic.split(",")[1][5:].strip() == depth.split()[1].strip():						
			find = 1
			depth_number = depth.split()[2].strip()

	line = pathogenic.split(",")[0][-5:] + "\t" + pathogenic.split(",")[1][5:] + "\t" + pathogenic.split(",")[2][5:] + "\t" + pathogenic.split(",")[3][6] + "\t" + pathogenic.split(")")[1].strip()
	if find == 1:
		print(line.strip() + "\t" + depth_number)
	else:
		print(line.strip() + "\t" + "NA")
