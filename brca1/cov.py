#!/usr/bin/env python

o = open('brca_depth.txt', 'w')

total = 0

with open('brca1.final.bed') as f:
	for line in f:
		array = line.split()
		chr = array[0]
		start = int(array[1])
		end = int(array[2])
		depth = array[4]
		total = total + (end - start)

		while (start < end):
			o.write(chr + "\t" + str(start) + "\t" + depth)
			o.write("\n")
			start = start + 1

print total
f.close()
o.close()
