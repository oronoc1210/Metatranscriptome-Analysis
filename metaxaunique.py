import sys
import argparse
import operator
import numpy as np

parser = argparse.ArgumentParser(description='more nuanced breakdown of metaxa stats.\n')
parser.add_argument("--input1", help="path to first input file")
parser.add_argument("--input2", help="path to second input file")
args = parser.parse_args()

in1 = open(args.input1, "r")
in2 = open(args.input2, "r")
sort1 = []
sort2 = []


lines1 = []
lines2 = []
for line in in1:
	id = line.split()[0]
	idparts = id.split(':')
	lines1.append(idparts)
for line in in2:
	id = line.split()[0]
	idparts = id.split(':')
	lines2.append(idparts)
	
sort1 = sorted(lines1, key=operator.itemgetter(4,5,6))
sort2 = sorted(lines2, key=operator.itemgetter(4,5,6))

in1.close()
lines1 = None
in2.close()
lines2 = None


t2 = len(sort2)
start = round(t2 / 2)
moves = 0
matches = []


def binary(line1, list, n):
	first = 0
	last = len(list) - 1
	found = False
	while first <=last and not found:
		midpoint = (first + last) // 2
		#print("midpoint: " + str(midpoint) + "\nfirst: " + str(first) + "\nlast: " + str(last) + "\n")
		if list[midpoint][n] == line1[n]:
			found = True
		else:
			if line1 < list[midpoint]:
				last = midpoint - 1
			else:
				first = midpoint + 1
	if not found:
		return
	if found and n == 4:
		binary(line1, list, 5)
	if found and n == 5:
		binary(line1, list, 6)
	if found and n == 6:
		matches.append(line1)
	return found
		
for value in sort1:
	binary(value, sort2, 4)
	

pct1 = round(len(matches) / len(sort1) * 100, 3)
pct2 = round(len(matches) / len(sort2) * 100, 3)
print("total in1: " + str(len(sort1)) + "\ntotal in2: " + str(len(sort2)) + "\n")
print("# matches: " + str(len(matches)) + "\npercent in1: " + str(pct1) + "%\npercent in2: " + str(pct2) + "%\n")
print("some examples: \n")
for i in range (0, 20):
	print(":".join(matches[i]))
sys.exit(0)
