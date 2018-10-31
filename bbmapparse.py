import argparse
import sys

parser = argparse.ArgumentParser(description="Takes bbmap results and reports taxonomic classifications. Takes queries of things you want to be in the string as well as things you don't want to be in the string. Summarizes stats.")

parser.add_argument("-i", "--input", help="path to input file")
parser.add_argument("--test", help="string to test for in reads; if in string, it's reported; separate multiple queries with a single comma ','")
parser.add_argument("--untest", help="string to test against in reads; if in string, don't report even if --test arg is also in string ; separate multiple queries with a single comma ','")
parser.add_argument("--nolsu", default="False", help="Ignore LSU results. Default false.")
parser.add_argument("--nossu", default="False", help="Ignore SSU results. Default False.")
parser.add_argument("--reportnum", default=50, help="Number of reads to report before stopping. Default 50.")


args = parser.parse_args()
infile = args.input
teststr = args.test
try:
	teststrings = teststr.split(",")
except:
	teststrings = False
unteststr = args.untest
try:
	unteststrings = unteststr.split(",")
except:
	unteststrings = False
if args.nolsu.lower() == "true" or args.nolsu.lower() == "t":
	ilsu = True
else:
	ilsu = False
if args.nossu.lower() == "true" or args.nossu.lower() == "t":
	issu = True
else:
	issu = False
teststop = int(args.reportnum)

inf = open(infile, "r")
total = 0
rnatotal = 0
lsutotal = 0
ssutotal = 0
testtot = 0

print('Search string: "%s" ; Unsearch string: "%s" ; Ignore LSU: "%s" ; Ignore SSU: "%s"' % (teststrings, unteststrings, str(ilsu), str(issu)))
for line in inf:
	report = True
	if testtot >= teststop:
		report = False
	if line.startswith( "@" ) or line == "\n":
		continue
	else:
		total += 1
		try:
			region = line.split("\t")[2]
		except:
			print(line)
			exit(1)
		if ";" not in region:
			continue
		rnatotal += 1
		if "LSU" in region:
			if ilsu == True:
				continue
			rna = "LSU"
			lsutotal += 1
		elif "SSU" in region:
			if issu == True:
				continue
			rna = "SSU"
			ssutotal += 1
		parts = region.split()
		taxa = " ".join(parts[1:])

		test = True
		untest = True

		if teststrings == False:
			pass
		else:
			for item in teststrings:
				if item not in taxa:   # If ANY of the test strings aren't there, test fails.
					test = False
		if unteststrings == False:
			pass
		else:
			for item in unteststrings:
				if item in taxa:       # If ANY of the test strings are there, untest fails.
					untest = False

		if test == True and untest == True:    # If all test and no untest present, counter increments.
			if report == True:
				print(rna + "\t" + taxa)
			testtot += 1

print("Total reads: %d\nTotal rRNA: %d\nTotal LSU: %d\nTotal SSU:  %d\nTotal match: %d\n" % (total, rnatotal, lsutotal, ssutotal, testtot))
