import argparse
import sys

parser = argparse.ArgumentParser(description='more nuanced breakdown of metaxa stats.\n')
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

required.add_argument("--taxlevel", type=int, default=2, help="Pick a taxonomic level from 1 (kingdom) to 7 (species). Must be an integer. Defaults to 2 (Phylum). Domain and Kingdom always reported, so both 0 and 1 turn this off.")
optional.add_argument("--filter", type=float, default=0.1, help="If a particular taxonomic classification makes up less than a given percent of the total reads, it doesn't show up in the results.\
  A decimal value (<1) is recommended, as this option is primarily designed to remove misclassified reads that clutter the output. Default = 0.2, set to 0 if you want the filter off entirely.")
required.add_argument("-i", "--input", help="path to input file")
optional.add_argument("-o", "--output", default="taxstats.txt", help="output name/location\nDefault is taxstats.txt")
required.add_argument("--rna", default="lsu", help="rRNA type; ssu or lsu")

args = parser.parse_args()

infile = args.input
outfile = args.output
level = args.taxlevel
filter = args.filter
rna = args.rna


try:
		L = int(level)
except:
		print("You must input a tax level! See help info ('metaxa2stats.py -h') for details")
		sys.exit(0)
if level not in range (0,8):
		print("Tax level must be between 0 and 7!")
		sys.exit(0)

if rna != "lsu" and rna != "ssu":
		print('rRNA type must be "ssu" or "lsu"!')
		sys.exit(0)
rna=rna.upper()


if level == 2:
	tclass = "Phylum"
if level == 3:
	tclass = "Class"
if level == 4:
	tclass = "Order"
if level == 5:
	tclass = "Family"
if level == 6:
	tclass = "Genus"
if level == 7:
	tclass = "Species"


inf = open(infile, "r")
total = 0
categories = dict()

for line in inf:
	total+=1

	###
	#Break down line by whitespace, then only look at the 1st and SOMETIMES 2nd index to get an array with taxonomic info
        #
        #This is because lines normally look like this:
        #MISEQ06:640:000000000-AWW91:1:2114:22959:3582   Chloroplast;;;;;        100     150     100
        #With all of the taxonomic information in the second column, with levels broken down by ";"
        #But sometimes they put spaces within the taxonomic information like this:
        #MISEQ06:640:000000000-AWW91:1:2108:13127:26069  Chloroplast;Unknown Eukaryote   N/A      N/A    N/A
        #So I check to see if the third column is only letters, and if it is I merge it with the second before splitting by ";"
	###

	summary=line.split()
	taxa=summary[1].split(';')
	if summary[2].isalpha() == True:
		fulltaxa=summary[1]+summary[2]
		taxa=fulltaxa.split(';')
	else:
		taxa=summary[1].split(';')

	###
	#Label the read's domain, kingdom, and input taxlevel for counting.
        #
        #If encounters just "Mitochondria", list will be ['Mitochondria'], with no other classification
        #In these cases I want to set the other taxonomic levels to "Unknown".
        #Hence the try/except statement: If I look for an index that doesn't exist, write that level as "Unknown".
        #
        #However, some are entered as "Mitochondria;;;;", resulting in a list of ['Mitochondia','','','','']
        #I also want to set the other tax levels to "Unknown" in this case, but the indexes DO exist here, they're just empty
        #Hence the if/else in the try statement: if the index exists but the string inside it is empty, write as "Unknown"
        #
        #Dots added so that we can concatenate levels to count 
        #(want to count "Mitochondria.Unknown" and "Chloroplast.Unknown" separately)
        #And then have something to split by when returning/printing later
        #(I.e. be able to easily turn "Mitochondria.Unknown" back into "Unknown" when reporting)
	###

	domain=taxa[0]
	try:
		if not taxa[1]:
			kingdom=".Unknown"
		else:
			kingdom="."+taxa[1]
	except:
		kingdom=".Unknown"

	try:
		if level not in range (2,8):
			pass
		elif not taxa[level]:
			family=".Unknown"
		else:
			family="."+taxa[level]
	except:
		family=".Unknown"

	###
	#Counting mechanism: dictionary.
        #Again, dots used to count the "Unknown"s separately, while being able to split back up later
	###

	if (domain+"..") not in categories: 
		categories[domain+".."]=[domain," "," ",1]
	else:
		categories[domain+".."][3]+=1
	if (domain+kingdom+".") not in categories:
		categories[domain+kingdom+"."]=[domain,kingdom," ",1]
	elif (domain+kingdom+".") != (domain+".."):
		categories[domain+kingdom+"."][3]+=1
	if level not in range(2,8):
		continue
	if (domain+kingdom+family) not in categories:
		categories[domain+kingdom+family]=[domain,kingdom,family,1]
	elif (domain+kingdom+family) != (domain+"..") and (domain+kingdom+family) != (domain+kingdom+"."):
		categories[domain+kingdom+family][3]+=1


outf = open(outfile, "w+")
if level not in range (2,8):
	outf.write("\nNo indent: Domain\n\tFirst indent: Kingdom\n\n")
else:
	outf.write("\nNo indent: Domain\n\tFirst indent: Kingdom\n\t\tSecond indent: " + tclass + "\n\n")
outf.write("Total " + rna + " rRNA: " + str(total) + "\n")

###Sort first by domain alphabetically, then by kingdom alphabetically, then by input tax level by frequency from highest to lowest.

sorted_categories=sorted(categories.items(), key=lambda x: (x[1][0],x[1][1],-x[1][3]))

###
#Filter used to remove entries with only a handful of (probably accidental) counts
#Split by the dots from earlier, then check to see if the entry is classified just by domain, to kingdom, or up to the input level.
#If input level, tabbed twice. If kingdom, tabbed once. If domain, not tabbed. Gives clear organization in results.
###

for category in sorted_categories:
	name=category[0]
	number=categories[name][3]
	percentage = round(float(number*100/total),2)
	percentagestr = str(percentage)+"%"
	levelfilter = filter * level

	if categories[name][2] != " " and level in range (2,8):
		if "Unknown" in categories[name][1] and "Unknown" in categories[name][2]:
			continue
		if percentage < levelfilter:
			continue
		else:
			outf.write("\t\t%s: %s (%s)\n" % (name.split(".")[2], number, percentagestr))
	elif categories[name][1] != " ":
		if percentage < filter:
			continue
		else:
			outf.write("\t%s: %s (%s)\n" % (name.split(".")[1], number, percentagestr))
	else:
		outf.write("%s: %s (%s)\n" % (name.split(".")[0], number, percentagestr))
inf.close
outf.close

print("Done!")

