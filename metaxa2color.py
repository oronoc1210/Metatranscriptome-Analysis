import argparse
import sys
import re

from sty import fg, bg, ef, rs

parser = argparse.ArgumentParser(description='Takes both lsu and ssu metaxa runs of the same dataset and merges results. If just one of the two provided it will still report its stats.')
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')

optional.add_argument("--filter", type=float, default=1, help="If a particular taxonomic classification makes up \
\	less than this percent of the total reads, iit doesn't report that classification in the results.\
\	Exception for the Domain level, which will always be reported. \
\	Default = 1, set to 0 if you want the filter off entirely (not recommended).")
required.add_argument("--lsu", help="path to lsu input file")
required.add_argument("--ssu", help="path to ssu input file")
optional.add_argument("-o", "--output", default="taxmerge.txt", help="output name/location. Default is taxmerge.txt")

args = parser.parse_args()
lsupath = args.lsu
ssupath = args.ssu
outfile = args.output
filter = args.filter


def makecolor(percent):
    if percent < 2:
        color = 230
    elif percent < 8:
        color = 229
    elif percent < 20:
        color = 228
    elif percent < 40:
        color = 227
    elif percent < 60:
        color = 221
    elif percent < 70:
        color = 215
    elif percent < 80:
        color = 209
    elif percent < 90:
        color = 202
    else:
        color = 196
    return color

"""Takes the Metaxa2 metaxa.taxonomy.txt output and summarizes results by taxonomic level"""
def TaxParser(path, rnatype, mode):
    inf = open(path, "r")
    if mode == "w": 
        outf = open(outfile, "w+")
        parseruns = 1
    elif mode == "a":
        outf = open(outfile, "a+")
    else:
        print("mode must be 'w' or 'a'!")
        sys.exit(0)

    rna = rnatype.upper()

    total = 0
    categories = dict()

    for line in inf:
            total+=1

            summary=line.split()
            newsum = []
            for item in summary:
                    if re.search(r"[0-9/]+", item):
                            pass
                    else:
                            newsum.append(item)
            taxa = ' '.join(newsum)
            if taxa == '':
                    continue
            taxalist = taxa.split(";")

            # The species (and ONLY species) level also ENDS in ";",
            # So if a species level exists, then there will be an empty string at the end of our list.
            # Hence, this removes the last item in our list only if it's an empty string.

            if taxalist[len(taxalist)-1] == '':
                    del taxalist[len(taxalist)-1]

            for i in range(1,len(taxalist)+1):
                    fullname = '.'.join(taxalist[0:i])
                    levelname = taxalist[i-1]
                    if fullname not in categories:
                            categories[fullname] = [1, levelname]
                    else:
                            categories[fullname][0] += 1
            if taxalist[0] == "Eukaryota":
                try:
                    test = taxalist[1]
                except:
                    if "Eukaryota.Unknown Eukaryote" not in categories:
                        categories["Eukaryota.Unknown Eukaryote"] = [1, "Unknown Eukaryote"]
                    else:
                        categories["Eukaryota.Unknown Eukaryote"][0] += 1
                            

    sortedcategories=sorted(categories)

    outf.write("\nTotal " + rna + " rRNA: " + str(total) + "\n")
    for category in sortedcategories:
            tabs = category.count('.')
            number = categories[category][0]
            name = categories[category][1]
            percentage = round(number / total * 100 , 2)
            color = makecolor(percentage)
            
            if name != "" and percentage >= filter or tabs == 0:
                    outf.write(fg(color) + "    "*tabs + name + ": " + str(number) + " (" + str(percentage) + "%)\n" + fg.rs)
            
    outf.write("\n" + "-"*50 + "\n")
    if rna == "LSU":
        global lsuCategories
        global lsuTotal
        lsuCategories = categories
        lsuTotal = total
    if rna == "SSU":
        global ssuCategories
        global ssuTotal
        ssuCategories = categories
        ssuTotal = total
    inf.close()
    outf.close()
  


"""Takes the results of lsu and ssu parsing and merges them, where possible"""
def rRNA_Merge(reportfilter):
    outf = open(outfile, "a+")
    mergedTotal = lsuTotal + ssuTotal
    outf.write("\nTotal rRNA: " + str(mergedTotal) + "\n")
    rfil = reportfilter
    includefilter = 0.01
    
    global mergedPiedict
    mergedPiedict = dict()

    sortedLSU = sorted(lsuCategories)
    sortedSSU = sorted(ssuCategories)
    for litem in sortedLSU:
        ldomain = litem.split(".")[0]
        lname = lsuCategories[litem][1]
        lnumber = lsuCategories[litem][0]
        lpct = round(lnumber / lsuTotal * 100, 2)

        # There are a few weird cases where you see the same names at a WAY more specific taxonomic level
        # One such example (and why I included this) is when Incertae sedis is listed
        # One fungus was listed as Ascomycota and then categorized up to the species level,
        #   but then followed by "Incertae sedis; Ascomycota", to indicate that it's not certain if it belongs there.
        # We don't want BOTH "Ascomycota"s in "Eukaryota;Fungi;Ascomycota;...;Incertae sedis;Ascomycota" to be added.
        #   and due to sorting this script goes through the broader tax levels first,
        #   so we just add together the first match
        # In larger datasets, this will also help with runtime by ending the loop immediately after a match occurs
        #   (But compared to the reads or even Metaxa.taxonomy files, these dictionaries should be a lot smaller;
        #    This is also the reason why I opted for a nested list for simplicity instead of a binary search)
        match = False

        for sitem in sortedSSU:
            sdomain = sitem.split(".")[0]
            sname = ssuCategories[sitem][1]
            snumber = ssuCategories[sitem][0]
            spct = round(snumber / lsuTotal * 100, 2)

            if sname == lname and sdomain == ldomain: #and lpct >= includefilter and spct >= includefilter:
                match = True
                mergednumber = snumber + lnumber
                mergedpct = round(mergednumber / mergedTotal * 100, 2)
                color = makecolor(mergedpct)
                if sname != "" and lname != "" and mergedpct >= rfil or sdomain == sname:
                    tabs = litem.count('.')
                    outf.write(fg(color) + "    "*tabs + lname + ": " + str(mergednumber) + " (" + str(mergedpct) + "%)\n" + fg.rs)
                
                mergedPiedict[litem] = [lname, mergedpct]
            if match == True:
                break
    outf.close()                


"""Function manually called to do all the reporting work"""
def main_results(lsupath, ssupath, filter):
    if lsupath:
        print("lsu provided, summarizing results...")
        TaxParser(lsupath, "lsu", "w")
    if ssupath:
        print("ssu provided, summarizing reuslts...")
        TaxParser(ssupath, "ssu", "a")
    if lsupath and ssupath:
        print("merging lsu and ssu...")
        rRNA_Merge(filter)
    print("Done!")

if __name__ == "__main__":
    main_results(lsupath, ssupath, filter)
    exit(0)
