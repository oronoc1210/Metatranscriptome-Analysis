import argparse
import sys
import os
import logging

def parser_gen():
    parser = argparse.ArgumentParser(
        description="""Takes output of bbmap mapping a metagenome
        to silva database and parses the data to generate a taxonomic summary.
        Requires Python 3 and krona perl module in PATH.""",
        epilog="""Example command: bbmaptax.py -i example.txt"""
        )
    parser.add_argument(
        "-i", "--input", dest='infile', required=True
        help="path to input file")
    parser.add_argument(
        "-o", "--output", dest='outbase', default='metatax'
        help="output base name")
    parser.add_argument(
        "--filter", type=float, dest='filter', default=0.2, 
        help="""If a particular taxa makes up less than <filter> percent 
        of the total reads, it isn't reported. Exception for the Domain level,
        which will always be reported.""")
    parser.add_argument(
        "--krona", dest='krona', default="True", 
        help="If True, will generate an interactive chrona chart from data.")
    args = parser.parse_args()
    return args

total = 0
lsuTotal = 0
ssuTotal = 0
lsuCategories = dict()
ssuCategories = dict()
sortedLsuCategories = dict()
sortedSsuCategories = dict()

taxOut = outbase + ".taxstats.txt"
kronaOut = outbase + ".krona.txt"

inf = open(infile, "r")
taxOutf = open(taxOut, "w+")

"""
 I wish I could do something a bit more sophisticated with the colors, but the terminal doesn't support truecolor,
 so any option that uses (r,g,b) doesn't work. 
 All the color modules supported by anaconda either use truecolor or just have about 8 basic colors,
 so I just made my own janky 8-bit colormap.
"""
reset = '\033[0m'
def makeColor(percent):
    if percent < 2:
        color = '\033[38;5;230m'
    elif percent < 8:
        color = '\033[38;5;229m'
    elif percent < 20:
        color = '\033[38;5;228m'
    elif percent < 40:
        color = '\033[38;5;227m'
    elif percent < 60:
        color = '\033[38;5;221m'
    elif percent < 70:
        color = '\033[38;5;215m'
    elif percent < 80:
        color = '\033[38;5;209m'
    elif percent < 90:
        color = '\033[38;5;202m'
    else:
        color = '\033[38;5;196m'
    return color

""" 
 Takes the bbmap .sam output and summarizes results by taxonomic level 
"""
def TaxParser():
    global total      # This is probably the most embarrassing part of this script.
    global lsuTotal   # I need the totals and dictionaries in other functions, so I just make them global.
    global ssuTotal
    for line in inf:
        if line.startswith( "@" ) or line == "\n":  # .sam files include many non-read lines that start with '@' (including the entire reference database), or empty lines. Want to ignore these.
            continue
        else:
            total += 1
            idRegion=line.split('\t')[2] # read lines are reported as follows: [readID]\t[score]\t[taxID]\t[lots of other things separated by more tabs]. We want the taxID ragion.
            if "tid" not in idRegion:    # taxonomic classifications start with 'tid'. If no classification, it just moves on.
                continue
            if "LSU" in idRegion:
                rna = "LSU"
                lsuTotal += 1
            if "SSU" in idRegion:
                rna = "SSU"
                ssuTotal += 1
            parts=idRegion.split()       # idregion has two parts. the taxid and rRNA type separated by '|'s, then the tax classification separated by ';'s
            taxa = ' '.join(parts[1:])   # Annoyingly, these parts are separated by spaces even though there are spaces in the tax classification.
            taxaList = taxa.split(";")   # So, need to just remove everything before the first space to get the full tax list.

            for i in range(1,len(taxaList)+1):      # Want Eukaryota, Eukaryota;Fungi, Eukaryota;Fungi;Dikarya, etc to have different categories to be able to count broader tax levels effectively.
                fullName = ';'.join(taxaList[0:i])  # e.g. if function sees Eukaryota;Fungi and Eukaryota;Chloroplastida, I want "Eukaryota" to go up by 2 and the two kingdoms to go up 1 each.
                levelName = taxaList[i-1]           # Having the FULL classification rather than just the last name is important for sorting and reporting later. Hence fullName as well as levelName.
                if rna == "LSU":
                    if fullName not in lsuCategories:
                        lsuCategories[fullName] = [1, levelName]
                    else:
                        lsuCategories[fullName][0] += 1
                if rna == "SSU":
                    if fullName not in ssuCategories:
                        ssuCategories[fullName] = [1, levelName]
                    else:
                        ssuCategories[fullName][0] += 1

    global sortedLsu
    sortedLsu=sorted(lsuCategories)
    global sortedSsu
    sortedSsu=sorted(ssuCategories)
    lsuPct = round(lsuTotal / total * 100, 2)
    ssuPct = round(ssuTotal / total * 100, 2)

    taxOutf.write("\nTotal LSU rRNA: " + str(lsuTotal) + " (" + str(lsuPct) + "% of reads)\n\n")
    for category in sortedLsu:
            tabs = category.count(';')
            number = lsuCategories[category][0]
            name = lsuCategories[category][1]
            percentage = round(number / lsuTotal * 100 , 2)
            color = makeColor(percentage)
            
            if name != "" and percentage >= filter or tabs == 0:
                    taxOutf.write(color + "    "*tabs + name + ": " + str(number) + " (" + str(percentage) + "% of rRNA)\n" + reset)
    taxOutf.write("\n" + "-"*50 + "\n")

    taxOutf.write("\nTotal SSU rRNA: " + str(ssuTotal) + " (" + str(ssuPct) + "% of reads)\n\n")
    for category in sortedSsu:
            tabs = category.count(';')
            number = ssuCategories[category][0]
            name = ssuCategories[category][1]
            percentage = round(number / ssuTotal * 100 , 2)
            color = makeColor(percentage)

            if name != "" and percentage >= filter or tabs == 0:
                    taxOutf.write(color + "    "*tabs + name + ": " + str(number) + " (" + str(percentage) + "% of rRNA)\n" + reset)
    taxOutf.write("\n" + "-"*50 + "\n")


"""
Takes the results of lsu and ssu parsing and merges them, where possible
NOTE: this will only report taxa that are present in BOTH LSU and SSU
i.e. even if something is 20% of LSU, if it's not in SSU at all it isn't reported here.
This is why the program also produces LSU and SSU results separately.
See function below for why it's not good to just add together LSU and SSU and why it HAS to be only matches
"""


def rRNA_Merge(reportfilter):
    global mergedCategories
    mergedCategories = dict()
    global mergedTotal
    mergedTotal = lsuTotal + ssuTotal
    mergedPct = round(mergedTotal/total * 100, 2)
    taxOutf.write("\nTotal rRNA: " + str(mergedTotal) + " (" + str(mergedPct) + "% of reads)\n\n")
    rfil = reportfilter
    includefilter = 0.01
    
    for lItem in sortedLsu:
        lDomain = lItem.split(";")[0]
        lName = lsuCategories[lItem][1]
        lNumber = lsuCategories[lItem][0]
        lPct = round(lNumber / lsuTotal * 100, 4)
        match = False

        # First, see if lItem matches with anything in SSU. If it does, add counts together and append to merged dict.
        for sItem in sortedSsu:
            sDomain = sItem.split(";")[0]
            sName = ssuCategories[sItem][1]
            sNumber = ssuCategories[sItem][0]
            sPct = round(sNumber / ssuTotal * 100, 4)
             
            if sName == lName and sDomain == lDomain and lPct >= includefilter and sPct >= includefilter:
                match = True
                mergedNumber = sNumber + lNumber
                mergedPct = round(mergedNumber / mergedTotal * 100, 2)
                mergedCategories[lItem] = [mergedNumber, lName]
                color = makeColor(mergedPct)
                if sName != "" and lName != "" and mergedPct >= rfil or sDomain == sName:
                    tabs = lItem.count(';')
                    taxOutf.write(color + "    "*tabs + lName + ": " + str(mergedNumber) + " (" + str(mergedPct) + "% of rRNA)\n" + reset)
                
            if match == True:  #If match already occurred, no reason to keep searching. Move on.
                break
    global sortedMerged
    sortedMerged = sorted(mergedCategories)

"""
!!! NOT USED !!!.
This truly adds together SSU and LSU into one dictionary
i.e. even if something in SSU isn't in LSU at all, it gets added.
Now merged is everything in LSU and SSU, with matches having their totals added together.
This is a bad solution, because LSU and SSU are classified differently.
Many of the 'middle' taxa levels are present in only SSU and not LSU. 
So species from only SSU and only LSU get put in different places when they shouldn't be.
"""
def True_rRNA_Merge(repFil):
    global mergedCategories
    mergedCategories = dict()
    global mergedTotal
    mergedTotal = lsuTotal + ssuTotal
    mergedPct = round(mergedTotal/total * 100, 2)
    taxOutf.write("\nTotal rRNA: " + str(mergedTotal) + " (" + str(mergedPct) + "% of reads)\n\n")
    matchFil = 0.01

    for lItem in sortedLsu:
        lDomain = lItem.split(";")[0]
        lName = lsuCategories[lItem][1]
        lNumber = lsuCategories[lItem][0]
        lPct = round(lNumber / lsuTotal * 100, 4)
        match = False

        # First, see if lItem matches with anything in SSU. If it does, add counts together and append to merged dict.
        for sItem in sortedSsu:
            sDomain = sItem.split(";")[0]
            sName = ssuCategories[sItem][1]
            sNumber = ssuCategories[sItem][0]
            sPct = round(sNumber / ssuTotal * 100, 4)
            if sName == lName and sDomain == lDomain and lPct >= matchFil and sPct >= matchFil:
                match = True
                mergedNumber = sNumber + lNumber
                mergedPct = round(mergedNumber / mergedTotal * 100, 2)
                mergedCategories[lItem] = [mergedNumber, lName]
                color = makeColor(mergedPct)
                #if sName != "" and lName != "" and mergedPct >= repFil or sDomain == sName:
                 #   tabs = lItem.count(';')
                  #  taxOutf.write(color + "    "*tabs + lName + ": " + str(mergedNumber) + " (" + str(mergedPct) + "% of rRNA)\n" + reset)
            
            if match == True:    # if match already occurred, no need to search through the rest of the dict. Move on to the next LSU item.
                break
        if match == False:    # if after searching through all of SSU there was no match, 
            if lPct > repFil:     # if the count is actually relatively high
                newName = lName + ' (only LSU)'    #Add it to merged dict as well, but mention that it only comes from LSU.
                mergedCategories[lItem] = [lNumber, newName]

    # we need to do this again for SSU, to add the sizable elements from it that aren't in LSU. This involves matching again, so that we DON'T add matches (which were already added)
    # That said, we don't need to search through all of lsuCategories, just the mergedCategories. It already has all the matches there and should be a good bit smaller.
    for sItem in sortedSsu:
        sDomain = sItem.split(";")[0]
        sName = ssuCategories[sItem][1]
        sNumber = ssuCategories[sItem][0]
        sPct = round(sNumber / ssuTotal * 100, 4)
        match = False
        if sName == "Spermatophyta":
            print(ssuCategories[sItem])
            print(sPct)
            break

        # First, see if sItem is already in mergedCategories. If it is, move on.
        for mItem in mergedCategories:
            mDomain = mItem.split(";")[0]
            mName = mergedCategories[mItem][1]
            if sName == mName and sDomain == mDomain:
                match = True
                break
        if match == False:    # if not already in mergedCategories
            if sPct > repFil:     # if the count is actually relatively high
                newName = sName + ' (only SSU)'    #Add it to merged dict as well, but mention that it only comes from LSU.
                mergedCategories[sItem] = [sNumber, newName]

    global sortedMerged
    sortedMerged = sorted(mergedCategories)
    for mItem in sortedMerged:
        mName = mergedCategories[mItem][1]
        mNumber = mergedCategories[mItem][0]
        mPct = round(mNumber / mergedTotal * 100, 4)
        color = makeColor(mPct)
        tabs = mItem.count(';')
        taxOutf.write(color + "    "*tabs + mName + ": " + str(mNumber) + " (" + str(mPct) + "% of rRNA)\n" + reset)


def kronaGen(taxDict, sortDict, rRNA_type, typeTotal):
    kronaDict = dict()
    kronaOut = outbase + ".krona." + rRNA_type + ".txt"
    kronaOutf = open(kronaOut, "w+")

    for sortItem in sortDict:
        name = taxDict[sortItem][1]
        kronaNum = taxDict[sortItem][0]
        sortTab = sortItem.count(';')
        for numItem in taxDict:
            numTab = numItem.count(';')
            if sortItem in numItem and numTab == sortTab + 1:
                kronaNum = kronaNum - taxDict[numItem][0]
        kronaDict[sortItem] = kronaNum
    if rRNA_type == "lsu" or rRNA_type == "ssu":
        kronaOutf.write(str(total - typeTotal) + "\t" + "non " + rRNA_type + " ribosomal RNA" + "\n")
    else:
        kronaOutf.write(str(total - typeTotal) + "\t" + "non-ribosomal RNA" + "\n")
    for item in sortDict:
        count = kronaDict[item]
        if round(count / typeTotal * 100 , 2) > 0.05:
            kronaOutf.write(str(kronaDict[item]) + "\tribosomal RNA\t" + item.replace(";", "\t") + "\n")
    kronaOutf.close()

"""Function manually called to do all the reporting work"""
def main_taxa(filter):
    print("Parsing data...")
    TaxParser()
    print("Merging lsu and ssu...")
    rRNA_Merge(filter)

def main_krona():
    print("Generating interactive krona charts...")
    kronaGen(mergedCategories, sortedMerged, "merged", mergedTotal)
    os.system("ktImportText %s -o %s.krona.%s.html"
              % (infile, outbase, rrnatype))
    kronaGen(ssuCategories, sortedSsu, "ssu", ssuTotal)
    os.system("ktImportText %s -o %s.krona.%s.html"
              % (infile, outbase, rrnatype))
    kronaGen(lsuCategories, sortedLsu, "lsu", lsuTotal)
    os.system("ktImportText %s -o %s.krona.%s.html"
              % (infile, outbase, rrnatype))
    print("All krona files created!")

if __name__ == "__main__":
    main_taxa(filter)
    if krona == True:
        main_krona()
    print("bbmap analysis complete!")
    exit(0)
