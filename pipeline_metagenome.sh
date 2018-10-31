#!/bin/sh

usage(){
echo "
Written by Conor O'Donoghue
Last modified Aug 16, 2018

Description: Pipeline for analyzing taxonomy of metagenomic data.
             Steps: Subsample --> mapping --> taxonomic summary --> visualization
             Requires Python 3
             Python dependencies: krona

Usage: test.sh [-i <filepath>] [-o <string>] [OPTIONS]

Input Parameters:
-i|--in <file>		REQUIRED. Path to input file.
--ref <file>		REQUIRED unless skipmapping=true. Path to reference database.

--filter 0.2		If a taxa makes up less than <filter>% of the sample, it isn't reported. 
			Default 0.2

--subsample t		Subsample reads to a given value before running pipeline. Highly recommended as a raw library
			has many more reads than is necessary to generate an accurate taxonomic breakdown.
			Subsamples as default. Set to false to skip. See samplenum option to set exact number of reads to use. 
--samplenum 400000	Number of reads desired in subsample. Default set to 400000. Ignore if subsample turned off.

--reformatopts <args>	reformat.sh has a number of additional options, pass them through here as a string to add them to that command.
                        ex) --bbmapopts 'int=t maxhistlen=10000 bhist=bhist.log' etc
                        if you want to look at reformat.sh options, enter --help as the option, and all this will do is display
                        bbmap options and then exit.
                        ex) --reformatopts '--help' --> reformatopts.sh usage info and exit

--bbmapopts <args>	bbmap.sh has a number of additional options, pass them through here as a string to add them to that command.
			ex) --bbmapopts 'minid=0.82 refstats=refstats.log machineout=t' etc
			if you want to look at bbmap.sh options, enter --help as the option, and all this will do is display
			bbmap options and then exit.
			ex) --bbmapopts '--help' --> bbmap.sh usage info and exit
			
			By default, this pipeline runs bbmap with maxindel=100 as the only optional command.
			the input is the same as the input for the pipeline, te output uses the out base name,
			and the reference used is the silva132 database ref located at
			/global/projectb/scratch/cmodonog/metatranscriptome

--skipmapping f		If set to true, will skip the bbmap step.
			For projects where bbmap has already been run, but the rest of the analysis still needs to be done.
			If you intend to use this, use the bbmap .sam output as the input file.

Output Parameters:
-o|--out <file>		REQUIRED. Output base name
--krona t		Whether or not to produce an interactive 'out.krona.html' visualization of data.

-h|--help		Displays usage info, then exits
"
}

# Setting defaults
FILTER=0.2
KRONA=true
SUBSAMPLE=true
SKIPMAPPING=false
SAMPLENUM=400000
DIR=$(pwd)

# Parse options and turn their arguments into variables.
while true ; do
    case "$1" in
        -i|--in)
            INFILE=$2
            shift 2
            ;;
        -o|--out)
            OUTBASE=$2
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --ref)
            REF=$2
            shift 2
            ;;
        --krona)
            case "$2" in
                true|t|1|y|yes)
                    KRONA=true
                    shift 2
                    ;;
                false|f|0|n|no)
                    KRONA=false
                    shift 2
                    ;;
                *)
                    echo "boolean value expected for --krona"
                    exit 1
                    ;;
            esac
            ;;
        --filter)
            FILTER=$2
            shift 2
            ;;
        --subsample)
            case "$2" in
                true|t|1|y|yes)
                    SUBSAMPLE=true
                    shift 2
                    ;;
                false|f|0|n|no)
                    SUBSAMPLE=false
                    shift 2
                    ;;
                *)
                    echo "boolean value expected for --subsample"
                    exit 1
                    ;;
            esac
            ;;
        --samplenum)
            SAMPLENUM=$2
            shift 2
            ;;
        --bbmapopts)
            case "$2" in
                --help)
                    shifter --image=bryce911/bbtools bbmap.sh --help
                    exit 0
                    ;;
                *)
                    BBMAPOPTS=$2
                    shift 2
                    ;;
            esac
            ;;
        --reformatopts)
            case "$2" in
                --help)
                    shifter --image=bryce911/bbtools reformat.sh --help
                    exit 0
                    ;;
                *)
                    REFORMATOPTS=$2
                    shift 2
                    ;;
            esac
            ;;

        --skipmapping)
            case "$2" in
                true|t|1|y|yes)
                    SKIPMAPPING=true
                    shift 2
                    ;;
                false|f|0|n|no)
                    SKIPMAPPING=false
                    shift 2
                    ;;
                *)
                    echo "boolean value expected for --skipmapping"
                    exit 1
                    ;;
            esac
            ;;

        --) shift ; break ;;
        ''|' ') break ;;
        *) echo "$1 is not an accepted option!" ; exit 1 ;;
    esac
done


if [[ ${INFILE} = '' ]] || [[ ${OUTBASE} = '' ]]; then
    echo "Infile and output base name required! If mapping, reference also required! See --help for usage details."
    exit 1
fi

if [ ${REF} = ''] && [ ${SKIPMAPPING} = false ]; then
    echo "If mapping, reference is required! See --help for usage details."
    exit 1
fi

# Subsampling step
if [[ ${SUBSAMPLE} = true ]]
then
    echo "Subsampling reads..."
    SUBSAMPLEOUT="${OUTBASE}.subsample.${SAMPLENUM}.fastq"
    (( SAMPLENUM = SAMPLENUM / 2 ))
    time shifter --image=bryce911/bbtools reformat.sh in=${INFILE} out=${SUBSAMPLEOUT} samplereadstarget=${SAMPLENUM} ${REFORMATOPTS}
else
    echo "Subsampling will be skipped."
    SUBSAMPLEOUT=${INFILE}
fi 

# Mapping step
if [[ ${SKIPMAPPING} = true ]]
then
    echo "Mapping will be skipped."
    BBMAPOUT=${INFILE}
else
    echo "Mapping reads to silva database..."
    BBMAPOUT="${OUTBASE}.bbmap.sam"
    time shifter --image=bryce911/bbtools bbmap.sh in=${SUBSAMPLEOUT} out=${BBMAPOUT} path=${REF} maxindel=100 ${BBMAPOPTS}
    BBMAPOUT="${OUTBASE}.bbmap.sam"
fi

# Taxonomic summarization step
echo "Generating taxonomic summary..."
python bbmaptax.py --in ${BBMAPOUT} --out ${OUTBASE} --filter ${FILTER} --krona ${KRONA}

KRONAMERGED=${OUTBASE}.krona.merged.txt
KRONASSU=${OUTBASE}.krona.ssu.txt
KRONALSU=${OUTBASE}.krona.lsu.txt

# Visualization step
if [[ ${KRONA} = true ]]
then
    echo "Generating krona charts..."
    ktImportText ${KRONAMERGED} -o ${OUTBASE}.krona.merged.html
    #/global/projectb/scratch/cmodonog/scripts/makeKrona.sh -i ${KRONAMERGED} -o ${OUTBASE}.krona.merged.html --rrnatype merged
    ktImportText ${KRONASSU} -o ${OUTBASE}.krona.ssu.html
    #/global/projectb/scratch/cmodonog/scripts/makeKrona.sh -i ${KRONASSU} -o ${OUTBASE}.krona.ssu.html --rrnatype ssu
    ktImportText ${KRONALSU} -o ${OUTBASE}.krona.lsu.html
    #/global/projectb/scratch/cmodonog/scripts/makeKrona.sh -i ${KRONALSU} -o ${OUTBASE}.krona.lsu.html --rrnatype lsu
else
    echo "Krona chart generation will be skipped."
fi

echo "Done!"
