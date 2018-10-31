#!/bin/sh
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
        --rrnatype)
            RRNATYPE=$2
            shift 2
            ;;
        --) shift ; break ;;
        ''|' ') break ;;
        *) echo "$1 is not an accepted option!" ; exit 1 ;;
    esac
done

ktImportText ${INFILE} -o ${OUTBASE}.krona.${RRNATYPE}.html
