#!/bin/sh

# generate pseudoalignment using Kallisto

HOMEDIR="<data dir>"

KALLISTO="<kallisto tool>"

DD="<trimmed data subdirectory>"

find ./$DD -maxdepth 1 -name "ERR*" -print |\
cut -d "/" -f3 |\
cut -d "." -f1 |\
cut -d "_" -f1 |\
sort |\
uniq > filelist.txt

for F in $(cat filelist.txt); do
	FULLSTRING=$F
	SUBSTRING=$(echo $F| cut -d '.' -f2)
	echo "$KALLISTO/kallisto quant --plaintext --bias --threads=2 \
	--index=$HOMEDIR/kindex_ercc \
	-o \"$HOMEDIR/outdir_final_data/${SUBSTRING}\" \
	\"$HOMEDIR/$DD/${FULLSTRING}_1.fastq.gz\" \
	\"$HOMEDIR/$DD/${FULLSTRING}_2.fastq.gz\""

done

