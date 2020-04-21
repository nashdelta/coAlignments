#!/bin/bash

# $1 is the name of the target directory

# This script takes an input directory containing text files with two columns
# of matched accession id's, host/virus, and generates alignments for the id's in each column.
# the output are returned to separate directories. The files in the directory are assumed
# to be named in the convention of 1_*,...,n_*.

cd $1
FILES=*.txt
count=0
for i in $FILES
do

count=$(($count+1))
mkdir ali_$count

cut -f 1 $i > tmp_1.txt
sort -u tmp_1.txt > tmp_2.txt
blastdbcmd -db nr -target_only -entry_batch tmp_2.txt > tmp_1.fa
makeblastdb -in tmp_1.fa -dbtype prot -parse_seqids -out tmp_db
nemp tmp_1.fa -r=">(\S+)" > tmp_1.lst

prof_align tmp_1.lst -db= tmp_db -bare -w=0 -n= host_$count -l= 0.5 -u=2.3 -gcut=0.667 -i=10

mv host_$count* ali_$count

rm tmp_*

cut -f 2 $i > tmp_1.txt
sort -u tmp_1.txt > tmp_2.txt
blastdbcmd -db nr -target_only -entry_batch tmp_2.txt > tmp_1.fa
makeblastdb -in tmp_1.fa -dbtype prot -parse_seqids -out tmp_db
nemp tmp_1.fa -r=">(\S+)" > tmp_1.lst

prof_align tmp_1.lst -db= tmp_db -bare -w=0 -n= virus_$count -l= 0.5 -u=2.3 -gcut=0.667 -i=10

mv virus_$count* ali_$count

rm tmp_*

done

