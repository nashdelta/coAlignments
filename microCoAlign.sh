#!/bin/bash

# NOTE: This script should be run after MATLAB processing (processFamilies) of the ouput from blast2cluster.sh
#
# $1 is the plant shortlisted interaction file (host family index, host protein, viral protein, tab-delimited)
# $2 is the vertebrate shortlisted interaction file
# $3 is the name of the total host .sr library
# $4 is the name of the total virus .sr library
#
# This script takes shortlisted interactions files and outputs co-alignments (in separate files) to review
# cases where there are more than one viral protein associated with a given host family. This information
# determines whether to do a single psiblast search for the viral family, or separate searches for each
# individual viral protein query (which consequently will "downgrade" the host protein - virus protein pairs).

#get master tab file
cat $1 > tmp_tab.txt
cat $2 > tmp_tab2.txt
cut -f 1 tmp_tab2.txt > tmp_1.txt
cut -f 2,3 tmp_tab2.txt > tmp_2.txt
sed -i 's/[^0-9]//g' "tmp_1.txt"
sed -i 's/^/b/g' "tmp_1.txt"
paste tmp_1.txt tmp_2.txt > tmp_3.txt
cat tmp_3.txt >> tmp_tab.txt
cut -f 1 tmp_tab.txt | sort -u -V > tmp_index.txt

#get master blast db
cat $3 > tmp_sr.txt
cat $4 >> tmp_sr.txt
sr2fa tmp_sr.txt -w=0 > tmp_1.fa
makeblastdb -in tmp_1.fa -input_type fasta -dbtype prot -parse_seqids -out tmp_Blast

#create .cls files
for i in $(cat tmp_index.txt); do
	echo $i > tmp_r.txt
	tab_select tmp_tab.txt -t=tmp_r.txt > tmp_3.txt
	cut -f 2 tmp_3.txt | sort -u > tmp_4.txt
	cut -f 3 tmp_3.txt | sort -u > tmp_5.txt
	tr '\n' ' ' < tmp_4.txt > tmp_6.txt
	tr '\r\n' ' ' < tmp_5.txt > tmp_7.txt
	sed -i 's/  / /g' tmp_7.txt
	sed -i 's/.$//g' tmp_6.txt
	sed -i 's/.$//g' tmp_7.txt
	echo "" >> tmp_6.txt
	echo "" >> tmp_7.txt
	cat tmp_6.txt >> tmp_8.txt
	cat tmp_7.txt >> tmp_9.txt
done
paste tmp_index.txt tmp_8.txt > tmp_h.cls
paste tmp_index.txt tmp_9.txt > tmp_v.cls

#generate alignments
cls2ali tmp_h.cls -d=tmp_Blast -sing > tmp_hAli.fa
cls2ali tmp_v.cls -d=tmp_Blast -sing > tmp_vAli.fa

#text processing
fa2sr tmp_hAli.fa -w=0 > reviewHOSTsr.txt
fa2sr tmp_vAli.fa -w=0 > reviewVIRUSsr.txt

#clean up
rm tmp_*
