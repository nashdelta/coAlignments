#!/bin/bash
#
# $1 is the target host tab virus protein query list
# $2 is the name of the total host .sr library
# $3 is the name of the total virus .sr library

#create output directory

mkdir targetPsiBlast

#get master .sr and target list
cat $2 > tmp_sr.txt
cat $3 >> tmp_sr.txt
cut -f 1 $1 > tmp_list.txt
cut -f 2 $1 >> tmp_list.txt

tab_select tmp_sr.txt -t=tmp_list.txt > tmp_sr2.txt
sr
cut -f 1 tmp_sr2.txt > tmp_list.txt

cp ~/virusCoEvoProject/tmp_sr2.txt ~/virusCoEvoProject/targetPsiBlast/tmp_sr2.txt
cp ~/virusCoEvoProject/tmp_list.txt ~/virusCoEvoProject/targetPsiBlast/tmp_list.txt

#clean up
rm tmp_*

#enter target directory

cd targetPsiBlast

for i in $(cat tmp_list.txt); do
	echo $i > tmp_r.txt
	tab_select tmp_sr2.txt -t=tmp_r.txt > tmp_sr.txt
	sr2fa tmp_sr.txt > tmp_fa.txt
	psiblast -db nr -query tmp_fa.txt -outfmt "6 delim=@ qaccver saccver" -comp_based_stats 0 -seg 'no' -max_target_seqs 100 > tmp.txt
done

#clean up
#rm tmp_*



