#!/bin/bash

# $1 is the name/prefix of the output files
# $2 is the input .sr file
# $3 is the input similarity threshold for mmseqs
# $4 is whether to send to sge cluster during cls2ali, if no enter '', if yes enter '-sge'
# This script takes an input .sr file, clusters using mmseqs, aligns the clusters, and calculates consensus sequences
# After this, psiblast is run and families are constructed in the {prefix}PB.txt file.

#make output directory
mkdir $1'_mmseqSR'

#convert .sr to fasta and cluster
sr2fa $2 -w=0 > tmp_1.fa
run_mmclust tmp_1.fa -s=$3 -w=0 > tmp_MMC.txt

#text processing
cut tmp_MMC.txt -f 1 > tmp_mmNum.txt
sr2fa $2 -w=0 > tmp_orig.fa

#make db and align
makeblastdb -in tmp_orig.fa -input_type fasta -dbtype prot -parse_seqids -out tmp_Blast
cls2ali tmp_MMC.txt -d=tmp_Blast $4 > tmp_ali.fa

#text processing
fa2sr tmp_ali.fa -w=0 > tmp_ali.sr
cut tmp_ali.sr -f 1 | cut -f 2 -d '|' > tmp_aliIndex.txt
paste tmp_aliIndex.txt tmp_ali.sr > tmp_ali.txt
cut tmp_ali.sr -f 1 | cut -f 4 -d '|' > tmp_clusterID.txt
tab_select $2 -t=tmp_clusterID.txt -n > tmp_single.sr
sort -u tmp_aliIndex.txt > tmp_uIndex.txt

#separate alignments into individual .sr files and calculate consensus
for i in $(cat tmp_uIndex.txt); do
	tab_select tmp_ali.txt -kcut=$i -kbar=$i > tmp_tmp.txt
	cut tmp_tmp.txt -f 2 | cut -f 4 -d '|' > tmp_tmp2.txt
	cut tmp_tmp.txt -f 3 > tmp_tmp3.txt
	paste tmp_tmp2.txt tmp_tmp3.txt > $1'_mmseqSR'/$1.$i.sr

	#set -hcons=0 so that no 'x' characters appear in the consensus
	paste tmp_tmp2.txt tmp_tmp3.txt | sr_filter -hcon=0 -consens > tmp_gap.txt
	sed -i 's/-//g' tmp_gap.txt
	cat tmp_gap.txt >> tmp_c.txt
	echo $i >> tmp_i.txt
done

#text processing
cut tmp_c.txt -f 1 > tmp_c1.txt
cut tmp_c.txt -f 2 > tmp_c2.txt
paste tmp_c1.txt tmp_i.txt -d '.' > tmp_i2.txt

#save singletons and consensuses to .sr file, save separate singleton file
paste tmp_i2.txt tmp_c2.txt > $1'_mmseqSR'/$1'C'.sr
cat tmp_single.sr >> $1'_mmseqSR'/$1'C'.sr
cat tmp_single.sr > $1'_mmseqSR'/$1'S'.sr

#clean up
rm tmp_*

#move into new directory
cd $1'_mmseqSR'

#create blastdb of consensus/singleton sequences
sr2fa $1'C'.sr -w=0 > tmp_1.fa
makeblastdb -in tmp_1.fa -input_type fasta -dbtype prot -parse_seqids -out tmp_Blast

#run psiblast alignment vs. the above db for all alignments

for file in `ls -v $1.*.sr`; do
	sr2fa $file -w=0 > tmp_2.fa
	psiblast -db tmp_Blast -in_msa tmp_2.fa -outfmt "6 delim=@ qaccver saccver qlen slen qstart qend sstart send evalue bitscore" -comp_based_stats 0 -seg 'no' -max_target_seqs 50000 > tmp_out.txt
	cut -f 1 tmp_out.txt > tmp_out1.txt
	awk -F '\t' -v OFS='\t' '{ $(NF+1) = 12345; print }' tmp_out1.txt > tmp_outN.txt
	sed -i "s/12345/${file}/g" "tmp_outN.txt"
	cut -f 2 tmp_outN.txt > tmp_out1.txt
	cut -f 2,3,4,5,6,7,8,9 tmp_out.txt > tmp_out2.txt
	cut -f 10 tmp_out.txt > tmp_out3.txt
	
	#for -in_msa, the highest bit score may not be against self consensus
	selfSr=$(sed '1!d' tmp_out1.txt)
	selfNum=$(echo $selfSr | cut -d. -f2)
	echo CONSENSUS.$selfNum > tmp_s.txt
	tab_select tmp_out.txt -k1=2 -t=tmp_s.txt | cut -f 10 | sort -n -r > tmp_sB.txt
	selfBit=$(sed '1!d' tmp_sB.txt)
	
	awk -F '\t' -v OFS='\t' '{ $(NF+1) = 12345; print }' tmp_out3.txt > tmp_out4.txt
	sed -i "s/12345/$selfBit/g" "tmp_out4.txt"
	cut -f 1 tmp_out4.txt > tmp_aBit.txt
	cut -f 2 tmp_out4.txt > tmp_sBit.txt
	paste -d/ tmp_aBit.txt tmp_sBit.txt | bc -l > tmp_fBit.txt
	paste tmp_out1.txt tmp_out2.txt tmp_fBit.txt >> $1'PB'.txt
done

#run psiblast for singletons against the above db
sr2fa $1'S'.sr -w=0 > tmp_3.fa
psiblast -db tmp_Blast -query tmp_3.fa -outfmt "6 delim=@ qaccver saccver qlen slen qstart qend sstart send evalue bitscore" -comp_based_stats 0 -seg 'no' -max_target_seqs 50000 > tmp_single.txt
cut -f 1 $1'S'.sr | sort -u > tmp_sName.txt
for i in $(cat tmp_sName.txt); do
	echo $i > tmp_r.txt
	tab_select tmp_single.txt -t=tmp_r.txt > tmp_sLittle.txt
	cut -f 1 tmp_sLittle.txt > tmp_out1.txt
	cut -f 2,3,4,5,6,7,8,9 tmp_sLittle.txt > tmp_out2.txt
	cut -f 10 tmp_sLittle.txt > tmp_out3.txt
	selfBit=$(sed '1!d' tmp_out3.txt)
	awk -F '\t' -v OFS='\t' '{ $(NF+1) = 12345; print }' tmp_out3.txt > tmp_out4.txt
	sed -i "s/12345/$selfBit/g" "tmp_out4.txt"
	cut -f 1 tmp_out4.txt > tmp_aBit.txt
	cut -f 2 tmp_out4.txt > tmp_sBit.txt
	paste -d/ tmp_aBit.txt tmp_sBit.txt | bc -l > tmp_fBit.txt
	paste tmp_out1.txt tmp_out2.txt tmp_fBit.txt >> $1'PB'.txt
done

#clean up
rm tmp_*
