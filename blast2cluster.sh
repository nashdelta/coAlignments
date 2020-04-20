#!/bin/bash

# NOTE: This script should be run after mmseq2blast.sh
#
# $1 is the name/prefix of the input psiblast results
# $2 is the original input .sr file
# $3 is the minimum coverage for a good blast hit
# $4 is the minimum bit score for a good blast hit
# $5 is the minumum number of taxa in a family for further processing
# $6 is the minumum number of sequences in a family for further processing
#
# This script takes the input $1PB.txt file and calculates stats after the $3/4 thresholds are applied,
# #taxa/sequences in the family. These non-exclusive families are then sorted by the number of
# taxa contained. Exclusive families are then assembled through a "greedy" search where sequences appearing in 
# larger clusters (by tax count not sequence count) are removed from smaller clusters. The output
# are directories of the resulting tax-sorted .sr files and word count files *_taxSortFam,
# *_MasterWC, and $1TCWC Finally, it returns .cls files for tax sorted families, $1TC.cls.

#make output directory for tax-sorted results
mkdir $1'_'$3'_'$4'_taxSortFam'

#enter mmseqs results directory
cd $1'_mmseqSR'

#text processing
cut -f 1,2 $1'PB'.txt > tmp_names.txt
cut -f 3 $1'PB'.txt > tmp_qLen.txt
cut -f 4 $1'PB'.txt > tmp_sLen.txt
cut -f 5 $1'PB'.txt > tmp_qStart.txt
cut -f 6 $1'PB'.txt > tmp_qEnd.txt
cut -f 7 $1'PB'.txt > tmp_sStart.txt
cut -f 8 $1'PB'.txt > tmp_sEnd.txt
paste -d- tmp_qEnd.txt tmp_qStart.txt | bc -l > tmp_qFoot.txt
paste -d- tmp_sEnd.txt tmp_sStart.txt | bc -l > tmp_sFoot.txt
cut -f 10 $1'PB'.txt > tmp_bitScore.txt
paste -d/ tmp_qFoot.txt tmp_qLen.txt | bc -l > tmp_fOverQ.txt
paste -d/ tmp_sFoot.txt tmp_sLen.txt | bc -l > tmp_fOverS.txt
paste tmp_names.txt tmp_fOverQ.txt tmp_fOverS.txt tmp_bitScore.txt > tmp_PB.txt

#select rows with good coverage and bit score or bit score close to 1
#to ensure "self-matches" are included
tab_select tmp_PB.txt -k1=3 -kcut=$3 > tmp_gCov1.txt
tab_select tmp_gCov1.txt -k1=4 -kcut=$3 > tmp_gCov.txt
tab_select tmp_gCov.txt -k1=5 -kcut=$4 > tmp_gBit.txt
tab_select tmp_PB.txt -k1=5 -kcut=0.95 >> tmp_gBit.txt

#make families unique (i.e. associate matches to different regions of a single pair)
#for each family get the number of taxa, the number of sequences, and record information in a word count table
cut -f 1,2 tmp_gBit.txt | sort -u -k1,2 > tmp_match.txt
cut -f 1 tmp_match.txt | uniq > tmp_m1.txt
for i in $(cat tmp_m1.txt); do
	echo $i > tmp_r.txt
	tab_select tmp_match.txt -t=tmp_r.txt | cut -f 2 > tmp_1.txt
	grep -Ev 'CONSENSUS' tmp_1.txt > tmp_nConSeq.txt
	awk '/CONSENSUS/' tmp_1.txt > tmp_cons.txt
	sed -i "s/CONSENSUS/$1/g" "tmp_cons.txt"
	sed -i 's/$/.sr/' "tmp_cons.txt"
	cat tmp_nConSeq.txt > tmp_allSeq.txt
	for j in $(cat tmp_cons.txt); do
		cat $j >> tmp_allSeq.txt
	done
	cut -f 1 -d '.' tmp_allSeq.txt | sort -u > tmp_uTax.txt
	seqNum=$(wc -l < tmp_allSeq.txt)
	taxNum=$(wc -l < tmp_uTax.txt) 
	echo $i$'\t'$taxNum$'\t'$seqNum >> $1'WC'.txt
done

#get unique pairs indexed by query MSA or singleton
cut -f 1,2 tmp_gBit.txt | sort -u -k1,2 > tmp_uMatch.txt

#select families meeting taxa and sequence number thresholds
#sort them based on the number of taxa
#get lists of family input queries matching these criteria
tab_select $1WC.txt -k1=2 -kcut=$5 > tmp_1.txt
tab_select tmp_1.txt -k1=3 -kcut=$6 > tmp_2.txt
sort -n -r tmp_2.txt -k2 > tmp_wc.txt
cut -f 1 tmp_wc.txt > tmp_newList.txt

#copy original sequence database into folder containing mmseqs .sr files
cp ~/virusCoEvoProject/$2 ~/virusCoEvoProject/$1_mmseqSR/tmp_rSeq.txt

#construct greedy, exclusive clusters sorted by size i.e. sequences appearing in larger clusters
#cannot reappear in smaller clusters

#tax-sort
loopNum=0
for i in $(cat tmp_newList.txt); do
	echo $i > tmp_r.txt
	tab_select tmp_uMatch.txt -t=tmp_r.txt | cut -f 2 > tmp_1.txt
	grep -Ev 'CONSENSUS' tmp_1.txt > tmp_nConSeq.txt
	awk '/CONSENSUS/' tmp_1.txt > tmp_cons.txt
	sed -i "s/CONSENSUS/$1/g" "tmp_cons.txt"
	sed -i 's/$/.sr/' "tmp_cons.txt"
	cat tmp_nConSeq.txt > tmp_allSeq.txt
	for j in $(cat tmp_cons.txt); do
		cat $j >> tmp_allSeq.txt
	done
	tab_select tmp_rSeq.txt -t=tmp_allSeq.txt > tmp_1.sr
	tmpWC=$(wc -l < tmp_1.sr)
	if (( $tmpWC \> 0 )); then
		loopNum=$(echo $loopNum + 1 | bc)
	   	cp ~/virusCoEvoProject/$1_mmseqSR/tmp_1.sr ~/virusCoEvoProject/$1'_'$3'_'$4'_taxSortFam'/$1TC.$loopNum.sr
	   	tab_select tmp_rSeq.txt -t=tmp_allSeq.txt -n > tmp_rSeq2.txt
	   	cat tmp_rSeq2.txt > tmp_rSeq.txt
	   	cut -f 1 -d '.' tmp_1.sr | sort -u > tmp_uTax.txt
	   	seqNum=$(wc -l < tmp_1.sr)
	   	taxNum=$(wc -l < tmp_uTax.txt)
		echo $i$'\t'$1TC.$loopNum.sr$'\t'$taxNum$'\t'$seqNum >> tmp_bigWC.txt
	else
	   	echo $i$'\t'NA$'\t'0$'\t'0 >> tmp_bigWC.txt
	fi
done
sort -n -r tmp_bigWC.txt -k3 > tmp_WC.txt

#Get remaining input queries from $1WC.txt
tab_select $1WC.txt -n -t=tmp_WC.txt > tmp_1.txt
cut -f 1 tmp_1.txt > tmp_2.txt
cat tmp_2.txt > tmp_3.txt
sed -i 's/$/\tNA\t0\t0/' "tmp_3.txt"
cat tmp_3.txt >> tmp_WC.txt
cp ~/virusCoEvoProject/$1_mmseqSR/tmp_WC.txt ~/virusCoEvoProject/$1'_'$3'_'$4'_taxSortFam'/$1TCWC.txt

#Get composite word count file: blast tax#seq#, TC tax#&seq#, mmseqs query, TC label,

sort $1WC.txt > tmp_1.txt
sort ~/virusCoEvoProject/$1'_'$3'_'$4'_taxSortFam'/$1TCWC.txt > tmp_2.txt
cut -f 1 tmp_1.txt > tmp_name1.txt
cut -f 2 tmp_2.txt > tmp_name2.txt
cut -f 2,3 tmp_1.txt > tmp_n1.txt
cut -f 3,4 tmp_2.txt > tmp_n2.txt

paste tmp_n1.txt tmp_n2.txt tmp_name1.txt tmp_name2.txt > tmp_1.txt
sort -n -r tmp_1.txt -k2 > tmp_2.txt
cat tmp_2.txt > ../$1_$3_$4_MasterWC.txt

#clean up
rm tmp_*

#Create .cls files
cd ../$1'_'$3'_'$4'_taxSortFam'
for file in `ls -v $1TC.*.sr`; do
	cut -f 1 $file > tmp_1.txt
	cat tmp_1.txt >> tmp_5.txt
	tr '\n' ' ' < tmp_1.txt > tmp_2.txt
	echo "" >> tmp_2.txt
	cat tmp_2.txt >> tmp_3.txt
	echo $file >> tmp_4.txt
done

#append with remaining sequences if not all are clustered
cp ~/virusCoEvoProject/$2 ~/virusCoEvoProject/$1'_'$3'_'$4'_taxSortFam'/tmp_allSeq.txt
cut -f 1 tmp_allSeq.txt > tmp_allN.txt
tab_select tmp_allN.txt -t=tmp_5.txt -n > tmp_remain.txt
tr '\n' ' ' < tmp_remain.txt > tmp_r2.txt
echo "" >> tmp_r2.txt
echo remain >> tmp_4.txt
cat tmp_r2.txt >> tmp_3.txt
paste tmp_4.txt tmp_3.txt > $1TC.cls

#clean up
rm tmp_*