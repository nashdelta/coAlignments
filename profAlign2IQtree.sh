#!/bin/bash

# $1 is the number of virus/host pairs (of .sr files) in the directory
# $2 is the name of the target directory

# This script takes pairs of input .sr files, host and virus alignments, calculates the mean sequence length among
# both hosts and viruses and prints to a file for later use within INDELible, and prepares files to run IQtree,
# both the input .fa files and the .bat file to start the program (on PC - so note the files must be
# manually moved to the PC project directory).

cd $2
vals=$(seq 1 1 $1)

echo 'cd C:\Users\rochmannd\Desktop\commandLineResources\iqtree\bin' > $2.bat

for i in $vals
do

sr2fa host_$i.sr > host_$i.fa
sr2fa virus_$i.sr > virus_$i.fa
cat host_$i.fa > tmp_1.fa
cat virus_$i.fa > tmp_2.fa
sed -i 's/-//g' tmp_1.fa
sed -i 's/-//g' tmp_2.fa
fa2len tmp_1.fa -s | sort > tmp_1.txt
fa2len tmp_2.fa -s | sort > tmp_2.txt
num0=$(wc -l < tmp_1.txt)
num1=$(($num0 / 2 + 1))
num0=$(wc -l < tmp_2.txt)
num2=$(($num0 / 2 + 1))
sed "${num1}q;d" tmp_1.txt >> tmp_i1.txt
sed "${num2}q;d" tmp_2.txt >> tmp_i2.txt

echo 'iqtree -fast -nt AUTO -st AA -mset Poisson,JTT,JTTDCMut,Dayhoff,DCMut,WAG,mtMAM,mtART,mtREV,rtREV,cpREV,VT,Blosum62,LG,HIVb,HIVw -mrate I+G -redo -s C:\Users\rochmannd\Desktop\co_Alignment_PC\'$2'\host_'$i'.fa' >> $2.bat
echo 'iqtree -fast -nt AUTO -st AA -mset Poisson,JTT,JTTDCMut,Dayhoff,DCMut,WAG,mtMAM,mtART,mtREV,rtREV,cpREV,VT,Blosum62,LG,HIVb,HIVw -mrate I+G -redo -s C:\Users\rochmannd\Desktop\co_Alignment_PC\'$2'\virus_'$i'.fa' >> $2.bat

done

paste tmp_i1.txt tmp_i2.txt > indel.tmp.txt

rm tmp_*