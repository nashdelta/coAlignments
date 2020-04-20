# hostVirus
## This is a step-by-step guide for generating the hostVirus coevolution results. Original code written for this project is contained in the repository but preexisting NCBI or Koonin group resources are not.

## __________________________________________________________________________________

## First retrieve host-protein, virus-protein targets from the STRING interaction db.

First download the STRING protein-protein interaction [database](http://viruses.string-db.org/cgi/download.pl?UserId=D8Isj3SrFhq7&sessionId=JIsh1R8qrtbm) with the link, "protein.links.v10.5.txt.gz". In February 2020 it was about 650 million lines with three fields: protein1, protein2, and score. Score is an integer. Proteins are of the form taxID.proteinID. Note "proteinID" can take a variety of forms.

## Remove rows corresponding to "self-interactions", the vast majority

Run the MATLAB script [readSTRINGdb.m](readSTRINGdb.m) to get the [binaryInteract.mat](binaryInteract.mat) cell with five fields: taxID and ProteinID split plus the score. Only rows corresponding to pairs with two different taxID's are recorded (i.e. pairs within a given organism are excluded). Note that "binaryInteract" saves the original information row-for-row from the STRING download wherever the taxID are different in columns 1&2. It turns out that STRING allows for dupilicated rows insensitive to order. i.e if the order of columns 1 and 2 doesn't matter some rows are not unique.

## Get taxonomy library

[readSTRINGdb.m](readSTRINGdb.m) also returns a text file, allTaxID.txt. Take that file and run:

`taxid2name allTaxID.txt -l > taxLibrary.txt`

in linux. Then run the MATLAB script [processTaxLibrary.m](processTaxLibrary.m) to return the [speciesTax.mat](speciesTax.mat) variable which can then be used with [binaryInteract.mat](binaryInteract.mat) to continue.

## Retrieve host-species, viral-species interaction data

Now download all available protein data from [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?VirusLineage_ss=Viruses,%20taxid:10239&SeqType_s=Nucleotide) with only the fields: Species	Genus	Family Protein	Host. Then for the list of hosts, get the tax lineage from the [portal](https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi), choose "save in file" and return the "tax_report" text file with fields: code	|	name	|	preferred name	|	lineage. Run the MATLAB script [taxLineageReport.m](taxLineageReport.m) to generate the [lineageInfo.mat](lineageInfo.mat) variable for all hosts.

## Generate species interaction bipartite graph and summary table

Now we can proceed given the NCBI Virus [data](ncbiVirusDat.mat) (which was compressed using the MATLAB script [compressVDWN.m](compressVDWN.m)), [lineageInfo.mat](lineageInfo.mat), [binaryInteract.mat](binaryInteract.mat), and [speciesTax.mat](speciesTax.mat) to run the MATLAB script [virusGeneraStats.m](virusGeneraStats.m) which generates basic summary statistics for these datasets. From the NCBI dataset it yields a table containing: Genus	#Plant	#Metazoa ex Vert	#Vertebrates	#Unique Spec/Host	#Viral Spec and plots specifying the number of viral genera with N or more hosts in major clades as well as the number of genera with N or more viral species contained and N or more unique species-hosts pair within that genus. From the STRING database it yields bipartite graphs of viruses and hosts connected with edges based on the number of protein-protein interactions known between the host and the virus and a small summary table with number of species, number of interacting species pairs, and number of protein-protein interactions for major clades. The graphs are generated using the MATLAB script [plotBipartite1.m](plotBipartite1.m).

## Retrieve protein sequences

To proceed we want to retrieve the sequences for each protein appearing in a host-virus pair. To start, download "protein.sequences.v10.5.fa.gz" (STRINGdb.fa) from the STRING database. Some of the proteins are missing. To determine what's missing get a unique list of all protein names in [binaryInteract.mat](binaryInteract.mat) including taxID (for completeness, unused). The first block of the MATLAB script [writeTaxLookup.m](writeTaxLookup.m) will write a file with this information, "subSTRINGnames.txt". Next convert the downloaded .fa file to an .sr file to get the list of names from this file:

`fa2sr STRINGdb.fa -w=0 > STRINGdb.sr`

Now get the rows from the .sr file that match the gene names in "subSTRINGnames.txt" to yield:

`tab_select STRINGdb.sr -t=subSTRINGnames.txt > subSTRING.sr`

The second block of [writeTaxLookup.m](writeTaxLookup.m) uses "subSTRINGnames.txt" and "subSTRING.sr" to write "remainingGene.txt" and "remainingTaxDotGene.txt" for completeness containing the TaxID. this is the input to uniprot.

Next go to the uniprot id mapping [portal](https://www.uniprot.org/uploadlists/), choose "UniProtKB AC/ID to: UniProtKB" and on the return page, add the "sequence" column. Download as tab delimited to "uniprotSSr.txt" only keeping entry and sequence columns. Also download the list of unmapped id's to "uniprotUnmapped1.txt". At the time this was originally done, id's, Q66541_EBVB9 and SCAF_EBVB9 map to the same uniprot page, P03234. The second was manually removed. One id, 1031711.RSPO_c03240, was retrieved manually from the uniprot webpage and added to the .sr file. All others were chain id's of the form, PRO_0000373446.

For the chain id's present in "uniprotUnmapped1.txt" (first remove manually cases like RSPO_c03240), run the third block of [writeTaxLookup.m](writeTaxLookup.m) which will do a web scraping routine, calling each chain id in the https://www.uniprot.org/ search bar, finding the parent protein, opening the child page, and recording the start and end acid for the chain. This outputs two textfiles "uniprotWebscrapResults.txt", which contains the chain id search, the uniprot id returned, the start acid and the end acid where found or NaN's otherwise, and "uUniprotWebscrap.txt" which contains a unique list of the uniprot reference id's returned.

Return to the uniprot id mapping [portal](https://www.uniprot.org/uploadlists/), choose "UniProtKB AC/ID to: UniProtKB" and on the return page, add the "sequence" column. Download as tab delimited to "uniprotPROsr.tab" only keeping entry and sequence columns. There should be no unmapped cases here.

Now run Block 4 of [writeTaxLookup.m](writeTaxLookup.m). Given "uniprotSSr.txt", "subSTRING.sr", "uniprotPROsr.tab", and "uniprotWebscrapResults.txt" create a master .sr file by processing and matching the id's in "subSTRING.sr" to those in "subSTRINGnames.txt". Additionally, match the id's in "uniprotPROsr.tab" to those in "uniprotWebscrapResults.txt". Then take the relevant portion of the sequence from "uniprotPROsr.tab" and pair it with the query id in "uniprotWebscrapResults.txt". These id's are not unique and it's possible that some of the id's correspond to more than one taxID. This is true for the "subSTRINGnames" but does not appear to be true for those id's returned from uniprot. The third part of writeTaxLookup takes the first sequence for each repeated sequence (alphabetical order, but again does not seem to be applicable here) so that each id corresponds to a single sequence, removes any sequences which have id's matched to more than 1 taxID (which also did not happen originally) returns the combined .sr file including the orignal sequences from the STRING database, and returns the list of id's for which no sequence was found. This left 28 genes "lost", from 7 viruses three of which have other "found" genes with interactions and 4 which are totally absent from the final curated database (in the original run).

code	taxid	primary taxid	taxname

Gene lost, others present

1	11053	11053	Dengue virus 1
1	11060	11060	Dengue virus 2
1	11069	11069	Dengue virus 3

All genes lost

1	118140	118140		Teschovirus A
1	12110	12110		Foot-and-mouth disease virus
1	1330524	1330524		Salivirus A
1	312185	312185		Erbovirus A

Now run block 5 of writeTaxLookup to split the .sr files into host and virus documents. Futher split the host documents
by taxa. Proceed to run mmseq2blast for each host group.

# $1 is the name/prefix of the output files
# $2 is the input .sr file
# $3 is the input similarity threshold for mmseqs
# $4 is whether to send to sge cluster during cls2ali, if no enter '', if yes enter '-sge'
# $5 is -grcut for making the consensus sequences
# $6 is the minimum coverage for a good blast hit
# $7 is a minimum bit score for a good blast hit
# This script takes an input .sr file, clusters using mmseqs, aligns the clusters, and calculates consensus sequences
# After this, psiblast is run and families are constructed in the {prefix}PB.txt file. Statistics for each family are
# recorded in {prefix}WC.txt with the input query (MSA or singleton), number of taxa, and number of sequences for the family.

Given the mmseqs results, run blast2cluster.

# $1 is the name/prefix of the input psiblast and wordcount results
# $2 is the original input .sr file
# $3 is the minimum coverage for a good blast hit
# $4 is the minimum bit score for a good blast hit
# $5 is the minumum number of taxa in a family for further processing
# $6 is the minumum number of sequences in a family for further processing
#
# NOTE: $3 and $4 should be equal to $6 and $7 used in mmseq2blast.sh.
#
# This script takes an input .sr file and the previously constructed mmseqs clusters, psiblasts the clusters and remaining
# singletons against a blast database of consensus sequences and singletons, and sorts the resulting non-exclusive families
# by the number of sequences and taxa contained. Exclusive families are then assembled through a "greedy" search where
# sequences appearing in larger clusters (either by sequence or tax count) are removed from smaller clusters. The output
# are directories of the resulting tax-sorted and sequence sorted .sr files and word count files $1_taxSortFam, $1_seqSortFam,
# $1_MasterWC, $1_MatchWC, and $1TCWC/$1SCWC. Finally, it returns .cls files for tax sorted and seq sorted families, $1TC.cls/$1SC.cls.
# Note the last line in the .cls files lists all of the remaining sequences which were not clustered

Given the .cls file of choice, tax-sorted here, for a given host group and "binaryInteract", return to MATLAB to run processFamilies.

Tried the following:

./mmseq2blast.sh 'plant' 'plantSTRINGsr.txt' 0.5 '' 0.5 0.75 0.5
./blast2cluster.sh 'plant' 'plantSTRINGsr.txt' 0.75 0.5 2 0


Trimmed plant graph:


















onvert from .sr to .fa files and run mmseqs through wrap-around run_mmclust with -s=0.5 to get clusters.

e.g.

sr2fa virusSTRINGsr.txt -w=0 > tmp.fa
run_mmclust tmp.fa -s=0.5 -w=0 > virusMMC.txt

Move the *MMC.txt files to Desktop and delete. Run block 6 to get a plot of #proteins/#clusters.
Now run the mmseq2blast file. This script takes an input .sr file, clusters using mmseqs, aligns the clusters,
and calculates consensus sequences.

# $1 is the name/prefix of the output files
# $2 is the input .sr file
# $3 is the input similarity threshold for mmseqs
# $4 is whether to send to sge cluster during cls2ali, if no enter '', if yes enter '-sge'
# $5 is -grcut for making the consensus sequences
# $6 is the minimum coverage for a good blast hit
# $7 is a minimum bit score for a good blast hit

After this, run the and families are constructed in the $1PB.txt file.
Statistics for each family are recorded in $1WC.txt with the input query (MSA or singleton), number of taxa, and number
of sequences for the family.

psiblast -db STRING -query blastSTRING.fa -out blastSTRINGout.txt -outfmt 6 -comp_based_stats 0 
psiblast -db hostSTRING -query blastHostSTRING.fa -out blastHostSTRINGout.txt -outfmt 6 -comp_based_stats 0 
psiblast -db virusSTRING -query blastVirusSTRING.fa -out blastVirusSTRINGout.txt -outfmt 6 -comp_based_stats 0 

Given the blast results, run the first block of processBlastFamilies to generate the 'familyCell' variable which
has 6 fields virus taxID, virus tax.gene, host taxID, host tax.gene, virus family, host family. The "families" are
cell arrays. For each virus and host gene, the family cell contains a vector with indicies referencing the genes
(by index of the former fields) that were returned in the blast search given that gene as input.

Then run block 2 to incorporate the binaryInteract variable into family cell. The result is  a seventh field
with each row representing a binary interaction between the host indexed in column 1 and the virus indexed
in column 2.







