# hostVirus
## This is a step-by-step guide for generating the hostVirus coevolution results. Original code written for this project is contained in the repository but preexisting NCBI or Koonin group resources are not.

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

Now run Block 4 of [writeTaxLookup.m](writeTaxLookup.m). Given "uniprotSSr.txt", "subSTRING.sr", "uniprotPROsr.tab", and "uniprotWebscrapResults.txt" create a master .sr file by processing and matching the id's in "subSTRING.sr" to those in "subSTRINGnames.txt". Additionally, match the id's in "uniprotPROsr.tab" to those in "uniprotWebscrapResults.txt". Then take the relevant portion of the sequence from "uniprotPROsr.tab" and pair it with the query id in "uniprotWebscrapResults.txt". These id's are not unique and it's possible that some of the id's correspond to more than one taxID. This is true for the "subSTRINGnames" but does not appear to be true for those id's returned from uniprot. The third part of writeTaxLookup takes the first sequence for each repeated sequence (alphabetical order, but again does not seem to be applicable here) so that each id corresponds to a single sequence, removes any sequences which have id's matched to more than 1 taxID (which also did not happen originally) returns the combined .sr file including the orignal sequences from the STRING database, and returns the list of id's for which no sequence was found. This left 28 genes "lost", from 7 viruses three of which have other "found" genes with interactions and 4 which are totally absent from the final curated database [(in the original run)](lostSTRING.txt).

Now run block 5 of [writeTaxLookup.m](writeTaxLookup.m) to split the .sr files into host and virus [files](STRINGseq.zip).

## Cluster host sequences

Proceed to run the linux script [mmseq2blast.sh](mmseq2blast.sh) for each host group. This script takes an input .sr file, clusters using mmseqs, aligns the clusters, and calculates consensus sequences. After this, psiblast is run and families are constructed in the {prefix}PB.txt file. For the original run, $3, the input similarity threshold for mmseqs, was set to 0.5.

Given the mmseqs results, run the linux script [blast2cluster.sh](blast2cluster.sh) for various threshold choices. For the original run 0.75/0.5, 0.5/0.3, and 0/0 were tried for $3, the minimum coverage for a good blast hit, and $4, the minimum bit score for a good blast hit. $5, the minumum number of taxa in a family for further processing, and $6, the minumum number of sequences in a family for further processing, were both kept at 0. This script takes the input $1PB.txt file and calculates stats after the $3/4 thresholds are applied. These non-exclusive families are then sorted by the number of taxa contained. Exclusive families are then assembled through a "greedy" search where sequences appearing in larger clusters (by tax count not sequence count) are removed from smaller clusters. The output are directories of the resulting tax-sorted .sr files and word count files ...taxSortFam, ...MasterWC, and $1TCWC. Finally, it returns .cls files for tax sorted families, $1TC.cls.

Families change composition with decreasing thresholds and in general become larger and less discriminant, the desired families (see next section) are insensitive to changes in these parameters and 0.75/0.5 were used in the original run. Note, originally two different sortings for the next step, number of species in a family and the number of sequences in a family, were tried. Only tax-sorted results which in preliminary evaluation did not differ dramatically from seq-sorted results at reasonable threshold values were kept in the final version of the script; however, this study was conducted on the plant dataset and may differ in the vertebrate dataset.

Given the .cls file for each host group and "binaryInteract", return to MATLAB to run processFamilies.
It takes the above families and makes a shortlist including only those that
interact with at least one viral protein which in turn only interacts with that family (or below some other
number of families if desired). It additionally thresholds by number of such interactions for the host family.
Three levels of evidence 1) single host protein, viral protein. 2) many to 1 viral protein to host protein or
host protein to viral protein. 3) multiple parallel host to viral protein interactions. This
script returns two output, a binary interaction graph displaying these results is printed and a tab document
with three fields 1) host family index, 2) host protein name, and 3) viral protein name is returned.

Then align the families with microCoAlign which produces two .sr files (host/virus) which are manually reviewed
and annotated for further consideration.

./microCoAlign.sh 'plantReviewInteractions.txt' 'vertebrateReviewInteractions.txt' 'hostSTRINGsr.txt' 'virusSTRINGsr.txt'

When there are many viruses in a group, the cls2ali output may be difficult to parse, so move to prof_align for a repeat
clustering and alignment step:

cut -f 1 virus5.sr > tmp_1.txt
cut -f 4 -d'|' tmp_1.txt > tmp_2.txt
sr2fa virusSTRINGsr.txt -w=0 > tmp_1.fa
makeblastdb -in tmp_1.fa -input_type fasta -dbtype prot -parse_seqids -out tmp_Blast
prof_align tmp_2.txt -db=tmp_Blast -nosge -w=0 -n=virus5
rm tmp_*

These results (about 10 preferred targets returned) exclude cases where there may be a larger virus protein family
with many interactions with a given host protein family but also interactions with other host proteins outside that family.
To retrieve any additional cases of interest, move on to cluster viral proteins into families.

./mmseq2blast.sh 'virus' 'virusSTRINGsr.txt' 0.5 ''
./blast2cluster.sh 'virus' 'virusSTRINGsr.txt' 0.75 0.5 0 0
./blast2cluster.sh 'virus' 'virusSTRINGsr.txt' 0.5 0.3 0 0
./blast2cluster.sh 'virus' 'virusSTRINGsr.txt' 0.3 0.1 0 0
./blast2cluster.sh 'virus' 'virusSTRINGsr.txt' 0 0 0 0

Coronavirus spike proteins were used as a litmus test and were not clustered until low thresholds of 0.3/0.1. These results
are much more sensitive to threshold choice than the host proteins. Given the host and viral clusters, proceed to run the
MATLAB script "processFamilyPairs2". This script identifies host protein - viral protein family pairs that have N parallel
interactions between them (i.e. one host protein to one viral protein) such that the viral protein family has at most N-2
parallel interactions with any other host family. The number of parallel interactions is determined by making a host-virus
protein-protein interaction matrix for proteins within the cluster pair (binary, 1 interacting, 0 not interacting) and finding
the rank of that matrix. These results are saved in a *ReviewFamilyPairs.txt file with fields N; pair index for a given N;
number of host families with m parallel interactions with the given virus family m running from 1 to max # observed
separated by pipes; host protein; virus protein.


./blastTargetNR.sh 'hostTargetCLS.txt' 'virusTargetCLS.txt' 'hostSTRINGsr.txt' 'virusSTRINGsr.txt'

./blastTargetNR2.sh 'STRINGtargets.txt' 'hostSTRINGsr.txt' 'virusSTRINGsr.txt'

./blastTargetNR2.sh 'TESTtargets.txt' 'hostSTRINGsr.txt' 'virusSTRINGsr.txt'

web blast for individual proteins or command line for -in_msa?

cat hostSTRINGsr.txt > tmp_sr.txt
cat virusSTRINGsr.txt >> tmp_sr.txt
cut -f 1 STRINGtargets.txt > tmp_list.txt
cut -f 2 STRINGtargets.txt >> tmp_list.txt

tab_select tmp_sr.txt -t=tmp_list.txt > tmp_sr2.txt
sr2fa tmp_sr2.txt > nrQueries.txt

rm tmp_*

split -l 2 -d nrQueries.txt nrQuery



nonstandard psiblast values: maxtargetsequences=10000, compositional adjustments=no adjustment



