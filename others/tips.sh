
# to select tophits from blast output file
awk '!x[$1]++' gene2.KB.blast > tophits.txt

# if the headers are repeated in a fasta file
perl -i.bak -pe 's/$/_$seen{$_}/ if ++$seen{$_}>1 and /^>/; ' file.fa

# tool to print idl wich have/haven't ????
awk -F "\t"  'FNR==NR{a[$1]=$2FS$3; next} ($1 in a){print $0,a[$1]; next} {print $0,"NA","NA"}' EUGGR.annotation gene.idl > awk.out 


# Create a database of uniprotK (uniprot and swissprot)
makeblastdb -in uniprot_sprot.fasta -input_type fasta -dbtype prot -out uniprot

# Make a blast protein between proteins found in compomics and uniprotdb created (downloaded from uniprot)
blastp -query proteins.fasta -db /media/javier/Datas/uniprotdb/uniprot -out uniprot_blastp -outfm 6

blastp -query FINAL.aa_pub.fa -db /media/javier/Datas/uniprotdb/uniprot -out FINAL_blastp -max_hsps 1 -max_target_seqs 1 -num_threads 16 -outfmt 6 # evalue E -10

# Create a file with the correspondance column betwen evigene pub id and uniprot id. Extracting those from blatp results obtained from uniprotK.
paste <(cut -f1 uniprot.blastp | sed 's/sp|//' | sed 's/tr|//' | sed 's/|.*//') <(cut -f2 uniprot.blastp | sed 's/sp|//' | sed 's/tr|//' | sed 's/|.*//') > uniprot.col

## ?? paste <(cut -f2 uniprot.blastp | sed 's/sp|//' | sed 's/tr|//' | sed 's/|/ /' | cut -d ' ' -f1 ) <(cut -f1 uniprot.blastp | sed 's/sp|//' | sed 's/tr|//' | sed 's/|/ /' | cut -d ' ' -f1 ) > uniprot.idf

### --- ###
# map tophits from blast with Retrieve/ID mapping tool of uniprot, in order to get GO terms and KO terms (RefSeq Nucleotide or UniprotKB AC/ID)

tail -n +2 /home/javier/Downloads/uniprot-yourlist%3mapped.tab | cut -f1,3,4 | sed 's/; /,/g' | sed 's/;//' > gene2.nt.tomap.ann

# skip first line, take columns of interest and create tomap.ann file



# Extract results from blastp and add the correspondance between evigene pub id 
awk '(NR==FNR){unp[$2]=$2;id[$2]=$1;next}($1 in unp){print id[$1]"\t"$0}' uniprot.col uniprot.tab > file4.txt
awk -F "\t" '(NR==FNR) {a[$1]=$2;next} ($1 in a) {print $0 "\t" a[$1]}' table list
awk -F "\t" '(NR==FNR) {a[$1]=$2FS$3FS$4;next} ($1 in a) {print $0 "\t" a[$1]}' table list

	NR==FNR # NR is the current input line number and FNR the current file's line number. The two will be equal only while the 1st file is being read.
	unp[$2]=$2 # Save column 2 (unp) from the first file (argument) uniprot.col in hash-array using column 2 as the key
	id[$2]=$1 # Save column 1 (id) from the first file (argument) uniprot.col in hash-array using column 2 as the key
	next # Then, skip to the next line so that this is only applied on the 1st file.
	($1 in unp) # check whether field 1 of second file uniprot.tab is in array unp.
	{print id[$1]"\t"$0} # If that's true print the id column from uniprot.col and the entire uniprot.tab with tab separator.



# Para excluir lineas (por que corresponden a contaminantes), necesitamos una lista de ID a excluir y la tabla general con IDs buenos y malos.
awk -F "\t" 'FNR==NR{a[$1]=$1;next}!($1 in a){print $0}' conta.idl proteins.rat > proteins_clean.rat




###
#awk -F "\t" 'FNR==NR{a[$1]=$2 FS $3 FS $4 FS $5;next}($1 in a){print $1 a[$1]}' proteins.rat file4.txt > file3.txt
#https://askubuntu.com/questions/707843/merge-files-using-a-common-column
#deberia cambiar el orden de los archivos? file4 es mas pequeno
###

##########################
### EVIGENE ANNOTATION ###
##########################

# Annotation of full proteins from final of evigene set.
blastp -query FINAL.aa_pub.fa -db /media/javier/Datas/uniprotdb/uniprot -out FINAL_50_blastp -max_hsps 1 -evalue 1e-50 -word_size 3 -max_target_seqs 1 -num_threads 16 -outfmt 6

##############################################
#### PROTEINAS QUANTIFICADAS CON REPORTER ####
##############################################

### 712 proteins condensed and cleaned with log ratio
cut -f 1 proteins_cond_clean_log.rat > FINAL.idl
### create equivalence identifiers betwen public ID evigene-Uniprot
awk -F "\t" 'FNR==NR{a[$1]=$1;next}($1 in a){print $0}' FINAL.idl uniprot.col > FINAL.col
### Uniprot Annotation of 
awk '(NR==FNR){unp[$2]=$2;id[$2]=$1;next}($1 in unp){print id[$1]"\t"$0}' FINAL.col uniprot.tab > FINAL.ann


### we need uniprot ID with a tag (prefix) type UNIPROT:
paste -d '\t' <(cut -f 1 FINAL.unp) <(cut -f 2 FINAL.col) <(cut -f 1 FINAL.col) > ipath
cat ipath | sed 's/:\t/:/' > ipatH
cut -f 1 ipatH > ipath
# validated list in website https://pathways.embl.de/tools.cgi
# resulting file iPath_ID_validation_50e4olNrnJ.txt

cat iPath_ID_validation_50e4olNrnJ.txt | sed 's/:/:\t/' > ipath.val
awk -F "\t" 'FNR==NR{a[$2]=$2;next}($2 in a){print $0}' ipath.val FINAL.col > ipath.col
awk -F "\t" 'FNR==NR{a[$1]=$1;next}($1 in a){print $0}' ipath.col proteins_cond_clean_log.rat > ipath.rat
awk -F "\t" 'FNR==NR{a[$1]=$1;b[$1]=$2;c[$1]=$3;d[$1]=$4;e[$1]=$5;next}($1 in a){print $0"\t"b[$1]"\t"c[$1]"\t"d[$1]"\t"e[$1]}' ipath.rat ipath.col > ipath.full


### R
### CONDENSE DUPLICATE ID CALCULATING MEAN ###

directory = '/media/javier/Datas/Proteomics/Euglena/ALLclean_MTCwt/new_assembly2'
setwd(directory)

library(limma)
library(edgeR)
library(ggplot2)
library(reshape2)
library(data.table)

## how to condense rows from a table that has the same identifier. We want to calculate the mean of the duplicate ID
## we invoke the library (from data.camp)

## we read our table
dat <- read.table(file='/media/javier/Datas/Proteomics/Euglena/ALLclean_MTCwt/new_assembly2/ipath.full', sep='\t', header=FALSE)
## we set our key term, such as ID (in this case)
dat <- dat[,2:6]
keys <- colnames(dat[1])
## we transform our data into a table
X <- as.data.table(dat)
## we calculate the mean of duplicate rows (!!!how it works??)
ratio <- X[,lapply(.SD,mean),keys]
## we transform the output into dataframe
ipath.dat <- as.data.frame(ratio)
ipath.dat$V1 <- c(rep("UNIPROT:", length(ipath.dat[,1])))
## we export our new table

write.table(ipath.dat, file='/media/javier/Datas/Proteomics/Euglena/ALLclean_MTCwt/new_assembly2/ipath.dat', quote=FALSE, row.names=FALSE,  sep='\t')

### UNIX

paste <(cut -f 6 ipath.dat) <(cut -f 1,2,3,4,5 ipath.dat) > ipath2.dat
cat ipath2.dat | sed 's/:\t/:/' > ipath.dat
