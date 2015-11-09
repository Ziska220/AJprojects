 #!/bin/bash 

file=rmstart-CC5dox.bed
name=CC5
genefile=refGene.bed
bpwindow=5000
hits=3


####Unique Alignemnts
#Keep only unique reads from aligned, cleaned reads.

awk '{print $1 "\t" $2 "\t" $4 "\t" $3}' $file > re-$file
sort re-$file | rev | uniq -f 1 | rev > re-unique-$file
awk '{print $1 "\t" $4 "\t" $3 "\t" $2}' re-unique-$file > re-re-unique-$file
sort re-re-unique-$file | rev | uniq -f 1 | rev > re-re-unique-uniq-$file 
awk '{print $1 "\t" $4 "\t" $2 "\t" $3}' re-re-unique-uniq-$file > $name-uniq.bed
echo "Number of alignments per file:"
wc -l $file
wc -l $name-uniq.bed
rm re-*

###Genes with 'hits' number of unique alignments

bedtools window -a $genefile -b $name-uniq.bed -c -w $bpwindow > $name$ext-gene-$bpwindow.bed
awk -v var="$hits" '$13>var' $name$ext-gene-$bpwindow.bed > $name$ext-gene-$bpwindow-$hits.bed 
awk '{print $1 "\t" $2 "\t" $3 "\t" $13}' $name$ext-gene-$bpwindow-$hits.bed > $name$ext-gene-$bpwindow-$hits-clean.bed
sort $name$ext-gene-$bpwindow-$hits-clean.bed | rev | uniq -f 2 | rev > $name$ext-gene-$bpwindow-$hits-uniq.bed
awk '{print $1 "\t" $3 "\t" $2 "\t" $4}' $name$ext-gene-$bpwindow-$hits-uniq.bed > re-$name$ext-gene-$bpwindow-$hits-uniq.bed
sort re-$name$ext-gene-$bpwindow-$hits-uniq.bed | rev | uniq -f 2 | rev > re-uniq-re-$name$ext-gene-$bpwindow-$hits-uniq.bed
awk '{print $1 "\t" $3 "\t" $2 "\t" $4}' re-uniq-re-$name$ext-gene-$bpwindow-$hits-uniq.bed > $name$ext-gene-$bpwindow-$hits-uniq.bed
rm re-*
wc -l $name$ext-gene-$bpwindow-$hits.bed
wc -l $name$ext-gene-$bpwindow-$hits-uniq.bed
rm $name$ext-gene-$bpwindow.bed 
rm $name$ext-gene-$bpwindow-$hits.bed
rm $name$ext-gene-$bpwindow-$hits-clean.bed
rm $name-uniq.bed

