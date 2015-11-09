#!/bin/bash 

date '+%m/%d/%Y_%H:%M:%S'

file=rmstart-CC5dox.bed
ctrl1=rmstart-CC5nodox.bed
ctrl2=rmstart-CC5PB.bed
name=SU2
genefile=refGene.bed
bpwindow=5000
hits=3
chang=hglft_genome_267b_3af6e0.bed
changwin=250000


####Unique Alignemnts
#Keep only unique reads from aligned, cleaned reads.

echo "Number of alignments per file:"

function unique_alignments () {
awk '{print $1 "\t" $2 "\t" $4 "\t" $3}' $1 > re-$1
sort re-$1 | rev | uniq -f 1 | rev > re-unique-$1
awk '{print $1 "\t" $4 "\t" $3 "\t" $2}' re-unique-$1 > re-re-unique-$1
sort re-re-unique-$1 | rev | uniq -f 1 | rev > re-re-unique-uniq-$1 
awk '{print $1 "\t" $4 "\t" $2 "\t" $3}' re-re-unique-uniq-$1 > $name-$1-uniq.bed
wc -l $1
wc -l $name-$1-uniq.bed
rm re-*
}

unique_alignments "$file" 
unique_alignments "$ctrl1"

###Subtract nodox from dox

bedtools intersect -a $name-$file-uniq.bed -b $name-$ctrl1-uniq.bed -c > $name-dox-ctrl1.bed
awk '$5<2' $name-dox-ctrl1.bed > $name-dox-nodox-less2.bed
awk '{print $1 "\t" $2 "\t" $3 "\t" $4}' $name-dox-nodox-less2.bed > $name-dox-sub-nodox.bed
rm *less2.bed
rm $name-dox-ctrl1.bed
rm $name-$file-uniq.bed
rm $name-$ctrl1-uniq.bed
wc -l $name-dox-sub-nodox.bed

### Genes with over "hits" number of alignments in Dox minus nodox

bedtools window -a $genefile -b $name-dox-sub-nodox.bed -c -w $bpwindow > $name-dox-sub-nodox-gene-$bpwindow.bed
awk -v var="$hits" '$13>var' $name-dox-sub-nodox-gene-$bpwindow.bed > $name-dox-sub-nodox-gene-$bpwindow-ovr$hits.bed
wc -l $name-dox-sub-nodox-gene-$bpwindow-ovr$hits.bed
rm $name-dox-sub-nodox-gene-$bpwindow.bed
awk '{print $1 "\t" $2 "\t" $3 "\t" $13}' $name-dox-sub-nodox-gene-$bpwindow-ovr$hits.bed > $name-dox-sub-nodox-gene-$bpwindow-ovr$hits-clean.bed
sort $name-dox-sub-nodox-gene-$bpwindow-ovr$hits-clean.bed | rev | uniq -f 2 | rev > re-$name-dox-sub-nodox-gene-$bpwindow-ovr$hits-uniq.bed
awk '{print $1 "\t" $3 "\t" $2 "\t" $4}' re-$name-dox-sub-nodox-gene-$bpwindow-ovr$hits-uniq.bed > re-re-$name-dox-sub-nodox-gene-$bpwindow-ovr$hits-uniq.bed
sort re-re-$name-dox-sub-nodox-gene-$bpwindow-ovr$hits-uniq.bed | rev | uniq -f 2 | rev > re-re-$name-dox-sub-nodox-gene-$bpwindow-ovr$hits-uniq-uniq.bed
awk '{print $1 "\t" $3 "\t" $2 "\t" $4}' re-re-$name-dox-sub-nodox-gene-$bpwindow-ovr$hits-uniq-uniq.bed > $name-dox-sub-nodox-gene-$bpwindow-ovr$hits-uniq.bed
rm $name-dox-sub-nodox-gene-$bpwindow-ovr$hits.bed
rm $name-dox-sub-nodox-gene-$bpwindow-ovr$hits-clean.bed
rm $name-dox-sub-nodox.bed
rm re-*
wc -l $name-dox-sub-nodox-gene-$bpwindow-ovr$hits-uniq.bed


### Genes with over 4 PB
#control2-PB

unique_alignments "$ctrl2"

bedtools window -a $genefile -b $name-$ctrl2-uniq.bed -c -w $bpwindow > $name-PB-gene-$bpwindow.bed
awk '$13>4' $name-PB-gene-$bpwindow.bed > $name-PB-gene-$bpwindow-ovr4.bed
wc -l $name-PB-gene-$bpwindow-ovr4.bed
rm $name-$ctrl2-uniq.bed
rm $name-PB-gene-$bpwindow.bed

### Remove gene with over 4 PB from Dox minus NoDox

bedtools window -a $name-dox-sub-nodox-gene-$bpwindow-ovr$hits-uniq.bed -b $name-PB-gene-$bpwindow-ovr4.bed -v > $name-dox-controls-ovr$hits.bed 
rm $name-PB-gene-$bpwindow-ovr4.bed
rm $name-dox-sub-nodox-gene-$bpwindow-ovr$hits-uniq.bed
wc -l $name-dox-controls-ovr$hits.bed


### overlap with chang data

bedtools window -a $name-dox-controls-ovr$hits.bed -b $chang -c -w $changwin > $name-dox-controls-ovr$hits-chang-$changwin.bed
rm $name-dox-controls-ovr$hits.bed
awk '$5>0' $name-dox-controls-ovr$hits-chang-$changwin.bed > $name-dox-controls-ovr$hits-chang-$changwin-ovr0.bed
rm $name-dox-controls-ovr$hits-chang-$changwin.bed
awk '{print $1 "\t" $3 "\t" $2 "\t" $4 "\t" $5}' $name-dox-controls-ovr$hits-chang-$changwin-ovr0.bed > re-$name-dox-controls-ovr$hits-chang-$changwin-ovr0.bed
sort re-$name-dox-controls-ovr$hits-chang-$changwin-ovr0.bed | rev | uniq -f 3 | rev > re-$name-dox-controls-ovr$hits-chang-$changwin.bed
awk '{print $1 "\t" $3 "\t" $2 "\t" $4 "\t" $5}' re-$name-dox-controls-ovr$hits-chang-$changwin.bed > unsort-$name-dox-controls-ovr$hits-chang-$changwin.bed 
sort unsort-$name-dox-controls-ovr$hits-chang-$changwin.bed > $name-dox-controls-ovr$hits-chang-$changwin.bed
rm unsort*
rm $name-dox-controls-ovr$hits-chang-$changwin-ovr0.bed
rm re-*
wc -l $name-dox-controls-ovr$hits-chang-$changwin.bed

awk '{print $1 "\t" $2 "\t" $3}' $name-dox-controls-ovr$hits-chang-$changwin.bed > $name-dox-controls-ovr$hits-chang-$changwin-cln.bed
sort $name-dox-controls-ovr$hits-chang-$changwin-cln.bed > $name-dox-controls-ovr$hits-chang-$changwin-cln-s.bed 
rm $name-dox-controls-ovr$hits-chang-$changwin-cln.bed

