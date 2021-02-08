#!/bin/bash

# install: 
#sudo apt-get install embassy-phylip

reference=NC_002516.2; 

for ORF_num in 1 2; do 

if [ $ORF_num == 1 ] ; then 
OLG_pos1=793 ; OLG_pos2=1749 ; 
elif [ $ORF_num == 2 ] ; then 
OLG_pos1=265 ; OLG_pos2=1800 ; 
fi ;

#FUNCTIONS:
linear () { awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }'; } ;

#######

# positions:

adjusted_positions=$(cat ../blastp/mGene"$ORF_num"_all-homologs_pal2nal.fasta | grep -A1 $reference | head -2 | 
grep -v ">" | fold -w 1 | awk '{print $1 "\t" NR}' | awk '{if ($1!="-") print}' | 
awk '{if (NR==("'$OLG_pos1'"+0) || NR==("'$OLG_pos2'"+0)) print $2 }' | 
tr '\n' ' ' | awk '{print}')   ;
adj=($adjusted_positions) ;


#OLG clade

if [ $ORF_num == 1 ] ; then 
cat ../blastp/mGene"$ORF_num"_all-homologs_pal2nal.fasta | paste - - | awk -F "\t" '{split($2,a,""); 
if (a["'${adj[0]}'"]=="T" && a["'${adj[0]}'"+1]=="C" && a["'${adj[0]}'"+2]=="A" \
&& a["'${adj[1]}'"-2]=="C" && a["'${adj[1]}'"-1]=="A" && a["'${adj[1]}'"]=="T") 
print $1 "\n" $2}' > alignment_OLG"$ORF_num".fasta ;  
elif [ $ORF_num == 2 ] ; then 
cat ../blastp/mGene"$ORF_num"_all-homologs_pal2nal.fasta | paste - - | awk -F "\t" '{split($2,a,""); 
if (a["'${adj[0]}'"]=="C" && a["'${adj[0]}'"+1]=="T" && a["'${adj[0]}'"+2]=="A" ) 
print $1 "\n" $2}' > alignment_OLG"$ORF_num".fasta ;  
fi ;


# mGene frame prot sequence, OLG clade

cat alignment_OLG"$ORF_num".fasta | paste - - | awk -F "\t" '{split($2,a,""); 
print $1 "\n" substr($2,"'${adj[0]}'",("'${adj[1]}'"-"'${adj[0]}'"+1))}' | 
seqkit translate | linear | sed -e "s|[()]||g" | awk -F ":" '{print $1}' > mGene_"$ORF_num"_OLGclade_prot.fasta ;  
fprotdist -filter -sequence mGene_"$ORF_num"_OLGclade_prot.fasta -outfile mGene_"$ORF_num"_OLGclade_prot_dist.txt ;  
cat mGene_"$ORF_num"_OLGclade_prot_dist.txt | tail -n +2 | 
awk '/NC_0025/{flag=1;} /[A]/{flag=0}flag' | tr ' ' '\n' | 
awk '{if ($1!="" && $1!~/[A-Z]/) print}' > mGene_"$ORF_num"_OLGclade_prot_dist2.txt ;

# OLG frame prot sequence, OLG clade

cat alignment_OLG"$ORF_num".fasta | paste - - | awk -F "\t" '{split($2,a,""); 
print $1 "\n" substr($2,"'${adj[0]}'",("'${adj[1]}'"-"'${adj[0]}'"+1))}' | seqkit seq -p -r | 
seqkit translate | linear | sed -e "s|[()]||g" | awk -F ":" '{print $1}' > OLG_"$ORF_num"_OLGclade_prot.fasta ;  
fprotdist -filter -sequence OLG_"$ORF_num"_OLGclade_prot.fasta -outfile OLG_"$ORF_num"_OLGclade_prot_dist.txt ;  
cat OLG_"$ORF_num"_OLGclade_prot_dist.txt | tail -n +2 | 
awk '/NC_0025/{flag=1;} /[A]/{flag=0}flag' | tr ' ' '\n' | 
awk '{if ($1!="" && $1!~/[A-Z]/) print}' > OLG_"$ORF_num"_OLGclade_prot_dist2.txt ;


done ; 

