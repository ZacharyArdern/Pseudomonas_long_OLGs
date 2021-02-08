#!/bin/bash

#these files required: genomes-all_metadata.tsv genome_metadata.tsv

ORF_num="$1" ;

cd ../MAGs/ ;

data="../../Data/fna/" ; 

cat genomes-all_metadata.tsv | grep Pseudomonadales | awk '{print $NF}' | while read -r file_online ;
do 
genome=$(echo $file_online | awk -F "/" '{split($NF,a,"."); print a[1]}'  ) ;
if [[ -e "$data""$genome".gff.gz &&  ! -e "$data""$genome".fna ]] ; then wget $file_online -P $data ; 
unpigz -p 4 "$data"$genome.gff.gz ; 
cat "$data""$genome".gff | awk '/##FASTA/,0' | tail -n +2 > "$data""$genome".fna ;
pigz -p 4 "$data""$genome".gff ; 
fi ;
echo $genome ; 
done > Mgnify_genomes.txt ; 


# MAGs - GEM (JGI) and Mgnify/UHGG:

cat Mgnify_genomes.txt genome_metadata.tsv | awk '{if ($1~"GUT_GENOME" || $0~"Pseudomonadales") print $1}'  | 
awk '{print $1 "\t" NR}' |
while read -r genome number ; do

input="$genome".fna ; echo $genome $number ; 

test -e "$data"$input.gz && unpigz -p 4 "$data"$input.gz ; 

# 6 frame translation (to create database, to check seqs are homologous)
    test ! -e ${input%.*}_6frame.fasta &&  { for frame in 1 2 3 -1 -2 -3 ; do 
    cat "$data"$input | seqkit translate -f $frame | linear | 
    awk '{if ($1~">") print  ">" "'$genome'" "__" substr($1,2) "__"  "'$frame'" "\t" substr($1,2); else print}' ; 
    done ; } > ${input%.*}_6frame.fasta ;

    # Diamond BLAST

    diamond makedb --in ${input%.*}_6frame.fasta --db ${input%.*}_6frame.fasta &> /dev/null ; 

    diamond blastp --more-sensitive -q ../blastp/mGene"$ORF_num"_aa.fasta \
    -d ${input%.*}_6frame.fasta -o ${input%.*}_diamond.txt -p 4 -k 1 -e 1e-10 -b 4 -t ~/tmp/ \
    -f 6 pident ppos nident sseqid qseqid evalue slen sstart send qlen qstart qend stitle &> /dev/null ;

    # extract hit sequences
   
    # Note: shifted to bedtools getfasta and updated lengths calculation, for speed

    if [ -s ${input%.*}_diamond.txt ] ; then 
    
    echo "$input has homologue" ;

    cat "$data"$input | linear | paste - - | 
    awk -F "\t" '{split($1,a," "); print substr(a[1],2) "\t" length($2)}' > ${input%.*}_lengths.txt ;  


    cat ${input%.*}_diamond.txt | awk '{print $NF "\t" $4 "\t" $8 "\t" $9 "\t" $NF "\t" NR}' | 
    awk -F"\t" 'NR==FNR{a[$1]=$0;next}{ if (a[$1])print $0 "\t" a[$1]}' ${input%.*}_lengths.txt -  |
    while read stitle sseqid sstart send chr NR null lgth ; do 

    frame=$(echo $sseqid | awk -F "__" '{print $NF}' ) ; 
    if [ $frame == +1 ]; 
    then echo -e $chr"\t"$[$sstart*3-2]"\t"$[$send*3-3]"\t"$NR"\t""0""\t""+" ; 
    elif [ $frame == +2 ]; 
    then echo -e $chr"\t"$[$sstart*3-1]"\t"$[$send*3-2]"\t"$NR"\t""0""\t""+" ; 
    elif [ $frame == +3 ]; 
    then echo -e $chr"\t"$[$sstart*3]"\t"$[$send*3-1]"\t"$NR"\t""0""\t""+" ; 
    elif [ $frame == -1 ]; 
    then echo -e $chr"\t"$[$lgth-($send*3)+1]"\t"$[$lgth-($sstart*3)]"\t"$NR"\t""0""\t""-" ;  
    elif [ $frame == -2 ]; 
    then echo -e $chr"\t"$[$lgth-($send*3)+0]"\t"$[$lgth-($sstart*3)-1]"\t"$NR"\t""0""\t""-" ;
    elif [ $frame == -3 ]; 
    then echo -e $chr"\t"$[$lgth-($send*3)-1]"\t"$[$lgth-($sstart*3)-2]"\t"$NR"\t""0""\t""-" ;
    fi | 
    awk -F "\t" '{print $1 "\t" $2-1 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' |
    #> ${input%.*}_homologs.bed ;
    bedtools getfasta -s -fi "$data"$input -bed - | 
    awk -F "\t" '{split("'$sseqid'",a,"_"); if ($1~">") print ">"a[1]"_"a[2]":""'$sstart'""-""'$send'"; else print}' ; 
    done |
    tee mGene"$ORF_num"_"${input%.*}"_homologs_ntd.fasta |
    seqkit translate | linear > mGene"$ORF_num"_"${input%.*}"_homologs_aa.fasta ;

    fi ; 
done ; 

cat mGene"$ORF_num"*_homologs_ntd.fasta > ../blastp/metagenome_homologs_mGene"$ORF_num".fasta ;

