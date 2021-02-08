#!/bin/bash

#FUNCTIONS:
linear () { awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }'; } ;

gapless () { tee tmp_input.txt | grep ">" > headers.txt ; 
cat tmp_input.txt | grep -v ">" |  awk '{if ($1!~">") gsub(/.{3}/,"& ")}1' |  
 awk '{ for (i=1; i<=NF; i++)  {a[NR,i] = $i}}
NF>p { p = NF } END {for(j=1; j<=p; j++) {str=a[1,j]
        for(i=2; i<=NR; i++){str=str" "a[i,j];} print str }}' |
awk '{if ($0~/[A-Z]/) print}' | 
 awk '{ for (i=1; i<=NF; i++)  {a[NR,i] = $i}}
NF>p { p = NF } END {for(j=1; j<=p; j++) {str=a[1,j]
        for(i=2; i<=NR; i++){str=str" "a[i,j];} print str }}' | sed -e "s| ||g" > tmp.txt ;
paste headers.txt tmp.txt | tr "\t" "\n" ; } ;


#########

for ORF_num in 1 2 ; do ; 


#reference=AE004091.2 ;
reference=NC_002516.2; 


if [ $ORF_num == 1 ] ; then 
OLG_pos1=793 ; OLG_pos2=1749 ; 
simulation_outgroup_chr="LT629762.1" ;
elif [ $ORF_num == 2 ] ; then 
OLG_pos1=265 ; OLG_pos2=1800 ; 
simulation_outgroup_chr="3300013718_1" ;
fi

# files required: mGene"$ORF_num"_ntd.fasta ;
cat mGene"$ORF_num"_ntd.fasta | seqkit translate | linear > mGene"$ORF_num"_aa.fasta ;


#### Find any-frame amino acid homologs of OLG with tBLASTn (searching a translated nucleotide database)
mkdir tblastn ; 
cd tblastn ; 
for ORF_num in 1 2 ; do 
time tblastn -db nt -query ../blastp/mGene"$ORF_num"_aa.fasta -out OLG"$ORF_num"_tblastn.txt \
-evalue 1e-10 -max_hsps 1 -max_target_seqs 10000 -remote -entrez_query "Pseudomonadales [Organism]" \
-outfmt '6 qseqid sseqid pident sseq qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore ppos \
salltitles stitle sblastnames sacc saccver'; done ;


#### FIND HOMOLOGS OF MOTHER GENE WITH PROTEIN BLAST and IPG

# protein BLAST of mother gene against all genomes in the order Pseudomonadales (takes ~10mins)
time blastp -db nr -query mGene"$ORF_num"_aa.fasta -out OLG"$ORF_num"_blastp.txt \
-evalue 1e-10 -max_hsps 1 -max_target_seqs 10000 -remote -entrez_query "Pseudomonadales [Organism]" \
 -outfmt '6 qseqid sseqid pident sseq qlen slen length mismatch gapopen qstart qend sstart send evalue bitscore ppos \
  salltitles stitle sblastnames sacc saccver';


cat OLG"$ORF_num"_blastp.txt | awk -F "\t" '{split($2,a,"|"); print a[2]}' > accessions_OLG"$ORF_num".txt ; 

epost -input accessions_OLG"$ORF_num".txt -db protein | efetch -format ipg | 
awk '{if ($1!="Id") print}' > ipg_OLG"$ORF_num".txt ; 

echo "Downloading nucleotide sequences for BLASTp hits" ; 
time cat ipg_OLG"$ORF_num".txt | 
while read -r id db accession first last strand prot_acc names ; 
do #echo $accession $strand ; 
    if [ "$strand" == "+" ] ; then esearch -db nuccore -query $accession 2>/dev/null </dev/null | 
    efetch -seq_start $first -seq_stop $last -format fasta | linear ; 
    elif [ "$strand" == "-" ] ; then esearch -db nuccore -query $accession 2>/dev/null </dev/null | 
    efetch -seq_start $first -seq_stop $last -format fasta | linear | seqkit seq -p -r | linear ; 
    fi ; 
done > mGene"$ORF_num"_all_homologs_ipg_ntd.fasta ;

# Unique sequences, including reference seq header

cat mGene"$ORF_num"_all_homologs_ipg_ntd.fasta | linear | paste - - | tee tmp.txt |
awk -F "\t" '{print $2}' | sort | uniq | while read -r seq ; 
do cat tmp.txt | awk -F "\t" '{if ($2=="'$seq'") print $0}' | 
awk -F "\t" '$0 ~ /'$reference'/ {print} $0 !~ /'$reference'/ { a[f++] = $0 } 
END { for (i = 0; i < f; i++) { print a[i] } }' |
awk -F "\t" '{if (NR==1) print $1 "\n" $2}' ; done > mGene"$ORF_num"_homologs_ipg_ntd_ds.fasta ;


#### Get correct [forward] frame for ipg sequences (often incorrect due to contig endings)
input=mGene"$ORF_num"_homologs_ipg_ntd_ds.fasta ;
    test ! -s ${input%.*}_6frame.fasta &&  { 
    cat $input | seqkit translate -f 1  | linear | awk '{if ($1~">") print $1"__+1" "\t" substr($1,2); else print}' ;
    cat $input | seqkit translate -f 2  | linear | awk '{if ($1~">") print $1"__+2" "\t" substr($1,2); else print}' ;
    cat $input | seqkit translate -f 3  | linear | awk '{if ($1~">") print $1"__+3" "\t" substr($1,2); else print}' ;
    cat $input | seqkit translate -f -1 | linear | awk '{if ($1~">") print $1"__-1" "\t" substr($1,2); else print}' ;
    cat $input | seqkit translate -f -2 | linear | awk '{if ($1~">") print $1"__-2" "\t" substr($1,2); else print}' ;
    cat $input | seqkit translate -f -3 | linear | awk '{if ($1~">") print $1"__-3" "\t" substr($1,2); else print}' ;
    } > ${input%.*}_6frame.fasta ;
    # Diamond BLAST
    diamond makedb --in ${input%.*}_6frame.fasta --db ${input%.*}_6frame.fasta &> /dev/null ; 
    diamond blastp --more-sensitive -q mGene"$ORF_num"_aa.fasta \
    -d ${input%.*}_6frame.fasta -o ${input%.*}_"$ORF_num"_diamond.txt -p 4 -k 50000 -e 1e-10 -b 4 -t ~/tmp/ \
    -f 6 pident ppos nident sseqid qseqid evalue slen sstart send qlen qstart qend stitle &> /dev/null ;

    # extract hit sequences
    if [ -s ${input%.*}_"$ORF_num"_diamond.txt ] ; then 
    cat $input | linear | paste - - | 
    awk -F "\t" '{split($1,a," "); print substr(a[1],2) "\t" length($2)}' > ${input%.*}_lengths.txt ;  
    cat ${input%.*}_"$ORF_num"_diamond.txt | awk '{print $NF "\t" $4 "\t" $8 "\t" $9 "\t" $NF "\t" NR}' | 
    awk -F"\t" 'NR==FNR{a[$1]=$0;next}{ if (a[$1])print $0 "\t" a[$1]}' ${input%.*}_lengths.txt -  |
    while read stitle sseqid sstart send chr NR null lgth ; do 
    frame=$(echo $sseqid | awk -F "__" '{print $NF}' ) ; 
    if [ $frame == +1 ];   then echo -e $chr"\t"$[$sstart*3-2]"\t"$[$send*3-3]"\t"$NR"\t""0""\t""+" ; 
    elif [ $frame == +2 ]; then echo -e $chr"\t"$[$sstart*3-1]"\t"$[$send*3-2]"\t"$NR"\t""0""\t""+" ; 
    elif [ $frame == +3 ]; then echo -e $chr"\t"$[$sstart*3]"\t"$[$send*3-1]"\t"$NR"\t""0""\t""+" ; 
    fi ; done | awk -F "\t" '{print $1 "\t" $2-1 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' > ${input%.*}_homologs.bed ;
    fi ;     
bedtools getfasta -s -fi $input \
-bed ${input%.*}_homologs.bed > mGene"$ORF_num"_homologs_ipg_ntd_corrected.fasta ;


#### ADD IN HOMOLOGS FROM MAGS (metagenome assembled genomes)

bash metagenome_seqs.bash $ORF_num ; 
# > metagenome_homologs_mGene"$ORF_num".fasta

######################################################

#### PRODUCE TREE AND ALIGNMENT 

# Downsample sequences to unique mother gene sequences, making sure to include reference sequence header
# include reference sequence, otherwise take first seq among any identical ones 
# Note: could add warning somewhere if reference seq not present 
# only take sequences which are at least 50% length of the reference mGene sequence
reference_length=$(cat mGene"$ORF_num"_ntd.fasta | linear | awk '{if (NR==2) print length}') ;
cat mGene"$ORF_num"_homologs_ipg_ntd_corrected.fasta metagenome_homologs_mGene"$ORF_num".fasta | 
linear | paste - - | tee tmp.txt |
awk -F "\t" '{print $2}' | sort | uniq | while read -r seq ; 
do cat tmp.txt | awk -F "\t" '{if ($2=="'$seq'") print $0}' | 
awk -F "\t" '$0 ~ /'$reference'/ {print} $0 !~ /'$reference'/ { a[f++] = $0 } 
END { for (i = 0; i < f; i++) { print a[i] } }' |
awk -F "\t" '{if (NR==1) print $1 "\t" $2}' ; done | 
awk -F "\t" '{if (length($2)>=(0.5*"'$reference_length'")) 
print $1 "\n" $2}' > mGene"$ORF_num"_all-homologs_ntd_ds.fasta ;


# filter out sequences with ambiguous nucleotides and remove stop site at end

cat mGene"$ORF_num"_all-homologs_ntd_ds.fasta | paste - - |
awk -F "\t" '{if ($2!~"S" && $2!~"Y" && $2!~"M" && $2!~"R" && $2!~"K" && $2!~"W" && $2!~"N" ) 
print $1 "\n" substr($2,1,(length($2)-3))}' > mGene"$ORF_num"_all-homologs_ntd_ds2.fasta ;


# Align protein sequences and match nucleotides (Pal2Nal)

cat mGene"$ORF_num"_all-homologs_ntd_ds2.fasta | seqkit translate | linear | 
sed -e "s|*|X|g" > mGene"$ORF_num"_all-homologs_aa_ds2.fasta ; 

quickprobs-2.06-linux mGene"$ORF_num"_all-homologs_aa_ds2.fasta -o mGene"$ORF_num"_all-homologs_aa_ds2.aln ; 
perl ~/pal2nal.v14/pal2nal.pl mGene"$ORF_num"_all-homologs_aa_ds2.aln mGene"$ORF_num"_all-homologs_ntd_ds2.fasta \
 -output fasta | linear > mGene"$ORF_num"_all-homologs_pal2nal.fasta ;

# remove stop codon and for OLG1 exclude "3300012274" sequence which has early stop (not present in OLG2 alignment)

#cat mGene"$ORF_num"_all-homologs_pal2nal.fasta | 
#paste - - | grep -v 3300012274 | awk -F "\t" '{print $1 "\n" $2}' | 
#awk -F "\t" '{if ($1~">") print; 
#else print substr($1,1,(length($1)-3))}' > mGene"$ORF_num"_all-homologs_pal2nal2.fasta ; 

#iqtree -s mGene"$ORF_num"_all-homologs_pal2nal2.fasta -st CODON -asr -redo -bb 1000 ;


iqtree -s mGene"$ORF_num"_all-homologs_pal2nal.fasta -st NT2AA -redo -bb 1000 ;


#rename fasta 

cat mGene"$ORF_num"_all-homologs_pal2nal.fasta | awk -F "\t" '{if ($1~">") print}' | 
awk -F "\t" '{ print $1 ; gsub("[:+()]","_",$1); print $1}' | paste - - > mGene"$ORF_num"_renamefasta.txt ; 

cat mGene"$ORF_num"_renamefasta.txt | while read -r header renamed ; 
do 
cat mGene"$ORF_num"_all-homologs_pal2nal.fasta | paste - - |
awk -F "\t" '{if ($1=="'$header'") 
print "'$renamed'" "\n" $2}' ; done > mGene"$ORF_num"_all-homologs_pal2nal2.fasta ;
 
# trim to 20 representative genomes, including reference seq and outgroup seq for evolutionary simulation

cat mGene"$ORF_num"_all-homologs_pal2nal2.fasta | grep ">" | sed -e "s|>||g" |
sed -e "s|[(+:)]|_|g" | 
awk -F "\t" '{if ($1~"'$reference'") print $1 "," "reference";
else if ($1~"'$simulation_outgroup_chr'") print $1 "," "outgroup" }' > ref_info.txt ; 
cat ref_info.txt | awk -F "," '{print $2",1"}' > ref_info2.txt ;

python ~/Treemmer/Treemmer_v0.3.py -lm ref_info.txt -lmc ref_info2.txt \
mGene"$ORF_num"_all-homologs_pal2nal.fasta.contree -X 20 ;




#cat mGene"$ORF_num"_all-homologs_pal2nal.fasta.contree_trimmed_list_X_20 | awk -F "-" '{print $1}' | 
#awk '{if ($1~"[A-Z]") print}' | awk -F "_" '{if ($1~"\\.") print $1 ; else print $1 "_" $2}' |
#while read -r query ; do
#esearch -db nuccore -query "$query"  < /dev/null | elink -target taxonomy | efetch -format xml |
#xtract -pattern Taxon -element ScientificName ;
#done >  mGene"$ORF_num"_trimmed_list1.txt ; 
#add in MAGs

#species names
cat mGene"$ORF_num"_all-homologs_pal2nal.fasta.contree_trimmed_list_X_20 | awk -F "-" '{print $1}' |
awk -F "_" '{if ($1~"\\.") print $1 ; else print $1 "_" $2}' |
while read -r query ; do

if   [[ "$query" =~ [A-Z] ]] && [[ "$query" != *"GUT_GENOME"* ]]   ; then
esearch -db nuccore -query "$query"  < /dev/null | elink -target taxonomy | efetch -format xml |
xtract -pattern Taxon -element ScientificName | awk -F "\t" '{print $0 "\t" "'$query'"}' ;

elif [[ "$query" =~ [0-9] ]] && [[ "$query" != *"GUT_GENOME"* ]]  ; then 
cat ../MAGs/genome_metadata.tsv | grep $query | 
awk -F "\t" '{split($15,a,";"); 
if (a[7]!="s__") print a[7]; else print a[6]}' | awk -F "\t" '{print $0 "\t" "'$query'"}' ;

elif [[ "$query" =~ "GUT_GENOME" ]] ; then 
cat ../MAGs/genomes-all_metadata.tsv | grep $query | 
awk -F "\t" '{split($19,a,";"); print a[7]}' | awk -F "\t" '{print $0 "\t" "'$query'"}' ;

fi ; 
done > mGene"$ORF_num"_trimmed_list_X_20_species.txt ;  

cat mGene"$ORF_num"_trimmed_list_X_20_species.txt | 
awk -F "\t" '{split($1,a," "); if (a[1]=="Candidatus" || a[2]=="sp.") print a[1] "_" a[2] "_" a[3] ".."$2; 
else print a[1] "_" a[2] ".."$2}' | 
sed -e "s|s__||g" | sed -e "s|g__||g"  > mGene"$ORF_num"_trimmed_list_X_20_species2.txt ;

paste mGene"$ORF_num"_all-homologs_pal2nal.fasta.contree_trimmed_list_X_20 \
mGene"$ORF_num"_trimmed_list_X_20_species2.txt > rename.txt ;  



# extract Treemmer seqs from alignment

cat mGene"$ORF_num"_all-homologs_pal2nal.fasta.contree_trimmed_list_X_20 | awk '{print}' |
#awk -F "-" '{print $1}' | awk -F "_" '{if ($1~"\\.") print $1 ; else print $1 "_" $2}' | 
while read -r query ;
do 
cat mGene"$ORF_num"_all-homologs_pal2nal2.fasta | grep -A1 $query ; done > mGene"$ORF_num"_treemmer.fasta ;  

### add in out-of-Pseudomonadales outgroup to tree, to root

if [ $ORF_num == "1" ] ; then prot_query=WP_105240680.1 ; 
elif [ $ORF_num == "2" ] ; then prot_query=WP_124026310.1 ; 
fi ;

esearch -db protein -query "$prot_query"  < /dev/null | 
efetch -format fasta | linear > mGene_"$ORF_num"_outgroup_prot.fasta ;  

cat mGene"$ORF_num"_treemmer.fasta | seqkit translate | linear | paste - - |
awk -F "\t" '{gsub("-","",$2); print $1 "\n" $2}' > tmp.fasta ; 
cat mGene_"$ORF_num"_outgroup_prot.fasta tmp.fasta > tmp2.fasta ; 

quickprobs-2.06-linux tmp2.fasta -o mGene_"$ORF_num"_outgroup_prot2.fasta ; 
iqtree -s mGene_"$ORF_num"_outgroup_prot2.fasta -redo -bb 1000 ;

#reroot and then prune outgroup off 
nw_reroot mGene_"$ORF_num"_outgroup_prot2.fasta.contree $prot_query > mGene_"$ORF_num"_outgroup_prot2.fasta_rw.contree ;
nw_prune mGene_"$ORF_num"_outgroup_prot2.fasta_rw.contree $prot_query > mGene_"$ORF_num"_outgroup_prot2.fasta.contree ;

nw_rename mGene_"$ORF_num"_outgroup_prot2.fasta.contree rename.txt | 
sed -e "s|Pseudomonas|P.|g"  > renamed_"$ORF_num".tree ;



# 1 rename fasta file, 2 extract OLG region & frame, 3 translate and remove gaps found in all columns 
#1
cat rename.txt | 
#awk -F "\t" '{split($1,a,"-"); print a[1] "\t" $2}' | awk -F "\t" '{split($1,a,"_"); if (a[1]~"\\.") print a[1] "\t" $2; else print a[1] "_" a[2] "\t" $2}' | 
while read -r query new_name ; 
do 
cat mGene"$ORF_num"_treemmer.fasta | paste - - |
awk -F "\t" '{if ($1~"'$query'") print ">""'$new_name'" "\n" $2}' ; done > renamed.fasta ;
#2
adjusted_positions=$(cat renamed.fasta | grep -A1 $reference | head -2 | 
grep -v ">" | fold -w 1 | awk '{print $1 "\t" NR}' | awk '{if ($1!="-") print}' | 
awk '{if (NR==("'$OLG_pos1'"+0) || NR==("'$OLG_pos2'"+0)) print $2 }' | 
tr '\n' ' ' | awk '{print}')   ;
adj=($adjusted_positions) ;
cat renamed.fasta | awk -F "\t" '{if ($1~">") print; 
else print substr($1,"'${adj[0]}'",("'${adj[1]}'"-"'${adj[0]}'"+1))}' | 
#3
seqkit seq -p -r | seqkit translate | linear | gapless | sed -e "s|*|X|g" |
sed -e "s|Pseudomonas|P.|g" > renamed_OLGregion_"$ORF_num".fasta ; 


python3 ../scripts_long_OLGs/visualise_tree_aln.py renamed_OLGregion_"$ORF_num".fasta renamed_"$ORF_num".tree ;



done ;

