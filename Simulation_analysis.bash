#!/bin/bash

# Written by Zachary Ardern (zachary.ardern@tum.de / z.ardern@gmail.com)
# Implements a simplified version of the method developed by Cassan et al. (2016) PNAS
# This script rewrites the work-flow, makes it compatible with multiple input genes, 
# and replaces CodonPhyML and Alfsim with IQTREE and Pyvolve respectively

seq_set_size=20 ; 
trees=10 ; 
simulations=10 ;
model_type=ECMunrest ;


omega_value=0.5 ;
#omega_value=0.3 ;
#omega_value=0.7 ;

#ORF_num=1 ; 
ORF_num=2 ; 

#codonfreqs=customfreqs ; 
codonfreqs=defaultfreqs ; 

reference=NC_002516.2;

if [ $ORF_num == 1 ] ; then 
OLG_pos1=793 ; OLG_pos2=1749 ; 
#outgroup_chromosome=LT629762.1 ;
outgroup_chromosome=NZ_CP008696.1 ; #alternative outgroup
elif [ $ORF_num == 2 ] ; then 
OLG_pos1=265 ; OLG_pos2=1800 ; 
#outgroup_chromosome="3300013718_1" ;
outgroup_chromosome="NZ_AUIE01000001.1" ; #alternative outgroup
fi ; 


# FUNCTIONS:

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


### prepare sequences of interest: 

full_aln_positions=$(cat ../blastp/mGene"$ORF_num"_all-homologs_pal2nal.fasta | grep -A1 $reference | head -2 | 
grep -v ">" | fold -w 1 | awk '{print $1 "\t" NR}' | awk '{if ($1!="-") print}' | 
awk '{if (NR==("'$OLG_pos1'"+0) || NR==("'$OLG_pos2'"+0)) print $2 }' | 
tr '\n' ' ' | awk '{print}')   ;
aln=($full_aln_positions) ;

if [ $ORF_num == 1 ] ; then 
cat ../blastp/mGene"$ORF_num"_all-homologs_pal2nal.fasta | paste - - | awk -F "\t" '{split($2,a,""); 
if (a["'${aln[0]}'"]=="T" && a["'${aln[0]}'"+1]=="C" && a["'${aln[0]}'"+2]=="A" \
&& a["'${aln[1]}'"-2]=="C" && a["'${aln[1]}'"-1]=="A" && a["'${aln[1]}'"]=="T") 
print $1 "\n" $2}' > alignment_OLG"$ORF_num".fasta ;  
elif [ $ORF_num == 2 ] ; then 
cat ../blastp/mGene"$ORF_num"_all-homologs_pal2nal.fasta | paste - - | awk -F "\t" '{split($2,a,""); 
if (a["'${aln[0]}'"]=="C" && a["'${aln[0]}'"+1]=="T" && a["'${aln[0]}'"+2]=="A" ) 
print $1 "\n" $2}' > alignment_OLG"$ORF_num".fasta ;  
fi ;

cat alignment_OLG"$ORF_num".fasta | paste - - > tmp.txt; 
cat ../blastp/mGene"$ORF_num"_all-homologs_pal2nal.fasta | paste - - |
awk -F"\t" 'NR==FNR{a[$1]=$1;next}{ if (!a[$1])print ;}' tmp.txt - |
awk -F "\t" '{print $1 "\n" $2}' > rest_of_alignment_OLG"$ORF_num".fasta ; 



### find outgroup seq from within OLG clade as seq with highest sum of pairwise distances with other sequences
#if [ $ORF_num == 2 ] ; then 
#cat alignment_OLG"$ORF_num".fasta | paste - - | 
#grep -v "NC_022360.1" | grep -v "NZ_AYON01000044.1" | grep -v "NZ_ASQU01000116.1" | 
#grep -v "NZ_LLNE01000012.1" | grep -v "NZ_JTSM01000116.1" | grep -v "NZ_LLVD01000012.1" | 
#grep -v "NZ_JTRY01000003.1" | grep -v "NZ_JTVG01000048.1" | grep -v "NZ_JTUN01000004.1" | 
#grep -v "NZ_ATAH01000207.1" | grep -v "NZ_LZQG01000064.1" | grep -v "NZ_LLNW01000001.1" | 
#grep -v "NZ_JTUG01000023.1" | grep -v "NZ_KI518928.1" | grep -v "NZ_JTUM01000003.1" | 
#grep -v "NZ_CP007224.1" | grep -v "NZ_LLLQ01000034.1" | grep -v "NZ_ASQY01000040.1" | 
#grep -v "NZ_JVSG01001464.1" | 
#awk -F "\t" '{print $1 "\n" $2}' > tmp.fasta ;  
#FastTree tmp.fasta > alignment_OLG"$ORF_num"_fast.tree ;
#nw_distance -n -m m alignment_OLG"$ORF_num"_fast.tree > alignment_OLG"$ORF_num"_fast_matrix.txt ; 
#outgroup_chromosome=$(cat alignment_OLG"$ORF_num"_fast_matrix.txt | tail -n +2 | 
#awk '{for (i=2;i<=NF;i++) sum[i]+=$i; }; END{for (i in sum) print sum[i];}' > tmp.txt ;
#cat alignment_OLG"$ORF_num"_fast_matrix.txt | awk '{print $1}' | tail -n +2 > tmp2.txt ; 
#paste tmp2.txt tmp.txt | sort -rgk2,2 | head -1 | awk '{print $1}' ) ; 
#fi ; 

##############################################

# for each sample, and simulations, do:

cat ../blastp/mGene"$ORF_num"_all-homologs_pal2nal.fasta | 
grep -A1 "$outgroup_chromosome" | head -2 > outgroup_seq.fasta ;


seq $trees | while read -r i ; 
do 


# pick some sequences from full alignment, in 'intact OLG region'
# match full alignment with list of seqs in 'intact OLG region' and pick random seqs, of number $seq_set_size  
# ensure only unique chrms picked 
# remove gap-only columns

cat alignment_OLG"$ORF_num".fasta | paste - - | 
awk -F "\t" '{split($1,a,":"); print a[1]}' | sort | uniq | while read -r chr ; 
do cat alignment_OLG"$ORF_num".fasta | paste - - | grep "$chr"  | shuf -n 1 ; done |
shuf -n $seq_set_size | 
awk -F "\t" '{print $1 "\n" $2}' > OLG"$ORF_num"_sample_set"$i".fasta ; 

# generate tree, and reroot on "outgroup" 

cat OLG"$ORF_num"_sample_set"$i".fasta outgroup_seq.fasta | gapless | 
awk -F "\t" '{if ($1~">") print; else print substr($1,1,(length($1)-3))}' > OLG"$ORF_num"_tree_set"$i".fasta ; 

ts_aln_aa_lgth=$( cat OLG"$ORF_num"_tree_set"$i".fasta | linear | awk '{if (NR==2) print length/3}' );


#################
### Find OLG positions in tree_set"$i" 
# note - only works if outgroup chromosome does not have too many gaps 
# (i.e. contains OLG region without gaps at start and stop)

# work out outgroup OLG co-ordinates without gaps, from original alignment seq
OLG_p1_outgroup_gapless=$(cat ../blastp/mGene"$ORF_num"_all-homologs_pal2nal.fasta  | 
grep -A1 $outgroup_chromosome | head -2 | 
awk -F "\t" '{if ($1!~">") print substr($1,1,("'${aln[0]}'"))}' | sed -e "s|-||g" | awk '{print length}' ) ;
OLG_p2_outgroup_gapless=$(cat ../blastp/mGene"$ORF_num"_all-homologs_pal2nal.fasta  | 
grep -A1 $outgroup_chromosome | head -2 | 
awk -F "\t" '{if ($1!~">") print substr($1,1,("'${aln[1]}'"))}' | sed -e "s|-||g" | awk '{print length}' ) ;

# corresponding positions in tree_set"$i"

adjusted_positions=$(cat OLG"$ORF_num"_tree_set"$i".fasta | grep -A1 $outgroup_chromosome | head -2 | 
grep -v ">" | fold -w 1 | awk '{print $1 "\t" NR}' | awk '{if ($1!="-") print}' | 
awk '{if (NR==("'$OLG_p1_outgroup_gapless'"+0) || NR==("'$OLG_p2_outgroup_gapless'"+0)) print $2 }' | 
tr '\n' ' ' | awk '{print}')   ;
adj=($adjusted_positions) ;

test -z ${adj[0]} && adj[0]=$OLG_p1_outgroup_gapless ;
test -z ${adj[1]} && adj[1]=$OLG_p2_outgroup_gapless ;

#################

# pick random non-outgroup set as "ancestor" to evolve from
# store gap loci co-ordinates to replace with gaps in simulated output 

#cat tree_set"$i".fasta | paste - - | grep -v "$outgroup_chromosome" | shuf -n 1 | 
#awk -F "\t" '{print $1 "\n" $2}' | tee anc_seq"$i".fasta | grep -v ">" | 
#awk -F "" '{for (i=1;i<=NF;i++) if($i=="-") print i}' > anc_seq_gapped_sites_set_"$i".txt ;

#ancestor_seq_chromosome=$(cat anc_seq"$i".fasta | awk -F ":" '{if (NR==1) print $1}' | sed -e "s|>||g"  ) ;

ancestor_seq_chromosome=$outgroup_chromosome ;
cat OLG"$ORF_num"_tree_set"$i".fasta | grep -A1 "$ancestor_seq_chromosome" > OLG"$ORF_num"_anc_seq"$i".fasta ;

cat OLG"$ORF_num"_anc_seq"$i".fasta | seqkit seq -p -r | seqkit translate | linear | grep -v ">" | 
awk -F "" '{for (i=1;i<=NF;i++) if($i=="-" || $i=="*") 
print ("'$ts_aln_aa_lgth'"-i+1) }' | sort -gk1,1 > OLG"$ORF_num"_anc_seq_gap-or-stop_sites_set_"$i".txt ;


#iqtree -redo -s tree_set.fasta -m HKY -asr ; 
echo "calculating tree $i" ;
iqtree -redo -s OLG"$ORF_num"_tree_set"$i".fasta -st CODON -m KOSI07 >/dev/null 2>&1 ;
echo "tree $i calculated" ;
cat OLG"$ORF_num"_tree_set"$i".fasta.treefile ;

outgroup_header=$(nw_labels -I OLG"$ORF_num"_tree_set"$i".fasta.treefile | grep $outgroup_chromosome) ; 

#reroot on outgroup 
#nw_reroot tree_set"$i".fasta.treefile $outgroup_header > tree_set_rw"$i".treefile ;
cp OLG"$ORF_num"_tree_set"$i".fasta.treefile OLG"$ORF_num"_tree_set_rw"$i".treefile ;

#can prune outgroup, but too conservative, and not needed when 'local' seq used:
#nw_prune -  $outgroup_header


# sequence to use as ancestor: 
# replace any gaps with GGG 
#[least likely to result in a stop in antisense; these sites are later excluded, in any case] - 

anc_seq_ntds=$(cat OLG"$ORF_num"_anc_seq"$i".fasta | grep -v ">" | sed -e "s|-|G|g") ;
echo "anc seq chromosome:" $ancestor_seq_chromosome ;
echo "anc seq OLG region:" $(cat OLG"$ORF_num"_anc_seq"$i".fasta | awk '{if ($1~">") print; 
else print substr($1,"'${adj[0]}'",("'${adj[1]}'"-"'${adj[0]}'"+1))}' |
seqkit seq -p -r | seqkit translate | linear | awk '{if ($1!~">") print }' ) ;

# count original stops; exclude outgroup sequence ; 

echo "tree set example seq" ;
cat OLG"$ORF_num"_tree_set"$i".fasta | linear | paste - - | grep -v "$outgroup_chromosome" | 
awk -F "\t" '{print $1 "\n" $2}' | seqkit seq -p -r | seqkit translate | linear | grep -v ">" | head -1 | 
awk -F "\t" '{pos1=("'$ts_aln_aa_lgth'"-("'${adj[1]}'"/3)+1) ; pos2=("'$ts_aln_aa_lgth'"-("'${adj[0]}'"/3)+1) ;
if ($1~">") print ; else print substr($1, pos1, ((pos2-pos1)+1))}' ;

echo "tree set average stops per seq" ;
cat OLG"$ORF_num"_tree_set"$i".fasta | linear | paste - - | grep -v "$outgroup_chromosome" | 
awk -F "\t" '{print $1 "\n" $2}' | seqkit seq -p -r | seqkit translate | linear | grep -v ">" | 
awk -F "" '{for (i=1;i<NF; i++) {{total=0} {if($i=="*") {total++}}; {print i "\t" total}}}' | 
sort | uniq -c | awk '{if ($3=="1") print ("'$ts_aln_aa_lgth'"-$2+1) "\t" $1/("'$seq_set_size'"+1)}' | 
sort -gk1,1 | tee OLG"$ORF_num"_tree_set"$i"_stops.txt | 
awk -F"\t" 'NR==FNR{a[$1]=$1;next}{ if (!a[$1])print ;}' OLG"$ORF_num"_anc_seq_gap-or-stop_sites_set_"$i".txt - |
tee OLG"$ORF_num"_tree_set"$i"_stops_filtered.txt |
awk '{if ( $1>=(("'${adj[0]}'"+2)/3) && $1<=(("'${adj[1]}'")/3) ) print}' | awk '{sum+=$2} END {print sum}' ;


seq $simulations | while read -r j ; 
do 
# pyvolve sequences along this tree using the same model, with "$iterations" number of simulations

if [ $codonfreqs == "customfreqs" ] ; then 
python3 ../scripts_long_OLGs/pyvolve_simulation_custom-freqs.py \
OLG"$ORF_num"_tree_set_rw"$i".treefile $anc_seq_ntds \
$model_type $omega_value rest_of_alignment_OLG"$ORF_num".fasta >/dev/null 2>&1 ;

elif [ $codonfreqs == "defaultfreqs" ] ; then 
python3 ../scripts_long_OLGs/pyvolve_simulation_default-freqs.py \
OLG"$ORF_num"_tree_set_rw"$i".treefile $anc_seq_ntds \
$model_type $omega_value >/dev/null 2>&1 ;
fi ; 

#python3 ../scripts_long_OLGs/pyvolve_simulation.py tree_set_rw"$i".treefile \
#$anc_seq_ntds $model_type $omega_value GCF_000006765.1_ASM676v1_cds_from_genomic.fna >/dev/null 2>&1 ;
if [ -s simulated_alignment.fasta ]; then 
mv simulated_alignment.fasta OLG"$ORF_num"_simulated_alignment_"$i"_"$j".fasta ;
else echo "tree iteration $i simulation $j failed" ; fi ; 

# count evolved stops
# remove gap sites 

echo "iteration $i simulation $j average stops per seq" ;  
cat OLG"$ORF_num"_simulated_alignment_"$i"_"$j".fasta | linear | paste - - | 
grep -v "$outgroup_chromosome" | awk -F "\t" '{print $1 "\n" $2}' | 
seqkit seq -p -r | seqkit translate | linear | grep -v ">" | 
awk -F "" '{for (i=1;i<NF; i++) {{total=0} {if($i=="*") {total++}}; {print i "\t" total}}}' | 
sort | uniq -c | awk '{if ($3=="1") print ("'$ts_aln_aa_lgth'"-$2+1) "\t" $1/("'$seq_set_size'")}' | 
sort -gk1,1 | tee OLG"$ORF_num"_simulated_alignment_"$i"_"$j"_stops.txt | 
awk -F"\t" 'NR==FNR{a[$1]=$1;next}{ if (!a[$1])print ;}' OLG"$ORF_num"_anc_seq_gap-or-stop_sites_set_"$i".txt - |
tee OLG"$ORF_num"_simulated_alignment_"$i"_"$j"_stops_filtered.txt | 
awk '{if ( $1>=(("'${adj[0]}'"+2)/3) && $1<=(("'${adj[1]}'")/3) ) print}' | awk '{sum+=$2} END {print sum}' ;

echo "number of stop-free sequences"
#original stop excluded by adding +3 to ${adj[0]}
cat OLG"$ORF_num"_simulated_alignment_"$i"_"$j".fasta | linear | paste - - | grep -v "$outgroup_chromosome" | 
awk -F "\t" '{print $1 "\n" substr($2,("'${adj[0]}'"+3),("'${adj[1]}'"-"'${adj[0]}'"+1))}' |
seqkit seq -p -r | seqkit translate | linear | paste - - | awk -F "\t" '{if ($2!~"*") print}' | wc -l ; 


done ; done > OLG"$ORF_num"_summary_output_0.5.txt ; 

# from this summary output can calculate percentages of seqs with stops, 
# taking each simulation as a new data point, to plot on a histogram 

cat OLG"$ORF_num"_summary_output_0.5.txt | grep -A1 "stop-free" | grep -v "-" | sort | uniq -c | sort -gk2,2 | 
awk '{print $2 "\t" $1}' |  
#add in extra zero data row to avoid error of not plotting "zero" values
awk '{print "-1" "\t" "0" "\n" $0}' | sort -gk1,1 | uniq > sim.txt ; 


for i in {1..100} ; do 
cat alignment_OLG"$ORF_num".fasta | awk -F "\t" '{if ($1~">") print;
else print substr($1,("'${aln[0]}'"+3),("'${aln[1]}'"-"'${aln[0]}'"-2))}' | 
seqkit seq -p -r | seqkit translate | linear | grep -v ">" | shuf -n 20 | grep "*" | wc -l | 
awk '{print 20-$1}' ; done | sort | uniq -c | sort -gk2,2 | awk '{print $2 "\t" $1}' > real_sampled.txt ; 

mv real_sampled.txt real_sampled_OLG"$ORF_num"_"$codonfreqs"_omega"$omega_value".txt ;
mv sim.txt sim_OLG"$ORF_num"_"$codonfreqs"_omega"$omega_value".txt ;


#gnuplot



if [ $ORF_num == 1 ] ; then 
OLGcolor="#005AB5" ; 
elif [ $ORF_num == 2 ] ; then 
OLGcolor="#DC3220" ; 
fi


gnuplot -e \
'dpi=600;
width=93;
height=40;

in2mm=25.4 ; 
pt2mm=0.3528; 
mm2px=dpi/in2mm ;
ptscale=pt2mm*mm2px ; 
round(x) = x - floor(x) < 0.5 ? floor(x) : ceil(x) ;
wpx = round(width * mm2px) ; hpx = round(height * mm2px) ;

set terminal pngcairo size wpx,hpx fontscale ptscale linewidth ptscale pointscale ptscale ;

set output "mGene'$ORF_num'_simulation-results-'$codonfreqs'-O'$omega_value'_og_'$outgroup_chromosome'.png" ;
set lmargin 4 ; 
set bmargin 2 ;
set border 1 lw 0.5 ; 
set xlabel "intact full-length ORFs" offset 0,1.2 font "arial, 7" ;
set ylabel "% of samples" offset 3 font "arial, 7" ; 
set xrange [-0.5:20.5] ; 
set ytics font ", 5" ; 
set xtics font ", 5" ; 
set xtics scale 0.1 nomirror offset 0,0.5 ; 
set ytics scale 0 offset 0.5,0;
set tics front ;
plot "real_sampled_OLG'$ORF_num'_'$codonfreqs'_omega'$omega_value'.txt" using 1:2 with boxes fc rgb "'$OLGcolor'" fs solid noborder lc rgb "black" noti, "sim_OLG'$ORF_num'_'$codonfreqs'_omega'$omega_value'.txt" using 1:2 with boxes lw 0.5 fc rgb "white" fs transparent solid 0.5 border lc rgb "black" noti' 















### Prev

set terminal pdf size 6cm,2.5cm 
set output 'simulation_boxplots.pdf'
set lmargin 3
set rmargin 2
set boxwidth 1
set style fill solid
set ytics font ", 7"
set xtics font ", 7"
set xtics scale 0.3
set ytics scale 0.3
set tics back
set xrange [0:20.5]
set yrange [0:100] 
plot "real_sampled.txt" using 1:2 with boxes lc rgb "#005AB5" fs transparent solid 0.5 noborder noti,\
"sim.txt" using 1:2 with boxes lc rgb "#FFA500" fs transparent solid 0.5 noborder noti,
set output 

pdftoppm -png simulation_boxplots.pdf simulation_boxplots ;







# use tree1 simulations as an example: 

length=$(cat tree_set1.fasta | head -2 | grep -v ">" | awk '{print length/3}' ) ;
seq $length | while read -r position ; do 
cat simulated_alignment_1_*stops_filtered.txt | awk '{if ($1=="'$position'") print}' | 
awk '{sum+=$2} END{print sum/"'$seq_set_size'"}' ; done |  awk '{print NR "\t" $1}' > sim1_stop.txt ; 

set terminal pdf size 14cm,4cm 
set output 'sim1_stops_filtered.pdf'
set style fill solid 1.0
set xtics font ",8" 
set ytics font ",8" 
set xtics scale 0
set ytics scale 0
set yrange [0:0.1] 
set boxwidth 0.5
set style fill solid
set arrow from 265, graph 0 to 265, graph 1 nohead lw 1 lc rgb "#8b0000"
set arrow from 582, graph 0 to 585, graph 1 nohead lw 1 lc rgb "#8b0000"
plot "sim1_stop.txt" using 1:2 with boxes lc rgb "#00008b" notitle 
set output 

pdftoppm -png sim1_stops_filtered.pdf sim1_stops_filtered ;


full_aln_aa_length=$(cat alignment.fasta | seqkit translate | linear | awk '{if (NR==2) print length}' ) ; 
aln_set_size=$(cat alignment.fasta | grep ">" | wc -l );
cat alignment.fasta | linear | paste - - | awk -F "\t" '{print $1 "\n" $2}' | 
seqkit seq -p -r | seqkit translate | linear | grep -v ">" | 
awk -F "" '{for (i=1;i<NF; i++) {{total=0} {if($i=="*") {total++}}; {print i "\t" total}}}' | 
sort | uniq -c | awk '{if ($3=="1") print ("'$full_aln_aa_length'"-$2+1) "\t" $1/("'$aln_set_size'"+1)}' | 
sort -gk1,1 > aln_stops.txt 





################

# remove outgroup seq, create new alignment without gap only columns

cat tree_set.fasta | paste - - | grep -v "CP034783.1-1898421-1900293--2" | 
awk -F "\t" '{print $1"\n" $2}' > tree_set2.fasta ;

# create tree with IQTREE

iqtree -s tree_set2.fasta -m HKY ; 

# pyvolve sequences along this tree using the same model, with "$iterations" number of simulations

sudo python3 /mnt/c/Users/Zachary/Dropbox/OLGs/Scripts/MK/Pseudomonas_long_ORFs/pyvolve_working.py ;






### EXTRA WORKING:



cat alignment.fasta | linear | awk -F "\t" '{if($1~">")print; else print substr($1,949,1407)}' | 
sed -e "s|:||g" | sed -e "s|[()]||g" | sed -e "s|+||g" > tmp3.fasta




cat alignment.fasta | grep ">" | grep -v "GUT" | grep -v "a:Ga" | awk -F ":" '{print $1}' | sed -e "s|>||g" | 
while read -r query ; do
esearch -db nuccore -query "$query"  < /dev/null | elink -target taxonomy | efetch -format xml |
xtract -pattern Taxon -element ScientificName 


