!/bin/bash

#need "codon_table.rds" file in working directory

#FUNCTIONS:
linear () { awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }'; } ;


time for gene in mGene1_ntd.fasta mGene2_ntd.fasta ; do
geneseq=$(cat $gene | linear | awk '{if (NR==2) print }') ;
Rscript ../scripts_long_OLGs/Frameshift_20000_revcom0.r $geneseq | 
awk '{if (NF==7 || $1~"#") print }' > ${gene%.*}-frameshift.txt ;
done ;


for ORF_num in 1 2 ; do 
cat  mGene"$ORF_num"_ntd-frameshift.txt | head -1 | tr '\t' '\n' | grep -v "#" > mGene"$ORF_num"_perm.txt ; 
cat  mGene"$ORF_num"_ntd-frameshift.txt | head -2 | tail -1 | tr '\t' '\n' | grep -v "#" > mGene"$ORF_num"_syn.txt ; 
done ; 



# GNUPLOT: 



if [ $ORF_num == 1 ] ; then 
OLGcolor="#005AB5" ; 
overlap_length=1217 ; 
max_length=1800 ; 
elif [ $ORF_num == 2 ] ; then 
OLGcolor="#DC3220" ; 
overlap_length=1541 ; 
max_length=1800 ; 
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

set output "mGene'$ORF_num'_ORFlengths_combined.png" ;
set lmargin 4 ; 
set bmargin 2 ;
set border 1 lw 0.5 ; 
set xlabel "simulated ORF lengths (codons)" offset 0,1.2 font "arial, 7" ;
set ylabel "ORFs freq" offset 5 font "arial, 7" ; 
binwidth=20 ;
bin(x,width)=width*floor(x/width) ;
stats "mGene'$ORF_num'_perm.txt" ;
n_perm=STATS_records ;
stats "mGene'$ORF_num'_syn.txt" ;
n_syn=STATS_records ;
set xrange [0:'$max_length'] ; 
set ytics font ", 5" ; 
set xtics font ", 5" ; 
set xtics scale 0.1 nomirror offset 0,0.5 ; 
set ytics scale 0 offset 0.5,0;
set tics front ;
set arrow from '$overlap_length', graph 0 to '$overlap_length', graph 1 nohead lw 1.5 dt 3 lc rgb "'$OLGcolor'" front ;
plot "mGene'$ORF_num'_syn.txt" using (bin($1,binwidth)):(1.0/n_syn) smooth freq with boxes lc rgb "#FFA500" fs transparent solid 0.7 noborder notitle, "mGene'$ORF_num'_perm.txt" using (bin($1,binwidth)):(1.0/n_perm) smooth freq with boxes lc rgb "#228833" fs transparent solid 0.5 noborder notitle' ; 

