#!/bin/bash

# OLGenie 


#Variables:

if [ $ORF_num == 1 ] ; then 
OLG_pos1=793 ; OLG_pos2=1749 ; 
elif [ $ORF_num == 2 ] ; then 
OLG_pos1=265 ; OLG_pos2=1800 ; 
fi


adjusted_positions=$(cat ../blastp/mGene"$ORF_num"_all-homologs_pal2nal2.fasta | grep -A1 $reference | head -2 | 
grep -v ">" | fold -w 1 | awk '{print $1 "\t" NR}' | awk '{if ($1!="-") print}' | 
awk '{if (NR==("'$OLG_pos1'"+0) || NR==("'$OLG_pos2'"+0)) print $2 }' | 
tr '\n' ' ' | awk '{print}')   ;
adj=($adjusted_positions) ;



##############


# p values matrix (ASL p-values): 

cat ../blastp/mGene"$ORF_num"_all-homologs_pal2nal2.fasta | linear | awk -F "\t" '{if ($1~">") print; 
else print substr($1,"'${adj[0]}'",("'${adj[1]}'"-"'${adj[0]}'"+1))}' > OLG"$ORF_num"_reg_ntd.fasta ;


cat ../blastp/mGene"$ORF_num"_all-homologs_pal2nal.fasta.contree_trimmed_list_X_20 | 
awk '{print $1 "\t" NR}' | while read header rownr; 
do cat OLG"$ORF_num"_reg_ntd.fasta | grep -A1 $header > tmp.fasta ;  
cat ../blastp/mGene"$ORF_num"_all-homologs_pal2nal.fasta.contree_trimmed_list_X_20 | 
awk '{print $1 "\t" NR}' | while read header2 rownr2; 
do cat OLG"$ORF_num"_reg_ntd.fasta | grep -A1 $header2 > tmp2.fasta ;  
cat tmp.fasta tmp2.fasta > tmp3.fasta ; 

perl ~/OLGenie/OLGenie.pl --fasta_file=tmp3.fasta --frame=sas13 \
--output_file=OLGenie-out.tsv --verbose >/dev/null 2>&1 ;
#awk '{if ($1~"dNN/dNS") print}' | sed -e "s|dNN/dNS=||g"   ; 

Rscript ~/OLGenie/OLGenie_bootstrap.R OLGenie-out.tsv 2 1000 4 | tee bootstrap.txt |
awk -F "\t" '{if (NR==3 && $16=="NaN") print "1" ;
else if (NR==3 && $16>=1) print 1/$29 ;
else if (NR==3 && $16<1) print -1/$29 }' ;

done | tr '\n' '\t' | awk '{print}' ; 
done > OLGenie_OLG"$ORF_num"_matrix2.txt ; 
cat OLGenie_OLG"$ORF_num"_matrix2.txt | sed -e "s|NaN|1|g" > OLGenie_OLG"$ORF_num"_matrix_p-values.txt ;


gnuplot -e 'set terminal pdf size 7,6 ;
set output "ASL-matrix_OLG'$ORF_num'.pdf"; 
set lmargin 2 ;
set rmargin 10 ;
set size square ;
set xrange [-0.5:19.5] ;
set yrange [19.5:-0.5] ;
set palette  defined (-1000 "dark-green", 1 "white", 1000 "dark-yellow") ;
set cbrange [-100:100] ; 
set xtics offset -0.5,0 ;
set format x "" ; 
set format y "" ; 
set format cb "" ;
set colorbox vertical user origin .9,.3 size .04,.4 ;
plot "OLGenie_OLG'$ORF_num'_matrix_p-values.txt" matrix w image noti ' ;

pdftoppm -png ASL-matrix_OLG"$ORF_num".pdf ASL-matrix_OLG"$ORF_num" ;




####### 

# sliding window:

cat ../Simulation/alignment_OLG"$ORF_num".fasta | 
sed -e "s|:|_|g" | sed -e "s|\.|_|g" | sed '/^>/s/ .*//' | sed '/^>/s/-//' |
sed -e "s|[(+)]||g" > tmp_pal2nal_OLG"$ORF_num".fasta ;

perl ~/OLGenie/OLGenie.pl --fasta_file=tmp_pal2nal_OLG"$ORF_num".fasta --frame=sas13 \
--output_file=OLGenie-out_"$ORF_num"_sas13.tsv --verbose >/dev/null 2>&1 ;
Rscript ~/OLGenie/OLGenie_sliding_windows.R OLGenie-out_"$ORF_num"_sas13.tsv \
 NN NS 50 10 1000 2 NONE 4 > OLGenie_sliding_windows_OLG"$ORF_num"_sas13.out

#perl ~/OLGenie/OLGenie.pl --fasta_file=tmp_pal2nal_OLG"$ORF_num".fasta --frame=sas12 \
#--output_file=OLGenie-out_"$ORF_num"_sas12.tsv --verbose >/dev/null 2>&1 ;
#Rscript ~/OLGenie/OLGenie_sliding_windows.R OLGenie-out_"$ORF_num"_sas12.tsv \
# NN NS 50 1 1000 2 NONE 4 > OLGenie_sliding_windows_OLG"$ORF_num"_sas12.out

#perl ~/OLGenie/OLGenie.pl --fasta_file=tmp_pal2nal_OLG"$ORF_num".fasta --frame=sas11 \
#--output_file=OLGenie-out_"$ORF_num"_sas11.tsv --verbose >/dev/null 2>&1 ;
#Rscript ~/OLGenie/OLGenie_sliding_windows.R OLGenie-out_"$ORF_num"_sas11.tsv \
# NN NS 50 1 1000 2 NONE 4 > OLGenie_sliding_windows_OLG"$ORF_num"_sas11.out


cat OLGenie-out_"$ORF_num"_sas13_WINDOWS_dNNdNS.tsv | tail -n +2 | awk '{print NR "\t" $(NF-13)}' |
awk -F "\t" '{if ($2~"[0-9]" && $2>100 || $2=="Inf") print $1 "\t" "100"; 
else if ($2~"[0-9]" && $2<0.01) print $1 "\t" "0.01";
else print}' > tmp_sas13.dat ;
#cat OLGenie-out_"$ORF_num"_sas12_WINDOWS_dNNdNS.tsv | awk '{print NR "\t" $(NF-13)}' > tmp_sas12.dat ;
#cat OLGenie-out_"$ORF_num"_sas11_WINDOWS_dNNdNS.tsv | awk '{print NR "\t" $(NF-13)}' > tmp_sas11.dat ;


####
# control sliding windows

cat ../Simulation/rest_of_alignment_OLG"$ORF_num".fasta | 
sed -e "s|:|_|g" | sed -e "s|\.|_|g" | sed '/^>/s/ .*//' | sed '/^>/s/-//' |
sed -e "s|[(+)]||g" > tmp_pal2nal_OLG"$ORF_num"_control.fasta ;

perl ~/OLGenie/OLGenie.pl --fasta_file=tmp_pal2nal_OLG"$ORF_num"_control.fasta  --frame=sas13 \
--output_file=OLGenie-out_"$ORF_num"_sas13_control.tsv --verbose >/dev/null 2>&1 ;
Rscript ~/OLGenie/OLGenie_sliding_windows.R OLGenie-out_"$ORF_num"_sas13_control.tsv \
 NN NS 50 10 1000 2 NONE 4 > OLGenie_sliding_windows_OLG"$ORF_num"_control.out

cat OLGenie-out_"$ORF_num"_sas13_control_WINDOWS_dNNdNS.tsv | awk '{print NR "\t" $(NF-13)}' > tmp_sas13_control.dat ;

##########

#sliding window figures (repeat this, and recalculate input files, for each OLG (ORF_num) )

alg_length=$(cat ../Simulation/alignment_OLG"$ORF_num".fasta  | linear | awk '{if (NR==2) print length/3}') ; 


if [ $ORF_num == 1 ] ; then 
OLGcolor="#005AB5" ; 
mGene=tle3 ;
elif [ $ORF_num == 2 ] ; then 
OLGcolor="#DC3220" ; 
mGene=PA1383 ; 
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

set output "OLGenie_sliding_windows_'$ORF_num'.png"; 
set lmargin 10; 
set lmargin 5; 
unset colorbox ;
set logscale y ; 
set xlabel "'$mGene' codon number" offset 0,1.2 font "arial, 7" ;
set ylabel "dNN/dNS" offset 5 font "arial, 7" ; 
set xtics scale 0.5 ; 
set ytics scale 0.5 ; 
set ytics font "arial, 5" ;
set xtics font "arial, 5" ;
set xtics offset 0,0.5 ;
set xrange [0:'$alg_length'] ; 
set yrange [0.01:100] ; 
set arrow from 0,1 to '$alg_length',1 lw 1.5 dt 3 lc rgb "black" nohead ;
set style rectangle fillcolor rgb "'$OLGcolor'" ;
set object 1 rect from '${adj[0]}'/3,0 to ('${adj[1]}'/3),100 fs transparent solid 0.7 back ; 
set style rectangle fillcolor rgb "#888888" ;
set object 2 rect from 0,0 to '${adj[0]}'/3,100 fs transparent solid 0.7 back ; 
set object 3 rect from ('${adj[1]}'/3),0 to '$alg_length',100 fs transparent solid 0.7 back ; 
plot "tmp_sas13_control.dat" u 1:2 with lines lt 1 lw 1.5 lc rgb "white" notitle, "tmp_sas13.dat" u 1:2 with lines lt 1 lw 1.5 lc rgb "#FFFF33" notitle ' ;


######## 



# single dNN/dNS values, reported in manuscript text:

cat ../Simulation/alignment_OLG"$ORF_num".fasta | linear | awk -F "\t" '{if ($1~">") print; 
else print substr($1,"'${adj[0]}'",("'${adj[1]}'"-"'${adj[0]}'"+1))}' | sed -e "s|[+():]||g" > tmp3.fasta ;
perl ~/OLGenie/OLGenie.pl --fasta_file=tmp3.fasta --frame=sas13 --output_file=OLGenie-out.tsv --verbose ;
Rscript ~/OLGenie/OLGenie_bootstrap.R OLGenie-out.tsv 2 1000 4 | tee bootstrap.txt |
awk -F "\t" '{print $16 "\t" $23 "\t" $29 }' > alignment_OLG"$ORF_num"_bootstrap.txt ; 
# control genomes
cat ../Simulation/rest_of_alignment_OLG"$ORF_num".fasta | linear | awk -F "\t" '{if ($1~">") print; 
else print substr($1,"'${adj[0]}'",("'${adj[1]}'"-"'${adj[0]}'"+1))}' | sed -e "s|[+():]||g" > tmp3.fasta ;
perl ~/OLGenie/OLGenie.pl --fasta_file=tmp3.fasta --frame=sas13 --output_file=OLGenie-out.tsv --verbose ;
Rscript ~/OLGenie/OLGenie_bootstrap.R OLGenie-out.tsv 2 1000 4 | tee bootstrap.txt |
awk -F "\t" '{print $16 "\t" $23 "\t" $29}' > rest_of_alignment_OLG"$ORF_num"_bootstrap.txt ; 

