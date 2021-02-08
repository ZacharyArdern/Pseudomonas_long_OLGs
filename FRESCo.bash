#!/bin/bash


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


###############

# create generic FRESCo input file

echo -e "inputRedirect = {};""\n"\
"inputRedirect["01"]="ALIGNMENT";""\n"\
"inputRedirect["02"]="TREE";""\n"\
"inputRedirect["03"]="50";""\n"\
"ASSUME_REVERSIBLE_MODELS = -1;""\n"\
"ExecuteAFile ("FRESCO.bf", inputRedirect);""\n" > FRESCOinput.bf ; 





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



cat ../Simulation/alignment_OLG"$ORF_num".fasta | 
sed -e "s|:|_|g" | sed -e "s|\.|_|g" | sed '/^>/s/ .*//' | sed '/^>/s/-//' |
sed -e "s|[(+)]||g" > tmp_pal2nal_"$ORF_num".fasta ;

file=tmp_pal2nal_"$ORF_num".fasta ; 

#OR:
FastTree -nt "$file" > ${file%.*}.tree ; 
cat FRESCOinput.bf | sed -e "s|ALIGNMENT|${file%.*}.fasta|g" | sed -e "s|TREE|${file%.*}.tree|g" > ${file%.*}_FRESCO.bf ;

#/home/zach/miniconda3/bin/hyphy ${file%.*}_FRESCO.bf > ${file%.*}_fresco4.out #path/to/hyphy path/to/fresco 
#or copy FRESCOinput.bf and FRESCO.bf files into the working directory

/usr/local/bin/hyphy ${file%.*}_FRESCO.bf > ${file%.*}_fresco3.out ;

cat  ${file%.*}_fresco3.out | tail -n +3 | awk '{print $1 "\t" $2}' > data_"$ORF_num".dat ;

### control

cat ../Simulation/rest_of_alignment_OLG"$ORF_num".fasta | 
sed -e "s|:|_|g" | sed -e "s|\.|_|g" | sed '/^>/s/ .*//' | sed '/^>/s/-//' |
sed -e "s|[(+)]||g" > tmp_pal2nal2_control"$ORF_num".fasta ;

file=tmp_pal2nal2_control"$ORF_num".fasta  ; 

#OR:
FastTree -nt "$file" > ${file%.*}.tree ; 
cat FRESCOinput.bf | sed -e "s|ALIGNMENT|${file%.*}.fasta|g" | sed -e "s|TREE|${file%.*}.tree|g" > ${file%.*}_FRESCO.bf ;

#/home/zach/miniconda3/bin/hyphy ${file%.*}_FRESCO.bf > ${file%.*}_fresco4.out #path/to/hyphy path/to/fresco 
#or copy FRESCOinput.bf and FRESCO.bf files into the working directory

/usr/local/bin/hyphy ${file%.*}_FRESCO.bf > ${file%.*}_fresco3.out ;

cat  ${file%.*}_fresco3.out | tail -n +3 | awk '{print $1 "\t" $2}' > control_"$ORF_num".dat ;



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

set output "FRESCo_sliding_windows_'$ORF_num'.png"; 
set lmargin 10; 
set lmargin 5; 
unset colorbox ;
set xlabel "'$mGene' codon number" offset 0,1.2 font "arial, 7" ; 
set ylabel "syn. rate" offset 4 font "arial, 7" ; 
set xtics scale 0.5 ; 
set ytics scale 0.5 ; 
set ytics font "arial, 5" ;
set xtics font "arial, 5" ;
set xtics offset 0,0.5 ;
set xrange [0:'$alg_length'] ; 
set yrange [0:3.5] ; 
set arrow from 0,1 to '$alg_length',1 lw 1.5 dt 3 lc rgb "black" nohead ;
set style rectangle fillcolor rgb "'$OLGcolor'" ;
set object 1 rect from '${adj[0]}'/3,0 to ('${adj[1]}'/3),100 fs transparent solid 0.7 back ; 
set style rectangle fillcolor rgb "#888888" ;
set object 2 rect from 0,0 to '${adj[0]}'/3,100 fs transparent solid 0.7 back ; 
set object 3 rect from ('${adj[1]}'/3),0 to '$alg_length',100 fs transparent solid 0.7 back ; 
plot "control_'$ORF_num'.dat" u 1:2 with lines lt 1 lw 1.5 lc rgb "white" notitle, "data_'$ORF_num'.dat" u 1:2 with lines lt 1 lw 1.5 lc rgb "#FFFF33" notitle ' ;











######### Prev


gnuplot -e 'set terminal pngcairo size 350,150;
set output "FRESCo_sliding_windows_'$ORF_num'.png"; 
set lmargin 10; 
set lmargin 5; 
unset colorbox ;
set xlabel "codon number" offset 0,graph 0.1 font "arial, 8" ; 
set ylabel "synonymous rate" offset 4 font "arial, 8" ; 
set xtics scale 0.5 ; 
set ytics scale 0.5 ; 
set ytics font "arial, 6" ;
set xtics font "arial, 6" ;
set xtics offset 0,graph 0.05 ;
set xrange [0:'$alg_length'] ; 
set yrange [0:3.5] ; 
set arrow from 0,1 to '$alg_length',1 lw 2 dt 2 lc rgb "black" nohead ;
set style rectangle fillcolor rgb '$OLGcolor' ;
set object 1 rect from '${adj[0]}'/3,0 to ('${adj[1]}'/3),3.5 fs transparent solid 0.7 back ; 
set style rectangle fillcolor rgb "#888888" ;
set object 2 rect from 0,0 to '${adj[0]}'/3,3.5 fs transparent solid 0.7 back ; 
set object 3 rect from ('${adj[1]}'/3),0 to '$alg_length',3.5 fs transparent solid 0.7 back ; 
plot "control_'$ORF_num'.dat" u 1:2 with lines lt 1 lw 2 lc rgb "white" notitle, "data_'$ORF_num'.dat" u 1:2 with lines lt 1 lw 2 lc rgb "#FFFF33" notitle ' ;





gnuplot -e 'set terminal pdf size 8.5,2.5 ;
set output "fresco_OLG'$ORF_num'-combined.pdf"; 
unset colorbox ;
set xrange [0:854] ; 
set arrow from 0,1 to 854,1 lw 3 dt 2 lc rgb "black" nohead ;
set arrow from '${adj[0]}'/3,1 to '${adj[1]}'/3,1 lw 15 lc rgb "#002E7F" nohead ;
plot "control_'$ORF_num'.dat" u 1:2 with lines lt 1 lw 2 lc rgb "#92DFF3" notitle, "data_'$ORF_num'.dat" u 1:2 with lines lt 1 lw 2 lc rgb "#005AB5" notitle' ;


pdftoppm -png fresco_OLG"$ORF_num"-combined.pdf fresco_OLG"$ORF_num"-combined ;






#from approx. codons 263-582

cat tmp_pal2nal3.fasta | seqkit seq -p -r | linear | 
awk '{if ($1~">") print ; else print substr($1,409,957)}'

gnuplot -e 'set terminal pdf size 7,4 ;
set output "fresco_OLG'$ORF_num'.pdf"; 
unset colorbox ;
plot "data.dat" u 1:2 with lines pt 7 ps 0.3 lt 1 lc rgb "blue" notitle' ;

gnuplot -e 'set terminal pdf size 7,4 ;
set output "fresco_OLG'$ORF_num'-control.pdf"; 
unset colorbox ;
plot "control.dat" u 1:2 with lines pt 7 ps 0.3 lt 1 lc rgb "blue" notitle' ;

