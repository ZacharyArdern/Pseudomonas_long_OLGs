from ete3 import PhyloTree, SeqGroup, SequenceFace, TreeStyle, AttrFace, NodeStyle, Tree
import sys


alignment_input=sys.argv[1]
tree_input=sys.argv[2]


alg=( alignment_input )



t = PhyloTree( tree_input , format=1, quoted_node_names=True )
seqs = SeqGroup(alg, format="fasta")


nodestyle1 = NodeStyle()
nodestyle1["size"] = 0
nodestyle1["vt_line_width"] = 2
nodestyle1["hz_line_width"] = 2

for node in t.traverse():
    node.set_style(nodestyle1)


for leaf in t.iter_leaves():
    item=seqs.get_seq(leaf.name)
    name_face = AttrFace(item, fsize=24)
    Bars = SequenceFace(item, seqtype='aa', fsize=24, bg_colors={'G': 'Khaki', 'A': 'Khaki', 'S': 'Khaki', 'T': 'Khaki', 'C': 'LightGreen', 'V': 'LightGreen', 'I': 'LightGreen', 'L': 'LightGreen', 'P': 'LightGreen', 'F': 'LightGreen', 'Y': 'LightGreen', 'M': 'YellowGreen', 'W': 'LightGreen', 'N': 'Thistle', 'Q': 'Thistle', 'H': 'Thistle', 'D': 'DarkSalmon', 'E': 'DarkSalmon', 'K': 'SkyBlue', 'R': 'SkyBlue', 'X':'Black', '-':'White' }, fg_colors=None, codon=None, col_w=1.5, alt_col_w=3, special_col=None, interactive=False)
    leaf.add_face(Bars, 2, "aligned")    
    
t.render("tree_and_alignment.png", h=100, units="mm")
t.render("tree_and_alignment.svg", h=100, units="mm")





t2 = PhyloTree( tree_input , format=1, quoted_node_names=True )
for node in t2.traverse():
    node.set_style(nodestyle1)
t2.convert_to_ultrametric(tree_length=None, strategy='balanced')
cladogram = TreeStyle()
cladogram.scale = 20
t2.render("cladogram.png", h=75, units="mm", tree_style=cladogram)

for leaf in t2.iter_leaves():
    leaf.name = " "
t2.render("cladogram_nameless.png", h=75, units="mm", tree_style=cladogram)




#colour options here: http://etetoolkit.org/docs/latest/reference/reference_treeview.html?highlight=colors#ete3.SVG_COLORS

#colours adapted from: http://www.bioinformatics.nl/~berndb/aacolour.html
#Lesk
#The colour scheme in Lesk, Introduction to Bioinformatics, uses 5 groups (note Histidine):
#Small nonpolar	G, A, S, T	Orange
#Hydrophobic	C, V, I, L, P, F, Y, M, W	Green
#Polar	N, Q, H	Magenta
#Negatively charged	D, E	Red
#Positively charged	K, R	Blue

# with M changed to 'Blue '

# List here: http://etetoolkit.org/docs/latest/reference/reference_treeview.html#color-names

