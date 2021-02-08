#!/bin/python3

import pyvolve ; import sys

tree_variable=sys.argv[1]
anc_seq_variable=sys.argv[2]
model_type=sys.argv[3]
omega_value=float(sys.argv[4])



# Simulation:

my_tree = pyvolve.read_tree(file = tree_variable)

my_model = pyvolve.Model(model_type, {"omega": omega_value })

my_partition = pyvolve.Partition(models = my_model, root_sequence = anc_seq_variable)

my_evolver = pyvolve.Evolver(tree = my_tree, partitions= my_partition)
my_evolver()   

#pyvolve.print_tree(tree_variable)