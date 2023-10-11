# morpho_characters
Long overdue scripts for analyzing morphological characters

Saturation: cite https://doi.org/10.1111/j.1469-185X.2012.00240.x
Understanding phylogenetic incongruence: lessons from phyllostomid bats.
Liliana M. Dávalos, Andrea L. Cirranello, Jonathan H. Geisler, Nancy B. Simmons 
Biological Reviews 87(4), 991-1041

Simulation, Gower distances, and statistical scaffold: cite https://doi.org/10.1093/sysbio/syu022
Integrating Incomplete Fossils by Isolating Conflicting Signal in Saturated and Non-Independent Morphological Characters
Liliana M. Dávalos, Paúl M. Velazco, Omar M. Warsi, Peter D. Smits, Nancy B. Simmons 
Systematic Biology, Volume 63, Issue 4, July 2014, 582–600

Gower distances:

Running dissimilarity_morphology.R requires running simulation_morphology.R first. 

Statistical scaffold:

Step 1: Running rax_scaffold.pl estimates backbone constrained morphological trees and requires a morphology character file, outgroup specification, and a file with Newick sample of molecular trees. 

Step 2: Running rax_mulopt.pl estimates branch lengths from both morphological and molecular data on backbone constrained morphological trees (from step 1) and requires a combined morphology and molecular character file, outgroup specification, a partition file, and a file with Newick set of best trees obtained in the step 1.