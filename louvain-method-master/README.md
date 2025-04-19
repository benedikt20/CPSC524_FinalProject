-----------------------------------------------------------------------------

Louvain method for community detection implemented in parallel 

The parallel version is compared to the sequential version from the article cited below

-----------------------------------------------------------------------------

Based on the article "Fast unfolding of community hierarchies in large networks"
Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre

-----------------------------------------------------------------------------
Author of the sequential algorithm:

Author   : E. Lefebvre, adapted by J.-L. Guillaume
Email    : jean-loup.guillaume@lip6.fr
Location : Paris, France
Time     : February 2008

From https://github.com/riyadparvez/louvain-method 
-----------------------------------------------------------------------------

Disclaimer:
If you find a bug, please send a bug report to jean-loup.guillaume@lip6.fr
including if necessary the input file and the parameters that caused the bug.
You can also send me any comment or suggestion about the program.

Note that the program is expecting a friendly use and therefore does not make
much verifications about the arguments.

-----------------------------------------------------------------------------


To execute:

1. Compile the communityPar.cpp script with the following command:
g++-12 -O2 -o community_par communityPar.cpp graph.cpp graph_binary.cpp main_community.cpp


2. Execute with the run_louvain.py script with the command:
python3 run_louvain.py

3. Summary output of the runtimes of the graphs is saved in a .txt file louvain_summary_par.txt

4. Output communities from the lovain method are saved in the folder output_communities

-----------------------------------------------------------------------------

Version history:
- First version April 18th 2025

