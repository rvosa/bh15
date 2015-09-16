#!/bin/bash

# list of input trees from treefam
INTREES=`ls data/treefam/*.nhx.emf`

# directory to read/write fossils
FOSSILS=data/fc/

# directory to write ratograms
RATOS=data/ratograms/

# where the r8s executable is
R8S=deps/dist/r8s

# template to generate r8s commands
TMPL=data/r8s.tmpl

# iterate over input trees
for INTREE in $INTREES; do
	perl script/ratogrammer.pl -f $FOSSILS -r $RATOS -e $R8S -t $TMPL -v -i $INTREE
done