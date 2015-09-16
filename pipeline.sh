#!/bin/bash

# directory for treefam files
TREEFAM=data/treefam

# data dump
TFDUMP=$TREEFAM/treefam_family_data.tar.gz

# extract and pre-process trees and alignments from treefam
perl script/treefammer.pl -i $TFDUMP -t $TREEFAM -v

# list of input trees from treefam
INTREES=`ls $TREEFAM/*.nhx.emf`

# directory to read/write fossils
FOSSILS=data/fc/

# directory to write ratograms
RATOS=data/ratograms/

# where the r8s executable is
R8S=deps/dist/r8s

# template to generate r8s commands
TMPL=data/r8s.tmpl

# iterate over input trees
#for INTREE in $INTREES; do
#	perl script/ratogrammer.pl -f $FOSSILS -r $RATOS -e $R8S -t $TMPL -v -i $INTREE
#done
