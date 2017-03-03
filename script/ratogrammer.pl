#!/usr/bin/perl
use strict;
use warnings;
use YAML;
use Template;
use Getopt::Long;
use List::Util 'sum';
use File::Temp 'tempfile';
use Statistics::Descriptive;
use Bio::Phylo::IO 'parse';
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::Util::CONSTANT ':objecttypes';
use Bio::Phylo::PhyLoTA::Service::CalibrationService;

# This script calibrates gene trees using r8s and calibration points from
# FossilCalibrations.org. The steps are roughly as follows:
# 1. check to see if there is a FASTA alignment that corresponds with the
#    tree. If yes, compute the number of sites in the alignment. Otherwise
#    the analysis cannot continue and the script consequently quits.
# 2. parse the input tree. In this step, outlying terminal branches are
#    that are more than 8 standard deviations from the mean are pruned.
#    This is a conservative measure intended to avoid pathologies in the 
#    rate smoothing step.
# 3. traverse the internal nodes and fetch fossil calibrations for them,
#    by accessing cached files and by querying the FossilCalibrations.org
#    web service.
# 4. apply the fossil calibrations to the internal nodes. The main challenge
#    here is that TreeFam is a it promiscuous in how it labels internal nodes:
#    sometimes, multiple, nested nodes are given the same label. Logically, 
#    this means that we should pick the deepest (nearest to the root) of 
#    these. The second challenge/assumption is that we apply the same fossil
#    to the appropriate speciation event in all paralogous copies.
# 5. run r8s. This uses the r8s.tmpl file to set up the commands block for
#    the analysis. Once the analysis is done, the ratogram is parsed out of
#    the log and written to an output directory.

# process command line arguments
my $fcdir;            # directory for cached calibrations
my $intree;           # input tree
my $exe = 'r8s';      # r8s executable
my $verbosity = WARN; # log level
my $template;         # template to generate r8s infile
my $ratodir;          # output directory for ratograms
my $dev = 8.0;        # stdev * $dev for outliers
GetOptions(
	'fcdir=s'    => \$fcdir,
	'intree=s'   => \$intree,
	'exe=s'      => \$exe,
	'template=s' => \$template,
	'verbose+'   => \$verbosity,
	'ratodir=s'  => \$ratodir,
	'dev=f'      => \$dev,
);

# instantiate helper objects
my $cs = Bio::Phylo::PhyLoTA::Service::CalibrationService->new;
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main',
);

# main program logic, functions follow
$log->info("going to analyze $intree");
my $nchar = get_nchar($intree);
my ($project,$tree) = read_data($intree);
my @files = get_fossil_files($tree);
my @fossils = get_fossils(@files);
my $result = run_r8s($project,$nchar,@fossils);
parse_result($result);

# check if aln available, exit if not, otherwise compute nchar
sub get_nchar {
	my $aln = shift;
	my $nchar;
	$aln =~ s/\.nhx\.emf/.cds.fasta/;
	if ( ! -e $aln ) {
		$log->warn("no alignment for $intree");
		exit(0);
	}
	else {
		open my $fh, '<', $aln or die $!;
		my ( $seq, $header );
		while(<$fh>) {
			chomp;
			if ( />/ ) {
				$header++ ? last : next;
			}
			$seq .= $_;
		}
		$nchar = length($seq);
	}
	$log->info("alignment has $nchar characters");
	return $nchar;
}

# read input tree file
sub read_data {
	my $intree = shift;

	# parse the tree
	$log->info("going to read $intree");
	my $project = parse(
		'-format'     => 'nhx',
		'-file'       => $intree,
		'-as_project' => 1,
	);

	# fetch tree from project
	my ($tree) = @{ $project->get_items(_TREE_) };
	
	# prune outlying terminal branch lengths
	my $test = 1;
	while( $test ) {
		my $stat = Statistics::Descriptive::Sparse->new;
		my @tips = @{ $tree->get_terminals };
		for my $tip ( @tips ) {
			$stat->add_data($tip->get_branch_length);
		}
		my $stdev = $stat->standard_deviation;
		my $mean  = $stat->mean;
		my @prune;
		for my $tip ( @tips ) {
			my $l = $tip->get_branch_length;
			if ( $l > ( $mean + $stdev * $dev ) or $l < ( $mean - $stdev * $dev ) ) {
				$log->warn("outlier: " . $tip->get_name);
				push @prune, $tip;
			}
		}
		$tree->prune_tips(\@prune);
		$test = scalar(@prune);
	}
	
	return $project, $tree;
}

# fetch calibration points for nodes we needed to do the
# relabeling of the tips first because the $cs uses this,
# and we do a pre-order traversal because TreeFam labels
# multiple, nested nodes with the same higher taxon label.
# we want the deepest of the nested nodes, which we visit
# first in this traversal.
sub get_fossil_files {
	my $tree = shift;
	my @files;
	$log->info("going to fetch fossils from web service");
	$tree->visit_depth_first(
		'-pre' => sub {
			my $node = shift;
			my $name = $node->get_name;
			my $file = "${fcdir}/${name}.yml";
			if ( $node->is_internal ) {

				# store fossil file name for later
				push @files, $file;
	
				# not yet processed this taxon, either during
				# this traversal or during a previous run			
				if ( ! -e $file ) {
		
					# fetch fossils
					$log->info("fetching fossils for $name");
					my @records = $cs->fetch_fossil_dates( $name );
					open my $fh, '>', $file or die $!;
					binmode( $fh, ':utf8' ); # wide characters otherwise
					print $fh Dump(\@records);
				}
				else {
					$log->debug("already fetched fossils for $name");
				}					
			}
		}
	);
	return @files;
}

# read fossils from file. XXX should test topology. doesn't, yet.
sub get_fossils {
	my @files = @_;
	
	# read distinct fossils from files. let's assume that
	# the same file can occur multiple times in the list,
	# and the same fossil can occur multiple times among files.
	my ( @fossils, %seen );
	$log->info("going to read fossil files");
	for my $file ( @files ) {
		next if $seen{$file}++;		
		open my $fh, '<', $file or die $!;
		my @set = @{ Load( do { local $/; <$fh> } ) };
		for my $f ( @set ) {
			if ( not $seen{ $f->nfos }++ ) {				
				push @fossils, $f;
			}
		}
	}
	
	# here we have to ensure that the names of fossils
	# match up with existing names on the tree. XXX for
	# now we ditch fossils that don't map, and we ditch
	# node names that have no fossil
	$log->info("going to map fossils to tree");
	my %node_by_name;
	for my $node ( @{ $tree->get_internals } ) {
		if ( my $name = $node->get_name ) {
		
			# if nhx:D (duplication event) is N ("No", though
			# NHX actually specifies "F", false), it's
			# a speciation event, which we can calibrate
			my $nhx_D = $node->get_meta_object('nhx:D');
			if ( $nhx_D eq 'N' ) {
				$node_by_name{$name} = [] if not $node_by_name{$name};
				push @{ $node_by_name{$name} }, $node;
				$log->debug("added speciation node $name");
			}
			else {
				$log->debug("found duplication node (nhx:D=$nhx_D) for $name");
			}
		}
	}
	
	# here we must map the same fossil to all speciation 
	# events (in different gene copies) that it corresponds with
	my @cleaned_fossils;
	FOSSIL: for my $f ( @fossils ) {
		my $name = $f->calibrated_taxon;
		
		# FIXME: if it's a stem fossil we actually should 
		# apply a calibration to the parent. But if that's
		# a duplication event we'd have to traverse deeper
		# and it'll all get a bit messy. Skip for now.
		if ( $f->crown_vs_stem ne 'crown' ) {
			$log->warn("Fossil '$name' is not a crown fossil, skipping");
			next FOSSIL;
		}
		
		# We have a crown fossil, and one or more speciation
		# nodes that correspond with it.
		if ( my $nodes = $node_by_name{$name} ) {
		
			# add distinct suffix to each instance, both
			# in tree and on clones of the fossil
			my $i = 1;
			for my $node ( @{ $nodes } ) {
				my $instance = $name . '_' . $i++;
				$node->set_name($instance);
				my $clone = Load(Dump($f));
				$clone->calibrated_taxon($instance);
				push @cleaned_fossils, $clone;
			}
			$log->info("successfully mapped fossil $name");
		}
		else {
			$log->debug("couldn't map fossil $name");
		}
	}
	return @cleaned_fossils;
}

# run r8s
sub run_r8s {
	my ($project,$nchar,@fossils) = @_;
	
	# first write to memory (so we can log), 
	# then to a temporary file
	my $r8s_commands;
	my ( $fh, $filename ) = tempfile();	
	my $tmpl = Template->new;
	my $date = localtime();
	$tmpl->process(
		$template,
		{
			'project' => $project,
			'nchar'   => $nchar,
			'fossils' => \@fossils,
			'script'  => $0,
			'date'    => $date,
		},
		\$r8s_commands
	);
	$log->debug("\n".$r8s_commands."\n");
	print $fh $r8s_commands;
	
	# capture output for parsing later
	$log->info("going to run r8s, this may take a while. watch log: $filename.log");
	system("$exe -b -f $filename > $filename.log");
	open my $infh, '<', "$filename.log";
	my $result = do { local $/; <$infh> };
	unlink $filename, "$filename.log";
	$log->debug("\n".$result."\n");
	return $result;
}

sub parse_result {
	my $result  = shift;
	my $counter = 0;
	my $passed;
	my @type = qw(ratogram chronogram phylogram);
	LINE: for my $line ( split /\n/, $result ) {
	
		# check if we passed, set flag
		if ( $line =~ /^\s*PASSED\s*$/ ) {
			$log->info("analysis passed OK - reading trees");
			$passed++;
		}
		
		# if we have a passed flag and a tree description, 
		# it's first the ratogram, then the chronogram
		if ( $passed and $line =~ /tree Tree\d+ = (\(.+;)/ ) {
			my $tree = $1;
			my $type = $type[$counter];
			
			# write to file
			my $outfile = $intree;
			$outfile =~ s/.+\///;
			$outfile = join '', $ratodir, '/', $outfile, '.', $type;
			open my $fh, '>', $outfile or die $!;
			print $fh $tree;
			$log->info("$type written to $outfile");
			$counter++;
		}
	}
	$log->error("*** analysis failure in $intree") if not $passed;
}

