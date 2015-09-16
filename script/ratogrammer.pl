#!/usr/bin/perl
use strict;
use warnings;
use YAML;
use Template;
use Getopt::Long;
use File::Temp 'tempfile';
use Bio::Phylo::IO 'parse';
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::Util::CONSTANT ':objecttypes';
use Bio::Phylo::PhyLoTA::Service::CalibrationService;

# process command line arguments
my $fcdir;            # directory for cached calibrations
my $intree;           # input tree
my $exe = 'r8s';      # r8s executable
my $verbosity = WARN; # log level
my $template;         # template to generate r8s infile
my $ratodir;          # output directory for ratograms
GetOptions(
	'fcdir=s'    => \$fcdir,
	'intree=s'   => \$intree,
	'exe=s'      => \$exe,
	'template=s' => \$template,
	'verbose+'   => \$verbosity,
	'ratodir=s'  => \$ratodir,
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
	
	# lookup to go from ensembl ID prefixes to taxon names
	my %lookup;
	while(<DATA>) {
		chomp;
		my ( $code, $species ) = split /\t/, $_;
		$lookup{$code} = $species;
	}	

	# parse the tree
	$log->info("going to read $intree");
	my $project = parse(
		'-format'     => 'newick',
		'-file'       => $intree,
		'-as_project' => 1,
	);

	# fetch tree from project
	my ($tree) = @{ $project->get_items(_TREE_) };

	# relabel terminal nodes
	my ( @prune, %seen );
	for my $tip ( @{ $tree->get_terminals } ) {
		my $name = $tip->get_name;
		$name =~ s/P\d+$//;
		if ( $lookup{$name} ) {
			$tip->set_name($lookup{$name} . '.' . ++$seen{$name});
		}
		else {
			$log->error("no mapping for prefix $name - will prune (!!!)");
			push @prune, $tip;
		}
	}
	$tree->prune_tips(\@prune);
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
	{
		my %seen;
		$tree->visit_depth_first(
			'-pre' => sub {
				my $node = shift;
				my $name = $node->get_name;
				my $file = "${fcdir}/${name}.yml";
				if ( $node->is_internal ) {

					# store fossil file name for later
					push @files, $file;
		
					# already processed this taxon, either during
					# this traversal or during a previous run			
					if ( -e $file ) {
						$log->info("already fetched fossils for $name");
					}
					else {
			
			
						# fetch fossils
						$log->info("fetching fossils for $name");
						my @records = $cs->fetch_fossil_dates( $name );
						if ( scalar @records ) {
							open my $fh, '>', $file or die $!;
							binmode( $fh, ':utf8' );
							print $fh Dump(\@records);
							close $fh;
						}
						else {
							$log->info("no fossils for $name");
							pop @files;
						}
					}
					
					# remove nested identical node label
					if ( $seen{$name}++ ) {
						$log->warn("already seen $name in ancestor");
						$node->set_name('');
					}					
				}
			}
		);
	}
	return @files;
}

# read fossils from file. XXX should test topology. doesn't, yet.
sub get_fossils {
	my @files = @_;
	my ( @fossils, %seen );
	for my $file ( @files ) {
		next if $seen{$file}++;
		$log->info("going to read fossils in $file");
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
	my %node_by_name = map { $_->get_name => $_ } @{ $tree->get_internals };
	my @cleaned_fossils;
	for my $f ( @fossils ) {
		my $name = $f->calibrated_taxon;
		if ( my $node = $node_by_name{$name} ) {
			push @cleaned_fossils, $f;
			delete $node_by_name{$name};
			$log->info("successfully mapped fossil $name");
		}
		else {
			$log->warn("couldn't map fossil $name");
		}
	}
	$_->set_name('') for values %node_by_name;
	return @cleaned_fossils;
}

# run r8s
sub run_r8s {
	my ($project,$nchar,@fossils) = @_;
	my ( $fh, $filename ) = tempfile();
	my $r8s_commands;
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
	$log->info("going to run r8s, this may take a while");
	my $result = `$exe -b -f $filename`;
	$log->debug("\n".$result."\n");
	return $result;
}

sub parse_result {
	my $result = shift;
	my $passed;
	LINE: for my $line ( split /\n/, $result ) {
	
		# check if we passed, set flag
		if ( $line =~ /^\s*PASSED\s*$/ ) {
			$log->info("analysis passed OK - reading ratogram");
			$passed++;
		}
		
		# if we have a passed flag and a tree description, it's the ratogram
		if ( $passed and $line =~ /tree Tree\d+ = (\(.+;)/ ) {
			my $ratogram = $1;
			
			# write to file
			my $outfile = $intree;
			$outfile =~ s/.+\///;
			$outfile = $ratodir . '/' . $outfile;
			open my $fh, '>', $outfile or die $!;
			print $fh $ratogram;
			$log->info("ratogram written to $outfile");
		}
	}
	$log->error("*** analysis failure in $intree") if not $passed;
}

__DATA__
ENSAME	Ailuropoda melanoleuca
ENSAPL	Anas platyrhynchos
ENSACA	Anolis carolinensis
ENSAMX	Astyanax mexicanus
ENSBTA	Bos taurus
ENSCEL	Caenorhabditis elegans
ENSCJA	Callithrix jacchus
ENSCAF	Canis lupus familiaris
ENSCPO	Cavia porcellus
ENSCSA	Chlorocebus sabaeus
ENSCHO	Choloepus hoffmanni
ENSCIN	Ciona intestinalis
ENSCSAV	Ciona savignyi
ENSDAR	Danio rerio
ENSDNO	Dasypus novemcinctus
ENSDOR	Dipodomys ordii
FB	Drosophila melanogaster
ENSETE	Echinops telfairi
ENSECA	Equus caballus
ENSEEU	Erinaceus europaeus
ENSFCA	Felis catus
ENSFAL	Ficedula albicollis
ENSGMO	Gadus morhua
ENSGAL	Gallus gallus
ENSGAC	Gasterosteus aculeatus
ENSGGO	Gorilla gorilla gorilla
ENS	Homo sapiens
ENSSTO	Ictidomys tridecemlineatus
ENSLAC	Latimeria chalumnae
ENSLOC	Lepisosteus oculatus
ENSLAF	Loxodonta africana
ENSMMU	Macaca mulatta
ENSMEU	Macropus eugenii
ENSMGA	Meleagris gallopavo
ENSMIC	Microcebus murinus
ENSMOD	Monodelphis domestica
ENSMUS	Mus musculus
ENSMPU	Mustela putorius furo
ENSMLU	Myotis lucifugus
ENSNLE	Nomascus leucogenys
ENSOPR	Ochotona princeps
ENSONI	Oreochromis niloticus
ENSOAN	Ornithorhynchus anatinus
ENSOCU	Oryctolagus cuniculus
ENSORL	Oryzias latipes
ENSOGA	Otolemur garnettii
ENSOAR	Ovis aries
ENSPTR	Pan troglodytes
ENSPAN	Papio anubis
ENSPSI	Pelodiscus sinensis
ENSPMA	Petromyzon marinus
ENSPFO	Poecilia formosa
ENSPPY	Pongo abelii
ENSPCA	Procavia capensis
ENSPVA	Pteropus vampyrus
ENSRNO	Rattus norvegicus
ENSSCE	Saccharomyces cerevisiae
ENSSHA	Sarcophilus harrisii
ENSSAR	Sorex araneus
ENSSSC	Sus scrofa
ENSTGU	Taeniopygia guttata
ENSTRU	Takifugu rubripes
ENSTSY	Tarsius syrichta
ENSTNI	Tetraodon nigroviridis
ENSTBE	Tupaia belangeri
ENSTTR	Tursiops truncatus
ENSVPA	Vicugna pacos
ENSXET	Xenopus tropicalis
ENSXMA	Xiphophorus maculatus