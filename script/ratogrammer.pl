#!/usr/bin/perl
use strict;
use warnings;
use YAML;
use Getopt::Long;
use Bio::Phylo::IO 'parse';
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::Util::CONSTANT ':objecttypes';
use Bio::Phylo::PhyLoTA::Service::CalibrationService;

# lookup to go from ensembl ID prefixes to taxon names
my %lookup;
while(<DATA>) {
	chomp;
	my ( $code, $species ) = split /\t/, $_;
	$lookup{$code} = $species;
}

# process command line arguments
my $fcdir;            # directory for cached calibrations
my $ratodir;          # directory to write ratogram
my $intree;           # input tree
my $exe = 'r8s';      # r8s executable
my $verbosity = WARN; # log level
GetOptions(
	'fcdir=s'   => \$fcdir,
	'ratodir=s' => \$ratodir,
	'intree=s'  => \$intree,
	'exe=s'     => \$exe,
	'verbose+'  => \$verbosity,
);

# instantiate helper objects
my $cs = Bio::Phylo::PhyLoTA::Service::CalibrationService->new;
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main',
);
my $project = parse(
	'-format'     => 'newick',
	'-file'       => $intree,
	'-as_project' => 1,
);

# fetch tree from project
my ($tree) = @{ $project->get_items(_TREE_) };

# relabel terminal nodes
for my $tip ( @{ $tree->get_terminals } ) {
	my $name = $tip->get_name;
	$name =~ s/P\d+$//;
	if ( $lookup{$name} ) {
		$tip->set_name($lookup{$name});
	}
	else {
		$log->error("no mapping for prefix $name");
	}
}

# fetch calibration points for nodes we needed to do the
# relabeling of the tips first because the $cs uses this,
# and we do a pre-order traversal because TreeFam labels
# multiple, nested nodes with the same higher taxon label.
# we want the deepest of the nested nodes, which we visit
# first in this traversal.
my %seen;
$tree->visit_depth_first(
	'-pre' => sub {
		my $node = shift;
		my $name = $node->get_name;
		my $file = "${fcdir}/${name}.yml";
		if ( $node->is_internal ) {
		
			# already processed this taxon, either during
			# this traversal or during a previous run			
			if ( -e $file ) {
				$log->info("already fetched fossils for $name");
			}
			else {
			
				# fetch fossils
				$log->info("fetching fossils for $name");
				my @records = $cs->fetch_fossil_dates( $name );
				open my $fh, '>', $file or die $!;
				binmode( $fh, ':utf8' );
				print $fh Dump(\@records);
				close $fh;
				
				# remove nested identical node label
				if ( $seen{$name}++ ) {
					$node->set_name;
				}
			}
		}
	}
);

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