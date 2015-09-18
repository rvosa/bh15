#!/usr/bin/perl
package Arg;
use strict;
use warnings;
use Carp;
our $AUTOLOAD;

sub new {
	my ( $package, %args ) = @_;
	my $self = {
		'chrono' => undef, # dup node in chronogram, will traverse
		'rato'   => undef, # idem in ratogram
		'phylo'  => undef, # idem in phylogram
		'taxon'  => undef, # static
		'dist'   => undef, # static
		'hash'   => undef, # static
		'ntips'  => undef, # computed in first iteration
		'height' => undef, # computed in first iteration	
		'pdist'  => undef, # phylogram distance
	};
	bless $self, $package;
	for my $key ( keys %args ) {
		$self->$key( $args{$key} );
	}
	return $self;
}

sub AUTOLOAD {
	my ( $self, $arg ) = @_;
	my $field = $AUTOLOAD;
	$field =~ s/.+://;
	return if $field =~ /^[A-Z]+$/;
	if ( exists $self->{$field} ) {
		if ( defined $arg ) {
			$self->{$field} = $arg;
		}
		return $self->{$field};
	}
	else {
		croak "No field $field!";
	}
}

package main;
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use List::Util 'sum';
use Digest::MD5 'md5_base64';
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::Util::Logger ':levels';

# process command line arguments
my $chronogram;
my $ratogram;
my $phylogram;
my $verbosity = WARN;
GetOptions(
	'chronogram=s' => \$chronogram,
	'ratogram=s'   => \$ratogram,
	'phylogram=s'  => \$phylogram,
	'verbose+'     => \$verbosity,
);

# instantiate helper objects
my $log = Bio::Phylo::Util::Logger->new(
	'-level'  => $verbosity,
	'-class'  => 'main',
);

# parse trees
my @trees;
for my $file ( $chronogram, $ratogram, $phylogram ) {
	$log->info("parsing $file");
	push @trees, parse_tree( '-format' => 'newick', '-file' => $file )->ladderize;
	$log->info("tree has ".scalar(@{$trees[-1]->get_terminals})." tips");
}
my ( $ct, $rt, $pt ) = @trees;

# index the trees
my $map = {};
index_tree($ct,$map);
index_tree($rt,$map);
index_tree($pt,$map);
$log->info("indexed trees");

# parse TreeFam family ID: strip path and extension(s)
my $TFID = $chronogram;
$TFID =~ s/.+\///;
$TFID =~ s/\..+//;
$log->info("tree family ID: $TFID");

# iterate over duplication nodes
print "treefam\ttaxon\tdistance\trate\thash\tntips\theight\n";
for my $hash ( keys %$map ) {
	my ( $cn, $rn, $pn ) = @{ $map->{$hash} };
	my $arg = Arg->new(
		'chrono' => $cn, # dup node in chronogram, will traverse
		'rato'   => $rn, # idem in ratogram
		'phylo'  => $pn, # idem in phylogram
		'taxon'  => $cn->get_name, # static
		'dist'   => $cn->get_generic('dist'), # static
		'hash'   => $hash, # static
		'pdist'  => $pn->get_generic('dist'),
	);
	traverse( $arg );
}

sub traverse {
	my ( $s, $ntips, $height ) = @_;
	my @cc = @{ $s->chrono->get_children };
	my @rc = @{ $s->rato->get_children };
	my @pc = @{ $s->phylo->get_children };
	for my $i ( 0 .. $#cc ) {
		$s->chrono( $cc[$i] );
		$s->rato(   $rc[$i] );
		$s->phylo(  $pc[$i] );
		my ( $first_ntips, $first_height );
	
		# compute number of tips first time, then keep carrying over
		if ( not defined $ntips ) {
			$first_ntips = scalar @{ $s->chrono->get_generic('tips') };
			$log->info($s->taxon . ' => ' . $first_ntips);
		}	
		
		# compute average tip height first time, then keep carrying over
		if ( not defined $height ) {
			my @heights = map { $_->get_generic('dist') } @{ $s->phylo->get_terminals };
			$first_height = ( sum(@heights) / scalar(@heights) ) - $s->pdist;
		}
		
		# compute distance from dup node and rate on focal branch
		my $dist = $s->chrono->get_generic('dist') - $s->dist;
		my $rate = $s->rato->get_branch_length;
		
		# print record
		print join("\t",
			$TFID,      # GLOBAL: treefam family id (stem of input files)
			$s->taxon,  # STATIC: taxon label on duplication node
			$dist,      # VAR: distance (MYA) of focal node to duplication node
			$rate,      # VAR: substitution rate on the focal branch
			$s->hash,   # STATIC: dup node identifier (MD5 of sorted tip labels)
			$ntips  || $first_ntips,  # STATIC: number of tips subtended by children of dup node
			$height || $first_height, # STATIC: average tip height subtended by children of dup node
		), "\n";
		
		# stop traversing when reaching another duplication node
		return if $s->chrono->get_name and $s->chrono->get_name !~ /_\d+$/;	
				
		# traverse deeper
		traverse( $s, $ntips, $height );
	}
}

# pre-compute the node distance from root
sub pre_order_processing {
	my $node = shift;
	if ( my $parent = $node->get_parent ) {
		my $dist = $parent->get_generic('dist') + $node->get_branch_length;
		$node->set_generic( 'dist' => $dist );
	}
	else {
		$node->set_generic( 'dist' => 0 );
	}
}

sub index_tree {
	my ( $tree, $map ) = @_;
	$tree->visit_depth_first(
		'-pre'  => \&pre_order_processing,
		'-post' => sub {
			my $node = shift;
			my $name = $node->get_name;
			if ( $node->is_terminal ) {
			
				# start the tips array
				$node->set_generic( 'tips' => [ $name ] );
			}
			else {
			
				# extend the tips array, take an MD5 hash
				my %tips;
				for my $child ( @{ $node->get_children } ) {
					$tips{$_}++ for @{ $child->get_generic('tips') };
				}
				my $hash = md5_base64( join ',', sort { $a cmp $b } keys %tips );
				$node->set_generic( 'tips' => [ keys %tips ] );
				$node->set_generic( 'hash' => $hash );	
				
				# create mapping, if duplication node
				if ( $name and $name !~ /_\d+$/ ) {
					$map->{$hash} = [] if not $map->{$hash};
					push @{ $map->{$hash} }, $node;
				}			
			}
		}
	);
}