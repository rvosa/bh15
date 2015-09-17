#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Digest::MD5 'md5_base64';
use Bio::Phylo::IO 'parse_tree';
use Bio::Phylo::Util::Logger ':levels';

# process command line arguments
my $chronogram;
my $ratogram;
my $verbosity = WARN;
GetOptions(
	'chronogram=s' => \$chronogram,
	'ratogram=s'   => \$ratogram,
	'verbose+'     => \$verbosity,
);

# instantiate helper objects
my $log = Bio::Phylo::Util::Logger->new(
	'-level'  => $verbosity,
	'-class'  => 'main',
);
my $ct = parse_tree(
	'-format' => 'newick',
	'-file'   => $chronogram,
)->ladderize;
my $rt = parse_tree(
	'-format' => 'newick',
	'-file'   => $ratogram,
)->ladderize;

# index the trees
my $map = {};
index_tree($ct,$map);
index_tree($rt,$map);

# iterate over duplication nodes
print "taxon\tdistance\trate\n";
for my $hash ( keys %$map ) {
	my ( $cn, $rn ) = @{ $map->{$hash} };
	traverse( $cn, $rn, $cn->get_name, $cn->get_generic('dist') );
}

sub traverse {
	my ( $cn, $rn, $taxon, $root_dist ) = @_;
	my @cc = @{ $cn->get_children };
	my @rc = @{ $rn->get_children };
	for my $i ( 0 .. $#cc ) {
	
		# stop traversing when reaching another duplication node
		if ( $cc[$i]->get_name and $cc[$i]->get_name !~ /_\d+$/ ) {
			return;
		}
		
		# compute and print rate by dist
		my $dist = $cc[$i]->get_generic('dist') - $root_dist;
		my $rate = $rc[$i]->get_branch_length;
		print join("\t",$taxon,$dist,$rate), "\n";
		
		# traverse deeper
		traverse( $cc[$i], $rc[$i], $taxon, $root_dist );
	}
}

sub index_tree {
	my ( $tree, $map ) = @_;
	$tree->visit_depth_first(
		'-pre' => sub {
		
			# pre-compute the node distance from root
			my $node = shift;
			if ( my $parent = $node->get_parent ) {
				my $dist = $parent->get_generic('dist') + $node->get_branch_length;
				$node->set_generic( 'dist' => $dist );
			}
			else {
				$node->set_generic( 'dist' => 0 );
			}
		},
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