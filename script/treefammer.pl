#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Archive::Tar;
use Bio::Phylo::IO qw'parse_tree unparse';
use Bio::Phylo::Util::Logger ':levels';
use Bio::Phylo::Util::CONSTANT ':namespaces';

# process command line arguments
my $infile;
my $treefam;
my $verbosity = WARN;
GetOptions(
	'infile=s'  => \$infile,
	'treefam=s' => \$treefam,
	'verbose+'  => \$verbosity,
);

# instantiate helper objects
my $tar = Archive::Tar->new;
my $log = Bio::Phylo::Util::Logger->new(
	'-level' => $verbosity,
	'-class' => 'main',
);

# open the archive
$log->info("going to read $infile");
my $iter = $tar->iter( $infile, { 'filter' => qr/^.+(?:nhx|aln)\.emf$/ } );

# iterate over file entries
while( my $file = $iter->() ) {
	my $name = $file->name;
	$name =~ s/.+\///;
	
	# skip anything other than NHX trees and DNA sequence alignments (FASTA)
	next unless $name and $name =~ qr/^.+(?:nhx|cds)/;
	$log->info("processing $name");
	my $content;
	
	# deal with tree descriptions
	if ( $name =~ /nhx/ ) {
		my %map;
		for my $line ( split /\n/, $file->get_content ) {
		
			# capture accession to taxon mapping
			if ( $line =~ /^SEQ/ ) {
				my ( $seq, $taxon, $ac ) = split /\s+/, $line;
				$map{$ac} = $taxon;
			}
			elsif ( $line =~ /^(\(.+;)/ ) {
			
				# parse NHX tree
				my $nhx = $1;
				my $tree = parse_tree(
					'-format' => 'nhx',
					'-string' => $nhx,
				);

				# remap tip labels
				my ( %seen, @prune );
				$tree->set_namespaces( 'nhx' => _NS_NHX_ );
				$tree->visit(sub{
					my $node = shift;
					if ( $node->is_terminal ) {
						my $acc = $node->get_name;
						if ( my $name = $map{$acc} ) {
							$name = ucfirst $name;
							$node->set_name( $name . '_' . ++$seen{$name} );
							$node->set_meta_object( 'nhx:AC' => $acc );
						}
						else {
							$log->error("no mapping for accession '$acc' in $infile - will prune");
							push @prune, $node;
						}
					}
				});
				$tree->prune_tips(\@prune) if scalar @prune;
				
				# write back to NHX
				$content = unparse( 
					'-format'     => 'nhx', 
					'-phylo'      => $tree, 
					'-nodelabels' => 1 
				);
			}
		}		
	}
	else {
		$content = $file->get_content;
	}
	open my $fh, '>', "${treefam}/${name}" or die $!;
	print $fh $content;
}