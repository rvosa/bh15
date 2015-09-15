#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Archive::Tar;
use Bio::Phylo::Util::Logger ':levels';

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
	next unless $name and $name =~ qr/^.+(?:nhx|cds)/;
	$log->info("processing $name");
	my $content;
	if ( $name =~ /nhx/ ) {
		($content) = grep { /^\(/ } split /\n/, $file->get_content;		
	}
	else {
		$content = $file->get_content;
	}
	open my $fh, '>', "${treefam}/${name}" or die $!;
	print $fh $content;
}