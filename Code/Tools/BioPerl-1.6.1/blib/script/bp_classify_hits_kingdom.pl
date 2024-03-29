#!/opt/local/bin/perl5.12 -w

eval 'exec /opt/local/bin/perl5.12 -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell
# $Id: classify_hits_kingdom.PLS 15088 2008-12-04 02:49:09Z bosborne $

=head1 NAME

classify_hits_kingdom - classify BLAST hits by taxonomic kingdom

=head2 USAGE

classify_hits_kingdom [-i tab_file] [-i second_BLAST_file] [-e evalue_cutoff]
                      [-t dir_where_TAXONOMY_files_are] [-g gi2taxid] 
                      [-z PATH_TO_zcat] [-v]

=head2 DESCRIPTION

Will print out the taxonomic distribution (at the kingdom level) for a
set of hits against the NR database.  This script assumes you've done
a search against the protein database, you'll have to make minor
changes in the gi_taxid part to point to the gi_taxid_nuc.dump file.

This expects BLAST files in tabbed -m9 or -m8 format.  Output with -m
8 or use blast2table.pl to convert (or fastam9_to_table.PLS if using
FASTA).

  Input values:
   -t/--taxonomy  directory where the taxonomy .dmp files are (from NCBI)
   -g/--gi        Location of gi_taxid_prot.dmp (or gi_taxid_nucl.dmp if 
                  the search was against a NT db)
   -i/--in        The name of the tab delimited -m8/-m9 output files to 
                  process.

    -e/--evalue   Provide an E-value cutoff for hits to be considered
    -z/--zcat     Path to the 'zcat' executable, can also be 'gunzip -c'
                  if no zcat on your system.
   Flags
    -v/--verbose  To turn on verbose messages
    -h/--help     Display this helpful information

This is intended to be useful starting script, but users may want to
customize the output and parameters.  Note that I am summarizing the
kingdoms here and Eukaryota not falling into Metazoa, Viridiplantae,
or Fungi gets grouped into the general superkingdom Eukaryota. for
simplicity.  There are comments in the code directing you to where
changes can be made if you wanted to display hits by phylum for
example.  Note that you must wipe out the cache file 'gi2class' that
is craeed in your directory after making these changes.

=head2 AUTHOR

Jason Stajich jason_at_bioperl_dot_org

=cut

use strict;
use Bio::DB::Taxonomy;
use DB_File;
use Env;
use File::Spec;
use vars qw($SEP);
my $DEBUG = 0;
use Getopt::Long;
$SEP = '_';

my $evalue_filter = 1e-3;
my @files;
my $zcat = 'zcat'; # or gunzip -c 
my $prefix = File::Spec->catfile($HOME,'taxonomy');
my $gi2taxidfile = "$prefix/gi_taxid_prot.dmp.gz";
GetOptions(
	   'v|verbose|debug' => \$DEBUG,
	   'z|zcat:s'    => \$zcat,
	   'i|in:s'      => \@files,
	   'e|evalue:f'  => \$evalue_filter,
	   't|taxonomy:s'=> \$prefix,
	   'g|gi|gi2taxid:s' => \$gi2taxidfile,
	   'h|help'      => sub { system('perldoc', $0);
				  exit },
	   );

# insure idx location is created
mkdir(File::Spec->catfile($prefix,'idx')) 
    unless -d File::Spec->catfile($prefix,'idx');

# these files came from ftp://ftp.ncbi.nih.gov/pub/taxonomy
my $taxdb = Bio::DB::Taxonomy->new
    (-source => 'flatfile',
     -directory => File::Spec->catfile
     ($prefix, 'idx'), 
     -nodesfile => File::Spec->catfile($prefix,'nodes.dmp'),
     -namesfile => File::Spec->catfile($prefix,'names.dmp')
     );
my %query;

my (%taxid4gi,%gi2node);
my $dbh = tie(%gi2node, 'DB_File', 'gi2class');


my $giidxfile = File::Spec->catfile($prefix,'idx','gi2taxid');
my $done = -e $giidxfile;
my $dbh2 = tie(%taxid4gi, 'DB_File', $giidxfile);

if( ! $done ) {
    my $fh;
    # this file came from ftp://ftp.ncbi.nih.gov/pub/taxonomy
    # I'm interested in protein hits therefor _prot file.
    if( $gi2taxidfile =~ /\.gz$/ ) {
	open($fh, "$zcat $gi2taxidfile |" ) || die "$zcat $gi2taxidfile: $!";
    } else {
	open($fh, $gi2taxidfile ) || die $!;
    }
    my $i = 0;
    while(<$fh>) {
	my ($gi,$taxid) = split;
	$taxid4gi{$gi} = $taxid;
	$i++;
	unless( $DEBUG && $i % 100000  ) {
	    warn "$i\n";
	}
    }
    $dbh2->sync;
}

for my $file ( @files ) {
    warn("$file\n");
    my $gz;
    if( $file =~ /\.gz$/) {
	$gz = 1;
    }
    my ($spname) = split(/\./,$file); 
    my ($fh,$i);
    if( $gz ) {
	open($fh, "$zcat $file |")  || die "$zcat $file: $!";
    } else {
	open($fh, $file) || die "$file: $!";
    }
    while(<$fh>) {
	my ($qname,$hname,$pid,$qaln,$mismatch,$gaps,
	    $qstart,$qend,$hstart,$hend,
	    $evalue,$bits,$score) = split(/\t/,$_);	
	next if( $evalue > $evalue );
	if( ! exists $query{$spname}->{$qname} ) {
	    $query{$spname}->{$qname} = {};
	}
	
	if( $hname =~ /gi\|(\d+)/) {		
	    my $gi = $1;	    
	    if( ! $gi2node{$gi} ){ # see if we cached the results from before
		my $taxid = $taxid4gi{$gi};
		if( ! $taxid ) {
		    warn("no taxid for $gi\n");
		    next;
		}
		my $node = $taxdb->get_Taxonomy_Node($taxid);
		if( ! $node ) {
		    warn("cannot find node for gi=$gi ($hname) (taxid=$taxid)\n");
		    next;
		}
		my $parent = $taxdb->get_Taxonomy_Node($node->parent_id);

		# THIS IS WHERE THE KINGDOM DECISION IS MADE
		# DON'T FORGET TO WIPE OUT YOUR CACHE FILE
		# gi2class after you make changes here
		while( defined $parent && $parent->node_name ne 'root' ) { 
		    # this is walking up the taxonomy hierarchy
		    # can be a little slow, but works...
		    #warn( "\t",$parent->rank, " ", $parent->node_name, "\n");
		    # deal with Eubacteria, Archea separate from 
		    # Metazoa, Fungi, Viriplantae differently
		    # (everything else Eukaryotic goes in Eukaryota)
		    if( $parent->rank eq 'kingdom') {
			# caching in ... 
			($gi2node{$gi}) = $parent->node_name;
			last;
		    } elsif( $parent->rank eq 'superkingdom' ) {
			# caching in ... 
			($gi2node{$gi}) = $parent->node_name;
			$gi2node{$gi} =~ s/ \<(bacteria|archaea)\>//g;
			last;
		    }
		    $parent = $taxdb->get_Taxonomy_Node($parent->parent_id);
		}		
		$dbh->sync;
	    }
	    my ($kingdom) = $gi2node{$gi};
#		warn("$gi2node{$gi}\n");
	    unless( defined $kingdom && length($kingdom) ) {
#		    warn("no kingdom for $hname\n");
	    } else {
		$query{$spname}->{$qname}->{$kingdom}++;		
	    }	
	} else {
	    warn("no GI in $hname\n");
	}
    }
    last if ( $DEBUG && $i++ > 10000);
}

# print out the taxonomic distribution
while( my ($sp,$d) = each %query ) {
    my $total = scalar keys %$d;
    print "$sp total=$total\n";
    my %seen;
    for my $v ( values %$d ) {
	my $tag = join(",",sort keys %$v );
	$seen{$tag}++;
    }
    for my $t ( sort { $seen{$a} <=> $seen{$b} } keys %seen ) {
	printf " %-20s\t%d\t%.2f%%\n",
	$t,$seen{$t}, 100 * $seen{$t} / $total;
    }
    print "\n\n";
}
