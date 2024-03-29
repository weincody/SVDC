# -*-Perl-*- Test Harness script for Bioperl
# $Id: Flat.t 15178 2008-12-16 12:48:21Z sendu $

use strict;

BEGIN { 
	use lib '.';
	use Bio::Root::Test;
	
	test_begin(-tests => 24,
				  -requires_module => 'DB_File');
	
	use_ok('Bio::DB::Flat');
}

my $verbose = test_debug();

# First of all we need to create an flat db

my $tmpdir = test_output_dir();

my $db = Bio::DB::Flat->new(-directory  => $tmpdir,
                            -index      => 'bdb',
									 -dbname     => 'mydb',
									 -format     => 'fasta',
									 -verbose    => $verbose,
									 -write_flag => 1 );
ok($db);
my $dir = test_input_file('AAC12660.fa');
my $result = $db->build_index(glob($dir));
ok($result);

# Now let's get the sequence out again
my $seq = $db->get_Seq_by_id('AAC12660');
ok($seq);
is($seq->length,504);
undef $db;

$db = Bio::DB::Flat->new(-directory  => $tmpdir,
                         -index      => 'bdb',
                         -format     => 'embl',
								 -dbname     => 'myembl',
                         -verbose    => $verbose,
                         -write_flag => 1 );

$dir= test_input_file('factor7.embl');

$result = $db->build_index(glob($dir));
ok($result);
$seq = $db->get_Seq_by_id('HSCFVII');
ok($seq);
is($seq->length,12850);

# deal with wantarray conditions
$seq = $db->get_Seq_by_acc('J02933');
ok($seq && ref($seq));
is($seq->length,12850);


undef $db;

$db = Bio::DB::Flat->new(-directory  => $tmpdir,
			 -index      => 'binarysearch',
			 -format     => 'fasta',
			 -dbname     => 'mybinfa',
			 -verbose    => $verbose,
			 -write_flag => 1
			 );

$dir= test_input_file('dbfa', '1.fa');
$result = $db->build_index($dir);
ok($result);
$seq = $db->get_Seq_by_id('AW057119');
ok($seq);
is($seq->length,808);
undef $db;

$db = Bio::DB::Flat->new(-directory  => $tmpdir,
			 -index      => 'binarysearch',
			 -format     => 'swiss',
			 -dbname     => 'mybinswiss',
			 -verbose    => $verbose,
			 -write_flag => 1
			 );
$dir= test_input_file('swiss.dat');
$result = $db->build_index($dir);

ok($result);
$seq = $db->get_Seq_by_id('ACON_CAEEL');
ok($seq);
is($seq->length,788);

$seq = $db->get_Seq_by_id('ACON_CAEEL');
ok($seq && ref($seq));

undef $db;

$db = Bio::DB::Flat->new(-directory  => $tmpdir,
                         -index      => 'binarysearch',
                         -format     => 'fasta',
								 -dbname     => 'myfasta',
                         -verbose    => $verbose,
                         -write_flag => 1 );

$dir = test_input_file('tmp.fst');
$result = $db->build_index(glob($dir));
ok($result);
$seq = $db->get_Seq_by_id('TEST00004');
is($seq->length,98);

undef $db;

$db = Bio::DB::Flat->new(-directory  => $tmpdir,
                         -index      => 'bdb',
                         -format     => 'fasta',
								 -dbname     => 'mybfasta',
                         -verbose    => $verbose,
                         -write_flag => 1 );

$dir = test_input_file('tmp.fst');
$result = $db->build_index(glob($dir));
ok($result);
for my $id ( qw(TEST00001 TEST00002 TEST00003 TEST00004) ) {
	$seq = $db->get_Seq_by_id($id);
	is($seq->length,98);
}
