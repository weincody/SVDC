# -*-Perl-*- Test Harness script for Bioperl
# $Id: alf.t 15383 2009-01-16 16:45:08Z cjfields $

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests               => 8,
			   -requires_module     => 'Bio::SeqIO::staden::read');
	
	use_ok('Bio::SeqIO::alf');
}

my $verbose = test_debug();

TODO: {
	my $format = 'alf';
	todo_skip "No tests for $format format -- no sample file to test against", 7, if 1;

	my $seqio_obj = Bio::SeqIO->new(-file   => test_input_file("test.$format"),
							        -format => $format);
	
	isa_ok($seqio_obj, 'Bio::SeqIO');
	
	my @methods = qw(next_seq write_seq);
	foreach my $method (@methods) {
		can_ok($seqio_obj, $method) || 
			diag "$method method not implemented for $format";	
	}
	
	# checking the first sequence object
	my $seq_obj = $seqio_obj->next_seq();
	isa_ok($seq_obj, 'Bio::Seq::Quality');
	my %expected = ('seq'         => '' .
					'length'      => '',
					'primary_id'  => '',
					'description' => qr(),
				   );
	is   ($seq_obj->seq(),         $expected{'seq'},         'sequence');
	is   ($seq_obj->length(),      $expected{'length'},      'length');
	is   ($seq_obj->primary_id(),  $expected{'primary_id'},  'primary_id');
	like ($seq_obj->description(), $expected{'description'}, 'description');
}