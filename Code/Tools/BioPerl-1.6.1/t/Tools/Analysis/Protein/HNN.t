# -*-Perl-*- Test Harness script for Bioperl
# $Id: HNN.t 15112 2008-12-08 18:12:38Z sendu $

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 13,
			   -requires_modules => [qw(IO::String LWP::UserAgent)],
			   -requires_networking => 1);
	
	use_ok("Bio::Seq");
	use_ok("Bio::Tools::Analysis::Protein::HNN");
}

my $seq = Bio::Seq->new(-seq => 'MSADQRWRQDSQDSFGDSFDGDPPPPPPPPFGDSFGDGFSDRSRQDQRS',
                        -display_id => 'test2');
ok my $tool = Bio::Tools::Analysis::Protein::HNN->new(-seq=>$seq->primary_seq);

SKIP: {
	ok $tool->run();
	skip "Skipping tests since we got terminated by a server error", 9 if $tool->status eq 'TERMINATED_BY_ERROR';
	ok my $raw = $tool->result('');
	ok my $parsed = $tool->result('parsed');
	is $parsed->[0]{'coil'}, '1000';
	my @res = $tool->result('Bio::SeqFeatureI');
	if (scalar @res > 0) {
		ok 1;
	}
	else {
		skip 'No results - could not connect to HNN server?', 6;
	}
	
	ok my $meta = $tool->result('meta');
	ok my $seqobj = Bio::Seq->new(-primary_seq => $meta, display_id=>"a");
	ok $seqobj->add_SeqFeature($tool->result('Bio::SeqFeatureI'));
	
	test_skip(-tests => 2, -requires_module => 'Bio::Seq::Meta::Array');
	is $meta->named_submeta_text('HNN_helix',1,2), '0 111';
	is $meta->seq, 'MSADQRWRQDSQDSFGDSFDGDPPPPPPPPFGDSFGDGFSDRSRQDQRS';
}
