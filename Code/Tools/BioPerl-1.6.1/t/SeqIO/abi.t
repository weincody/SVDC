# -*-Perl-*- Test Harness script for Bioperl
# $Id: abi.t 15300 2009-01-06 01:00:43Z cjfields $

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 7,
			   -requires_module => 'Bio::SeqIO::staden::read 1.006');
	
	use_ok('Bio::SeqIO');
}

my $verbose = test_debug();

my $io = Bio::SeqIO->new(-format => 'abi',
								 -verbose => $verbose,
								 -file => test_input_file('test.abi',),
								 -get_trace_data => 1);
my $seq = $io->next_seq;
isa_ok($seq, 'Bio::PrimarySeqI');
is($seq->seq, "GCNTATGACGTGGATTNCGAATTCTNNNNNCGGTAGNNGAAAATCCCCGGNCAAGNTTNNCCCTGCAAANGGAANAANNTGGCCGAGCGCTACGGGCTGATCTGGGTGTGCCTGTTTCCCCCGGCCGGGGGGAGNGATGCAGGACATCCAAGTATCCCGCCNATGGNGGGCTGAGGACGAGGACGGCTTCCATCAGATCAGTGTGCCCGGNCTTCGACATCGGCGGCAGCGCCGCGCGCCAACTGGAAGGCTTCATCGACGTGNAGCATTTTGNCTTCNTGCGCACCGCTACCTTCACCCANCCGGACAAGCGCNAANTGCNGNCCTACACCACCACNGAAACACCGACCGGNTTNAATGCCGATTACCTGAGNNGCGTGGCAAATTATTCGGNGGACNTGCCGCTGNCGGACGTGGACCCGAACTTCCAATGGCTGCGTCATTNCTAGGTGAATCTGCCTTTCACCGCCACGCTCACCATCCACTTCCCGGTGCCGGGCAAGCGGTTGGTGATNATGAATGCCGCCAGACCGGTGTCCAAGCACACCANCCGCCTGNTGGTGCCGATCGNCCGCTAATTTCGACACCCATCTGCCNGNGGGAAGACGTACATGNGTTCAACCTTGCACNTNGTTCNAAAAAAACCNTGCCATGGTGGNAANCGCAAGCGGNCCGGAAATATCNGCCGGNTTGACCCGCNTGNTTGGAAAGTGCATATTCCCCNCCGATNCNCAATTTCGAT");
# trace data points, only added if get_trace_data is invoked
is($seq->get_trace_graph( -trace => 'a' ), 8793); 
is($seq->get_trace_graph( -trace => 't' ), 8793); 
is($seq->get_trace_graph( -trace => 'g' ), 8793); 
is($seq->get_trace_graph( -trace => 'c' ), 8793); 