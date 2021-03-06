# -*-Perl-*- Test Harness script for Bioperl
# $Id: phylip.t 14971 2008-10-28 16:08:52Z cjfields $

use strict;

BEGIN {
	use lib '.';
    use Bio::Root::Test;
    
    test_begin(-tests => 11);
	
	use_ok('Bio::AlignIO::phylip');
}

my $DEBUG = test_debug();

my ($str,$aln,$strout,$status);

# PHYLIP interleaved
$str = Bio::AlignIO->new(
    '-file' => test_input_file("testaln.phylip"),
    '-format' => 'phylip');
isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'Homo_sapie/1-45';

$strout = Bio::AlignIO->new(
    '-file'  => ">".test_output_file(),
    '-format' => 'phylip');
$status = $strout->write_aln($aln);
is $status, 1, "phylip output test";

# PHYLIP sequential/non-interleaved
$strout = Bio::AlignIO->new('-file'  => test_input_file('noninterleaved.phy'),
			    '-format' => 'phylip');
$aln = $strout->next_aln($aln);
isa_ok($aln,'Bio::Align::AlignI');
is($aln->get_seq_by_pos(2)->seq(), 'CCTCAGATCACTCTTTGGCAACGACCCCTCGTCACAATAA'.
   'AGGTAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGACATGAATT'.
   'TGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGGTTTATCAAAGTAAGACAGTATGATCAGA'.
   'TACCCATAGAGATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCCACACCTGTCAATATAATTG'.
   'GAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT' );

# PHYLIP interleaved with long Ids
$str = Bio::AlignIO->new(
    '-file' => test_input_file("protpars_longid.phy"),
    '-format' => 'phylip',
    'longid' => 1);

isa_ok($str,'Bio::AlignIO');
$aln = $str->next_aln();
isa_ok($aln,'Bio::Align::AlignI');
is $aln->get_seq_by_pos(1)->get_nse, 'S I N F R U  P 0 0 1 /1-84';
is $aln->get_seq_by_pos(2)->get_nse, 'SINFRUP002/1-84';

