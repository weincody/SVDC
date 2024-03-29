#!/usr/bin/perl -w
# $Id: seqretsplit.PLS 15088 2008-12-04 02:49:09Z bosborne $
#
=head1 NAME

seqretsplit - split a sequence (or stream) into a single file per sequence

=head1 SYNOPSIS

  seqretsplit file1 file2 .. 
  # or
  seqretsplit < file1 

=head1 DESCRIPTION 

The script will split all sequences from fasta file(s) (or stdin) to
individual files. The filename is the sequence ID (everything before
the 1st whitespace in a FASTA header).  Currently it doesn't check to
see that it isn't overwriting an existing file so IDs should be unique

This is inspired by EMBOSS seqretsplit tool.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
    http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

http://bugzilla.open-bio.org/

=head1 AUTHOR

Jason Stajich E<lt>jason_AT_bioperl_DOT_orgE<gt>

=cut

use strict;
use Bio::SeqIO;
my $in = Bio::SeqIO->new(-format => 'fasta',
			 -fh   => \*ARGV);
while( my $s = $in->next_seq ) {
    Bio::SeqIO->new(-format => 'fasta',
		    -file   => ">".$s->id.".fa")->write_seq($s);
}
    
