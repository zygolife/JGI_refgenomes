#!/usr/bin/env perl

eval 'exec /usr/bin/env perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
use strict;
use warnings;

=head1 NAME

bp_translate_seq - translates a sequence

=head1 SYNOPSIS

bp_translate_seq E<lt> cdna_cds.fa E<gt> protein.fa

=head1 DESCRIPTION 

The script will translate one fasta file (on stdin) to protein on stdout

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

  https://github.com/bioperl/bioperl-live/issues

=head1 AUTHOR

  Ewan Birney E<lt>birney@ebi.ac.ukE<gt>

=cut

use Bio::SeqIO;
use Getopt::Long;

my ($format) = 'fasta';
my $outfile;

GetOptions(
	   'format:s'  => \$format,
	   'o|out:s'   => \$outfile,
	   );

my $oformat = 'fasta';

# this implicity uses the <> file stream

my $seqin;
if ( @ARGV) {
	$seqin = Bio::SeqIO->new( -format => $format, -file => shift @ARGV);
} else {
	$seqin = Bio::SeqIO->new( -format => $format, -fh => \*STDIN);
}

my $seqout;
if ( $outfile ) {
	$seqout = Bio::SeqIO->new( -format => $oformat, -file => ">$outfile");
} else { $seqout = Bio::SeqIO->new( -format => $oformat); }



while( my $seq = $seqin->next_seq() ) {
	my $pseq = $seq->translate();
	$seqout->write_seq($pseq)
}

__END__
