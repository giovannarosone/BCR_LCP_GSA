#!/usr/bin/perl


use strict;
use warnings;


die("usage: ./fromFastQ2.pl fastqFile") unless ($#ARGV==0);

open(FILE,"$ARGV[0]");

open(FASTA,">$ARGV[0].fasta");
#open(SEQS,">$ARGV[0].seqs");
#open(MANZINI,">$ARGV[0].cat");
open(QUALSCORE,">$ARGV[0].fasta.qs");


my $cur_seq_header = "";
my $cur_seq_seq = "";
my $cur_qual_header = "";
my $cur_qual_seq = "";


while( <FILE> ) 
{
    $cur_seq_header = $_;
    $cur_seq_seq = <FILE>;
    $cur_qual_header = <FILE>;
    $cur_qual_seq = <FILE>;

    print FASTA ">" . $cur_seq_header;
    print FASTA $cur_seq_seq;
    
    print QUALSCORE "{" . $cur_seq_header;
    print QUALSCORE $cur_qual_seq;
	
    #print SEQS $cur_seq_seq;

    chomp($cur_seq_seq);
    #print MANZINI $cur_seq_seq . "0";

}



close(FILE);
