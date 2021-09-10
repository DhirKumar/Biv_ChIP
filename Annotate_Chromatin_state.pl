#H3K4me3_enriched_genes_0hr.txt H3K27me3_enriched_genes_0hr.txt

use strict;
use warnings;
use Cwd;
use Getopt::Std;

my $numArgs = $#ARGV + 1;
if($numArgs != 10)
{
  print "\nUSAGE: Annotate_Chromatin_state.pl -A <file> -B <file> -P <file> -o <file>\n";
  print "\t-A\tFile containing gene/transcript-IDS enriched for H3K4me3 \n";
  print "\t-B\tFile containing gene/transcript-IDS enriched for H3K27me3\n";
  print "\t-P\tFile containing promoter information with transcript IDs in the first column of the file\n";
  print "\t-o\t*name of the output file\n";
  die "\n";
}

my %Options;
my $ok = getopts('A:B:P:o:', \%Options);
die "\n\nInvalid options on the command line\n\n" if (!$ok);

my $K4_file = $Options{A};
my $K27_file = $Options{B};
my $promoter_file = $Options{P};
#my $lab=$Options{l};
my $outputFileName = $Options{o};


my %chromatin_state=();
open(IN,$promoter_file) || die("Cannot open $promoter_file");
while(my $line=<IN>)
{
  my @line_data=split(/\s+/,$line);
  $chromatin_state{$line_data[0]}=[];
}
my $i=0;
#while(my $line=<DATA>)
#{
#my @line_data=split(/\s+/,$line);
#my $K4_file=$line_data[0];
#my $K27_file=$line_data[1];
open(K4,$K4_file) || die("Cannot open $K4_file");
open(K27,$K27_file) || die("Cannot open $K27_file");
my %K4_enriched=();
my %K27_enriched=();
while(my $K4_line=<K4>)
{
  my @line_data=split(/\s+/,$K4_line);
  $K4_enriched{$line_data[0]}=1;
}
close K4;
while(my $K27_line=<K27>)
{
  my @line_data=split(/\s+/,$K27_line);
  $K27_enriched{$line_data[0]}=1;
}
close K27;

foreach my $key(keys(%chromatin_state))
{
  if((exists($K4_enriched{$key})) && (exists($K27_enriched{$key})))
  {
    $chromatin_state{$key}[$i]='Bivalent';
  }
  elsif(exists($K4_enriched{$key}))
  {
    $chromatin_state{$key}[$i]='K4';
  }
  elsif(exists($K27_enriched{$key}))
  {
    $chromatin_state{$key}[$i]='K27';
  }
  else
  {
    $chromatin_state{$key}[$i]='none';
  }
}

open(OUT,">$outputFileName")|| die("Cannot create output file name");
print OUT "GENE\tChromatin-State\n",;
foreach my $key(keys(%chromatin_state))
{
  print OUT $key,"\t",join("\t",@{$chromatin_state{$key}}),"\n";
}
close OUT;
