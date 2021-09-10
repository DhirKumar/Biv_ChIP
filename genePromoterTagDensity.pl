#!/usr/bin/perl

#### Adapted from geneTagDensity.pl on July 02, 2013 to enable generating tag density plots based on reads mapped
#### to positive or negative strands, i.e., GRO-Seq or RNA-Seq data.

### this script is generic enough to use for all types of data ChIP-Seq, RNA-Seq, etc

$numArgs = $#ARGV + 1;
if($numArgs != 16) 
{ 
  print "\nUSAGE: perl genePromoterTagDensity.pl -i <input-file> -b <input-file> -s <input-file> -n <1/0> -F <int> -r <1/0> -w <int> -o <output-file>\n\n"; 
  print "\t-i\tinput file fontaining list of genes in UCSC format\n";
  print "\t-b\tinput file containing ChIP-Seq reads/tags in BED format\n";
  print "\t-s\tinput file containing chromsome lengths in UCSC format(https://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/mm9.chrom.sizes)\n";
  print "\t-n\t1 for normalization by total number of reads in BED file, 0 otherwise\n";
  print "\t-F\tAverage ChIP-Seq fragment length for read shifting, e.g., 160, 200, 300 \n";
  print "\t-r\t1 (-1/0) if same (opposite, both resp.) strand reads need to be considered; e.g., Gro-Seq data\n";
  print "\t-w\tnumber of base pairs up- and down-stream of TSS, for promoter definition ; e.g., 500\n";
  print "\t-o\toutput file into which results will be stored\n";
  die "\n";
}

use Cwd;
use Getopt::Std;
my %Options;
$ok = getopts('i:b:s:n:o:F:r:w:', \%Options);
die "\n\nInvalid options on the command line\n\n" if (!$ok);

$windowSize = 100;
$windowPercent = 0.05;
my $max_bin_count=int(($upStreamDist+$downStreamDist)/$windowSize);###################

$geneFileName = $Options{i};
$bedFileName = $Options{b};
$outputFileName = $Options{o};
$normalize = $Options{n};
$readDirection = $Options{r};
$fragLength = $Options{F};
$shift = int($Options{F}/2+0.5);
$upStreamDist = $Options{w};
$downStreamDist = $Options{w};
$ChromLengthFileName= $Options{s};
#$t = localtime(); print "$t\tstart\n";
 
### HASHING CHROMOSOME LENGTHS
open(LEN, "$ChromLengthFileName") || die "can't open file $ChromLengthFileName\n";
%chromLength = ();
while(<LEN>)
{
  chomp;
  @words = split /\s+/, $_;
  $chromLength{$words[0]} = $words[1];
}
close LEN;

### HASHING BED COORDINATES
$numberTags = 0;
open(BED, "$bedFileName") || die "can't open file $bedFileName\n";
%bedStartsInChrom = ();
%posBedStartsInChrom = ();
%negBedStartsInChrom = ();
while(<BED>)
{ 
  $numberTags++;
  chomp; 
  @words = split /\t/, $_;      ### 0-chrom 1-bedStart 2-bedEnd 5-strand
  if ( !exists($bedStartsInChrom{$words[0]}) )
  {
    $coordPlus = $words[1]+$shift;
    $coordMinus = $words[2]-$shift;
    $bedStartsInChrom{$words[0]} = $coordPlus if ($words[5] eq "+");
    $bedStartsInChrom{$words[0]} = $coordMinus if ($words[5] eq "-");
    $posBedStartsInChrom{$words[0]} = $coordPlus if ($words[5] eq "+");
    $negBedStartsInChrom{$words[0]} = $coordMinus if ($words[5] eq "-");
  }
  else
  {
    $coordPlus = $words[1]+$shift;
    $coordMinus = $words[2]-$shift;
    $bedStartsInChrom{$words[0]} .= " ".$coordPlus if ($words[5] eq "+");
    $bedStartsInChrom{$words[0]} .= " ".$coordMinus if ($words[5] eq "-");
    $posBedStartsInChrom{$words[0]} .= " ".$coordPlus if ($words[5] eq "+");
    $negBedStartsInChrom{$words[0]} .= " ".$coordMinus if ($words[5] eq "-");
  }
}
close BED;


### HASHING GENE COORDINATES
open(GN, "$geneFileName") || die "can't open file $geneFileName\n";
%geneCoord = ();
%genesInChrom = ();
$head = <GN>;
while(<GN>)
{
  chomp;
  @words = split /\t/, $_;      ### 0-gene 2-chrom 3-txStart 4-txEnd 5-strand
				### 0-Gene 1-Chrom 2-Strand 3-TxStart 4-TxEnd 5-ExonCount 6-ExonStarts 7-ExonEnds
  $key = $words[0]."_".$words[1]."_".$words[3]."_".$words[4];
  $geneCoord{$key}= $words[1]."\t".$words[3]."\t".$words[4]."\t".$words[2];

  if (!exists($genesInChrom{$words[1]}) )
  {
    $genesInChrom{$words[1]} = $key;
  }
  else
  {
    $genesInChrom{$words[1]} .= " ".$key;
  }
}
close GN;

#$t = localtime(); print "$t\tdone hasing gene coordinates\n";


($upTagCount, $upSeqLength, $upTagDistr,$geneWiseTagCountref) = nearTSSTagDensities(-1);  # -1: upstream of TSS
###################################################################
#Added by dhirendra to get the gene wise tag count
###################################################################
my %geneWiseTagCount=%$geneWiseTagCountref;
my $gene_outfilename=$outputFileName;
$gene_outfilename.='_geneWise.txt';
open(GENE,">$gene_outfilename") or die ("Cannot create file $gene_outfilename");
foreach my $gene_d(keys(%geneWiseTagCount))
{
  my @gene_data=split(/_/,$gene_d);
  if(scalar(@{$geneWiseTagCount{$gene_d}})<$max_bin_count)
  {
    next;
  }
  if($gene_d=~/^NM|^NR/)
  {
    print GENE $gene_data[0],"_",$gene_data[1],"\t",join("\t",@gene_data[2..$#gene_data]),"\t";
  }
  else
  {
    print GENE join("\t",@gene_data),"\t";
  }
  #print join("\t",@{$geneWiseTagCount{$gene_d}}),"\n";
  foreach my $val(@{$geneWiseTagCount{$gene_d}})
  {
    my $val_RPM=($val*(10**9))/($numberTags*$windowSize);
    print GENE $val_RPM,"\t",
  }
  print GENE "\n";
}
close GENE;
###################################################################

#$t = localtime(); print "$t\tdone up and down stream\n";
open(OUT, ">$outputFileName") || die "can't open file $outputFileName\n";
$window = 0;
@upTagCountArr = @$upTagCount;
@upSeqLengthArr = @$upSeqLength;
@upTagDistrArr = @$upTagDistr;
@stdev;
for ($i = 0; $i <= $#upTagDistrArr; $i++)	### num windows
{
  #print "$upTagDistrArr[$i]\t";
  @numTags = split /\s/, $upTagDistrArr[$i];	### num genes (for eg) with window i
  $mean = 0;
  for ($k = 0; $k <= $#numTags; $k++)		### compute mean;
  {
     $mean += $numTags[$k]; 
  }
  $mean = $mean/($#numTags+1);
  #print "mean:$mean\t";
  $stdev[$i] = 0;
  for ($k = 0; $k <= $#numTags; $k++)		### compute stdev
  {
     $stdev[$i] += ( ($numTags[$k]-$mean)**2 );
  }
  $stdev[$i] = sqrt($stdev[$i]/($#numTags+1));
  #print "SD:$stdev[$i]\t";
  $stdev[$i] = $stdev[$i]/(sqrt($#numTags+1));	### confidence interval, CI = mean+/-(2*SE), where SE = stdev/sqrt(n)
  #print "CI:$stdev[$i]\n";
}
for ($i = 0; $i <= $#upTagCountArr; $i++)
{
  $tagDensity = 0;
  $tagDensity = $upTagCountArr[$i] / $upSeqLengthArr[$i] if ($upSeqLengthArr[$i] > 0);
  $stdErrorDensity = $stdev[$i] / $windowSize if ($upSeqLengthArr[$i] > 0);
  $tagDensity /= $numberTags if ($normalize == 1);
  $stdErrorDensity /= $numberTags if ($normalize == 1);
  print OUT "$window\t$upTagCountArr[$i]\t$upSeqLengthArr[$i]\t$tagDensity\t$stdErrorDensity\n";
  ###print OUT "$window\t$tagDensity\n";
  #print ">UP> $window\t$upTagCountArr[$i]\t$tagDensity\n";
  $window++;
}

close OUT;

####################### END OF MAIN ####################################

sub nearTSSTagDensities 	### call:nearTSSTagDensities(x); 
				### x=-1 upstream of bsite, 1 downstream of bsite 
{
  my($upDown)	= @_;		### upDown = -1 if upstream, 1 if downstream
  my(@tagStarts);		### array that stores the tag starts in a given chromosome
  my(@totalNumTags);		### one of the return arrays
  my(@totalSeqLength);		### another of the return arrays
  my(@tagDistribution);		### another of the return arrays
  my($nTags);
  my($start, $end, $windowNumIndex);
  my($numWindows);
  #$numWindows = $upStreamDist/$windowSize if ($upDown == -1);
  #$numWindows = $downStreamDist/$windowSize if ($upDown == 1); 
  $numWindows = ($upStreamDist + $downStreamDist)/$windowSize if ($upDown == -1); 

  ### INITIALIZATION
  for($i = 0; $i < $numWindows; $i++)
  {
    $totalNumTags[$i] = 0;
    $totalSeqLength[$i] = 0;
    $tagDistribution[$i] = "";
  }
my %genes_dist=();################
  foreach $chrom (sort keys%genesInChrom)	### for every chromosome with Genes 
  {
    #print "..$chrom\n";	
    @gene = split /\s/, $genesInChrom{$chrom};### all genes on this chromosome
    foreach my $key_d(@gene)
    {
      $genes_dist{$key_d}=[];#############
    }
    if ( exists($bedStartsInChrom{$chrom}) )	### if BED tags present in this chromosome
    {
      @tagStarts = split /\s/, $bedStartsInChrom{$chrom};
      @tagStarts = sort { $a <=> $b} @tagStarts;
      @posTagStarts = split /\s/, $posBedStartsInChrom{$chrom};
      @posTagStarts = sort { $a <=> $b} @posTagStarts;
      @negTagStarts = split /\s/, $negBedStartsInChrom{$chrom};
      @negTagStarts = sort { $a <=> $b} @negTagStarts;
	  
      for ($g = 0; $g <= $#gene; $g++) 		### for every gene on this chromosome
      {
	$key_d=$gene[$g];
        my (@terms) = split /\t/, $geneCoord{$gene[$g]}; ### 0-chrom 1-txStart 2-txEnd 3-strand
	if ($terms[3] eq "+")
	{
	  $start = $terms[1] - $upStreamDist if ($upDown == -1);
	  $start = $terms[1] if ($upDown == 1);
	  $end = $start + $windowSize - 1;	### -1:  end-start+1 = windowSize
	  $windowNumIndex = 0;
	  if ( ($start >= 0 && $upDown == -1) || ( ($terms[1]+$downStreamDist) <= $chromLength{$chrom} && $upDown == 1) )
	  {
	    if ($start > $tagStarts[$#tagStarts] || $end < $tagStarts[0] || ( ($start > $posTagStarts[$#posTagStarts]) && $readDirection == 1) || ($end < $posTagStarts[0] && $readDirection == 1) || ( ($start > $negTagStarts[$#negTagStarts]) && $readDirection == 1) || ($end < $negTagStarts[0] && $readDirection == -1))
	    {
	      $totalNumTags[$windowNumIndex] += 0; 
              $genes_dist{$key_d}->[$windowNumIndex]+=0;
	      $endIndex = -1 if ($end < $tagStarts[0] && $readDirection == 0); 
	      $endIndex = -1 if ($end < $posTagStarts[0] && $readDirection == 1); 
	      $endIndex = -1 if ($end < $negTagStarts[0] && $readDirection == -1); 
	      $tagDistribution[$windowNumIndex] .= "0" if ($tagDistribution[$windowNumIndex] eq "");
	      $tagDistribution[$windowNumIndex] .= " 0" if ($tagDistribution[$windowNumIndex] ne "");
	    }
	    else
	    {
              if ($readDirection == 0)
	      {
  	        $startIndex = binarySearch(1, $start, \@tagStarts);	### 1 cos' window start
  	        $endIndex = binarySearch(0, $end, \@tagStarts);		### 0 cos' window end
	      }
	      elsif ($readDirection == 1)
	      {
  	        $startIndex = binarySearch(1, $start, \@posTagStarts);	### 1 cos' window start
  	        $endIndex = binarySearch(0, $end, \@posTagStarts);	### 0 cos' window end
	      }
 	      else  ## $readDirection == -1
	      {
  	        $startIndex = binarySearch(1, $start, \@negTagStarts);	### 1 cos' window start
  	        $endIndex = binarySearch(0, $end, \@negTagStarts);	### 0 cos' window end
	      }

	      $totalNumTags[$windowNumIndex] += ($endIndex-$startIndex+1);
              $genes_dist{$key_d}->[$windowNumIndex]+=($endIndex-$startIndex+1);##################
	      $nTags = $endIndex-$startIndex+1;
	      $tagDistribution[$windowNumIndex] .= $nTags  if ($tagDistribution[$windowNumIndex] eq "");
	      $tagDistribution[$windowNumIndex] .= " ".$nTags if ($tagDistribution[$windowNumIndex] ne "");
	    }
	    $totalSeqLength[$windowNumIndex] += ($end - $start + 1);
	    for ($k = 1; $k <= $numWindows-1; $k++)
	    {
	      $start += $windowSize;
	      $end += $windowSize;
	      $windowNumIndex++;
	      if ($start > $tagStarts[$#tagStarts] || $end < $tagStarts[0] || ($start > $posTagStarts[$#posTagStarts] && $readDirection == 1) || ($end < $posTagStarts[0] && $readDirection == 1) || ($start > $negTagStarts[$#negTagStarts] && $readDirection == 1) || ($end < $negTagStarts[0] && $readDirection == -1) )
	      {
	        $totalNumTags[$windowNumIndex] += 0; 
                $genes_dist{$key_d}->[$windowNumIndex]+=0;
	        $endIndex = -1 if ($end < $tagStarts[0] && $readDirection == 0);  
	        $endIndex = -1 if ($end < $posTagStarts[0] && $readDirection == 1); 
	        $endIndex = -1 if ($end < $negTagStarts[0] && $readDirection == -1); 
	        $tagDistribution[$windowNumIndex] .= "0" if ($tagDistribution[$windowNumIndex] eq "");
	        $tagDistribution[$windowNumIndex] .= " 0" if ($tagDistribution[$windowNumIndex] ne "");
	      }
	      else
	      {
	        if ($readDirection == 0)
		{
  	          $startIndex = $endIndex + 1;	#binarySearch(1, $start, \@tagStarts);	### 1 cos' window start
  	          $endIndex = binarySearch(0, $end, \@tagStarts);		### 0 cos' window end
		}
		elsif ($readDirection == 1)
		{
  	          $startIndex = $endIndex + 1;	#binarySearch(1, $start, \@tagStarts);	### 1 cos' window start
  	          $endIndex = binarySearch(0, $end, \@posTagStarts);		### 0 cos' window end
		}
		else ## $readDirection == -1)
		{
  	          $startIndex = $endIndex + 1;	#binarySearch(1, $start, \@tagStarts);	### 1 cos' window start
  	          $endIndex = binarySearch(0, $end, \@negTagStarts);		### 0 cos' window end
		}
	        $totalNumTags[$windowNumIndex] += ($endIndex-$startIndex+1);
			$genes_dist{$key_d}->[$windowNumIndex]+=($endIndex-$startIndex+1);##################
	        $nTags = $endIndex-$startIndex+1;
	        $tagDistribution[$windowNumIndex] .= $nTags  if ($tagDistribution[$windowNumIndex] eq "");
	        $tagDistribution[$windowNumIndex] .= " ".$nTags if ($tagDistribution[$windowNumIndex] ne "");
	      }
	      $totalSeqLength[$windowNumIndex] += ($end - $start + 1);
	    }
	  }
	}
	elsif ($terms[3] eq "-")
	{
	  $end = $terms[2] + $upStreamDist if ($upDown == -1);
	  $end = $terms[2] if ($upDown == 1);
	  $start = $end - $windowSize + 1;	
	  $windowNumIndex = 0;
	  if ( ($end <= $chromLength{$chrom} && $upDown == -1) || ( ($terms[2]-$downStreamDist) >= 0 && $upDown == 1) )
	  {
	    if ($start > $tagStarts[$#tagStarts] || $end < $tagStarts[0] || ($start > $negTagStarts[$#negTagStarts] && $readDirection == 1) || ($end < $negTagStarts[0] && $readDirection == 1) || ($start > $posTagStarts[$#posTagStarts] && $readDirection == 1) || ($end < $posTagStarts[0] && $readDirection == -1) )
	    {
	      $totalNumTags[$windowNumIndex] += 0 ;
	      $genes_dist{$key_d}->[$windowNumIndex]+=0;
	      $startIndex = $#tagStarts+1 if ($start > $tagStarts[$#tagStarts] && $readDirection == 0);
	      $startIndex = $#negTagStarts+1 if ($start > $negTagStarts[$#negTagStarts] && $readDirection == 1);
	      $startIndex = $#posTagStarts+1 if ($start > $posTagStarts[$#posTagStarts] && $readDirection == -1);
	      $tagDistribution[$windowNumIndex] .= "0" if ($tagDistribution[$windowNumIndex] eq "");
	      $tagDistribution[$windowNumIndex] .= " 0" if ($tagDistribution[$windowNumIndex] ne "");
	    }
	    else
	    {
	      if ($readDirection == 0)
	      {
  	        $startIndex = binarySearch(1, $start, \@tagStarts);	### 1 cos' window start
  	        $endIndex = binarySearch(0, $end, \@tagStarts);		### 0 cos' window end
	      }
	      elsif ($readDirection == 1)
	      {
  	        $startIndex = binarySearch(1, $start, \@negTagStarts);	### 1 cos' window start
  	        $endIndex = binarySearch(0, $end, \@negTagStarts);		### 0 cos' window end
	      }
	      else ##$readDirection == -1
	      {
  	        $startIndex = binarySearch(1, $start, \@posTagStarts);	### 1 cos' window start
  	        $endIndex = binarySearch(0, $end, \@posTagStarts);		### 0 cos' window end
	      }
	      $totalNumTags[$windowNumIndex] += ($endIndex-$startIndex+1);
		  $genes_dist{$key_d}->[$windowNumIndex]+=($endIndex-$startIndex+1);##################
	      $nTags = $endIndex-$startIndex+1;
	      $tagDistribution[$windowNumIndex] .= $nTags  if ($tagDistribution[$windowNumIndex] eq "");
	      $tagDistribution[$windowNumIndex] .= " ".$nTags if ($tagDistribution[$windowNumIndex] ne "");
	    }
	    $totalSeqLength[$windowNumIndex] += ($end - $start + 1);
	    for ($k = 1; $k <= $numWindows-1; $k++)
	    {
	      $start -= $windowSize;
	      $end -= $windowSize;
	      $windowNumIndex++;
	      if ($start > $tagStarts[$#tagStarts] || $end < $tagStarts[0] || ($start > $negTagStarts[$#negTagStarts] && $readDirection == 1) || ($end < $negTagStarts[0] && $readDirection == 1) || ($start > $posTagStarts[$#posTagStarts] && $readDirection == 1) || ($end < $posTagStarts[0] && $readDirection == -1) )
	      {
	        $totalNumTags[$windowNumIndex] += 0; 
                $genes_dist{$key_d}->[$windowNumIndex]+=0;
	        $startIndex = $#tagStarts+1 if ($start > $tagStarts[$#tagStarts] && $readDirection == 0);
	        $startIndex = $#negTagStarts+1 if ($start > $negTagStarts[$#negTagStarts] && $readDirection == 1);
	        $startIndex = $#posTagStarts+1 if ($start > $posTagStarts[$#posTagStarts] && $readDirection == -11);
	        $tagDistribution[$windowNumIndex] .= "0" if ($tagDistribution[$windowNumIndex] eq "");
	        $tagDistribution[$windowNumIndex] .= " 0" if ($tagDistribution[$windowNumIndex] ne "");
	      }
	      else
	      {
	        if ($readDirection == 0)
		{
	          $endIndex = $startIndex - 1;
  	          $startIndex = binarySearch(1, $start, \@tagStarts);	### 1 cos' window start
		}
		elsif ($readDirection == 1)
		{
	          $endIndex = $startIndex - 1;
  	          $startIndex = binarySearch(1, $start, \@negTagStarts);### 1 cos' window start
		}
		else ##$readDirection == -1
		{
	          $endIndex = $startIndex - 1;
  	          $startIndex = binarySearch(1, $start, \@posTagStarts);### 1 cos' window start
		}

	        $totalNumTags[$windowNumIndex] += ($endIndex-$startIndex+1);
			$genes_dist{$key_d}->[$windowNumIndex]+=($endIndex-$startIndex+1);##################
	        $nTags = $endIndex-$startIndex+1;
	        $tagDistribution[$windowNumIndex] .= $nTags  if ($tagDistribution[$windowNumIndex] eq "");
	        $tagDistribution[$windowNumIndex] .= " ".$nTags if ($tagDistribution[$windowNumIndex] ne "");
	      }
	      $totalSeqLength[$windowNumIndex] += ($end - $start + 1);
	    }
	  }
	}
	
	
      }
    }
    undef @tagStarts;
    undef @posTagStarts;
    undef @negTagStarts;
  }
  return (\@totalNumTags, \@totalSeqLength, \@tagDistribution,\%genes_dist);
}


sub nearTESTagDensities 	### call:nearTESTagDensities(x); 
				### x=-1 upstream of bsite, 1 downstream of bsite 
{
  my($upDown)	= @_;		### upDown = -1 if upstream, 1 if downstream
  my(@tagStarts);		### array that stores the tag starts in a given chromosome
  my(@totalNumTags);		### one of the return arrays
  my(@totalSeqLength);		### another of the return arrays
  my($start, $end, $windowNumIndex);
  my($numWindows);
  #$numWindows = $upStreamDist/$windowSize if ($upDown == -1);
  #$numWindows = $downStreamDist/$windowSize if ($upDown == 1); 
  $numWindows = ($upStreamDist + $downStreamDist)/$windowSize if ($upDown == -1); 

  ### INITIALIZATION
  for($i = 0; $i < $numWindows; $i++)
  {
    $totalNumTags[$i] = 0;
    $totalSeqLength[$i] = 0;
  }

  foreach $chrom (sort keys%bedStartsInChrom)	### for every chromosome with BED tags
  {
    #print "..$chrom\n";
    @tagStarts = split /\s/, $bedStartsInChrom{$chrom};
    @tagStarts = sort { $a <=> $b} @tagStarts;
    if ( exists($genesInChrom{$chrom}) )	### if genes present in this chromosome
    {
      @gene = split /\s/, $genesInChrom{$chrom};### all genes on this chromosome
      for ($g = 0; $g <= $#gene; $g++) 		### for every gene on this chromosome
      {
        my (@terms) = split /\t/, $geneCoord{$gene[$g]}; ### 0-chrom 1-txStart 2-txEnd 3-strand
	if ($terms[3] eq "+")
	{
	  $start = $terms[2] - $upStreamDist if ($upDown == -1);
	  $start = $terms[2] if ($upDown == 1);
	  $end = $start + $windowSize - 1;	### -1:  end-start+1 = windowSize
	  $windowNumIndex = 0;
	  if ( ($start >= 0 && $upDown == -1) || ( ($terms[2]+$downStreamDist) <= $chromLength{$chrom} && $upDown == 1) )
	  {
	    ###$totalNumTags[$windowNumIndex] += countTagsInWindow($start, $end, @tagStarts);
	    ###$totalSeqLength[$windowNumIndex] += ($end - $start + 1);
	    #-----------------------
	    if ($start > $tagStarts[$#tagStarts] || $end < $tagStarts[0])
	    {
	      $totalNumTags[$windowNumIndex] += 0;
	      $endIndex = -1 if ($end < $tagStarts[0]);
	    }
	    else
	    {
  	      $startIndex = binarySearch(1, $start, \@tagStarts);	### 1 cos' window start
  	      $endIndex = binarySearch(0, $end, \@tagStarts);		### 0 cos' window end
	      $totalNumTags[$windowNumIndex] += ($endIndex-$startIndex+1);
	    }
	    $totalSeqLength[$windowNumIndex] += ($end - $start + 1);
	    #-----------------------
	    for ($k = 1; $k <= $numWindows-1; $k++)
	    {
	      $start += $windowSize;
	      $end += $windowSize;
	      $windowNumIndex++;
	      ###$totalNumTags[$windowNumIndex] += countTagsInWindow($start, $end, @tagStarts);
	      ###$totalSeqLength[$windowNumIndex] += ($end - $start + 1);
	      #-----------------------
	      if ($start > $tagStarts[$#tagStarts] || $end < $tagStarts[0])
	      {
	        $totalNumTags[$windowNumIndex] += 0; 
	        $endIndex = -1 if ($end < $tagStarts[0]);
	      }
	      else
	      {
  	        $startIndex = $endIndex + 1;	#binarySearch(1, $start, \@tagStarts);	### 1 cos' window start
  	        $endIndex = binarySearch(0, $end, \@tagStarts);		### 0 cos' window end
	        $totalNumTags[$windowNumIndex] += ($endIndex-$startIndex+1);
	      }
	      $totalSeqLength[$windowNumIndex] += ($end - $start + 1);
	      #-----------------------
	    }
	  }
	}
	elsif ($terms[3] eq "-")
	{
	  $end = $terms[1] + $upStreamDist if ($upDown == -1);
	  $end = $terms[1] if ($upDown == 1);
	  $start = $end - $windowSize + 1;	
	  $windowNumIndex = 0;
	  if ( ($end <= $chromLength{$chrom} && $upDown == -1) || ( ($terms[1]-$downStreamDist) >= 0 && $upDown == 1) )
	  {
	    ###$totalNumTags[$windowNumIndex] += countTagsInWindow($start, $end, @tagStarts);
	    ###$totalSeqLength[$windowNumIndex] += ($end - $start + 1);
	    #-----------------------
	    if ($start > $tagStarts[$#tagStarts] || $end < $tagStarts[0])
	    {
	      $totalNumTags[$windowNumIndex] += 0 ;
	      $startIndex = $#tagStarts+1 if ($start > $tagStarts[$#tagStarts]);
	    }
	    else
	    {
  	      $startIndex = binarySearch(1, $start, \@tagStarts);	### 1 cos' window start
  	      $endIndex = binarySearch(0, $end, \@tagStarts);		### 0 cos' window end
	      $totalNumTags[$windowNumIndex] += ($endIndex-$startIndex+1);
	    }
	    $totalSeqLength[$windowNumIndex] += ($end - $start + 1);
	    #-----------------------
	    for ($k = 1; $k <= $numWindows-1; $k++)
	    {
	      $start -= $windowSize;
	      $end -= $windowSize;
	      $windowNumIndex++;
	      ###$totalNumTags[$windowNumIndex] += countTagsInWindow($start, $end, @tagStarts);
	      ###$totalSeqLength[$windowNumIndex] += ($end - $start + 1);
	      #-----------------------
	      if ($start > $tagStarts[$#tagStarts] || $end < $tagStarts[0])
	      {
	        $totalNumTags[$windowNumIndex] += 0; 
	        $startIndex = $#tagStarts+1 if ($start > $tagStarts[$#tagStarts]);
	      }
	      else
	      {
	        $endIndex = $startIndex - 1;
  	        $startIndex = binarySearch(1, $start, \@tagStarts);	### 1 cos' window start
	        $totalNumTags[$windowNumIndex] += ($endIndex-$startIndex+1);
	      }
	      $totalSeqLength[$windowNumIndex] += ($end - $start + 1);
	      #-----------------------
	    }
	  }
	}
      }
    }
  }
  return (\@totalNumTags, \@totalSeqLength);
}

sub intragenicTagDensities 
{
  my(@tagStarts);		### array that stores the tag starts in a given chromosome
  my(@totalNumTags);		### one of the return arrays
  my(@totalSeqLength);		### another of the return arrays
  my(@tagDistribution);
  my($nTags);
  my($start, $end, $windowNumIndex);
  my($numWindows);
  my($startIndex, $endIndex);
  $numWindows = 1.0 / (1.0*$windowPercent);

  ### INITIALIZATION
  for($i = 0; $i < $numWindows; $i++)
  {
    $totalNumTags[$i] = 0;
    $totalSeqLength[$i] = 0;
    $tagDistribution[$i] = "";
  }

  foreach $chrom (sort keys%bedStartsInChrom)	### for every chromosome with BED tags
  {
    #print "..$chrom\n";
    @tagStarts = split /\s/, $bedStartsInChrom{$chrom};
    @tagStarts = sort { $a <=> $b} @tagStarts;
    if ( exists($genesInChrom{$chrom}) )	### if genes present in this chromosome
    {
      @gene = split /\s/, $genesInChrom{$chrom};### all genes on this chromosome
      for ($g = 0; $g <= $#gene; $g++) 		### for every gene on this chromosome
      {
        my (@terms) = split /\t/, $geneCoord{$gene[$g]}; ### 0-chrom 1-txStart 2-txEnd 3-strand
	#print "$terms[1]\t$terms[2]\t$terms[3]\t$numBasesInWindow";
        $geneLength = abs($terms[2]-$terms[1]+1);
        $numBasesInWindow = $geneLength * (1.0*$windowPercent);
	if ($terms[3] eq "+")
	{
	  $start = $terms[1];
	  $end = $start + int($numBasesInWindow + 0.5) - 1;	### +0.5: rounding
	  $windowNumIndex = 0;
	  if ( $start >= 0 )
	  {
	    if ($start > $tagStarts[$#tagStarts] || $end < $tagStarts[0])
	    {
	      $totalNumTags[$windowNumIndex] += 0;
	      $endIndex = -1 if ($end < $tagStarts[0]);
	      $tagDistribution[$windowNumIndex] .= "0" if ($tagDistribution[$windowNumIndex] eq "");
	      $tagDistribution[$windowNumIndex] .= " 0" if ($tagDistribution[$windowNumIndex] ne "");
	    }
	    else
	    {
  	      $startIndex = binarySearch(1, $start, \@tagStarts);	### 1 cos' window start
  	      $endIndex = binarySearch(0, $end, \@tagStarts);		### 0 cos' window end
	      $totalNumTags[$windowNumIndex] += ($endIndex-$startIndex+1);
	      $nTags = $endIndex-$startIndex+1;
	      $tagDistribution[$windowNumIndex] .= $nTags  if ($tagDistribution[$windowNumIndex] eq "");
	      $tagDistribution[$windowNumIndex] .= " ".$nTags if ($tagDistribution[$windowNumIndex] ne "");
	    }
	    $totalSeqLength[$windowNumIndex] += ($end - $start + 1);
	    #$dummy = $end - $start + 1; print "\t0 - $start\t$end\t$dummy\n";
	    for ($k = 1; $k <= $numWindows-1; $k++)
	    {
	      $start = $end + 1;
              #$tt = int( (($k+1)*$numBasesInWindow) + 0.5 ); 
              $end = $terms[1] + int( (($k+1)*$numBasesInWindow) + 0.5 );
	      $windowNumIndex++;
	      if ($start > $tagStarts[$#tagStarts] || $end < $tagStarts[0])
	      {
	        $totalNumTags[$windowNumIndex] += 0; 
	        $endIndex = -1 if ($end < $tagStarts[0]);
	        $tagDistribution[$windowNumIndex] .= "0" if ($tagDistribution[$windowNumIndex] eq "");
	        $tagDistribution[$windowNumIndex] .= " 0" if ($tagDistribution[$windowNumIndex] ne "");
	      }
	      else
	      {
  	        $startIndex = $endIndex + 1;	#binarySearch(1, $start, \@tagStarts);	### 1 cos' window start
  	        $endIndex = binarySearch(0, $end, \@tagStarts);		### 0 cos' window end
	        $totalNumTags[$windowNumIndex] += ($endIndex-$startIndex+1);
	        $nTags = $endIndex-$startIndex+1;
	        $tagDistribution[$windowNumIndex] .= $nTags  if ($tagDistribution[$windowNumIndex] eq "");
	        $tagDistribution[$windowNumIndex] .= " ".$nTags if ($tagDistribution[$windowNumIndex] ne "");
	      }
	      $totalSeqLength[$windowNumIndex] += ($end - $start + 1);
	      #$dummy = $end - $start + 1; print "\t$k - $start\t$end\t$dummy\n";
	    }
	  }
	}
	elsif ($terms[3] eq "-")
	{
	  $end = $terms[2];
	  $start = $end - int($numBasesInWindow + 0.5) + 1; 	### +0.5: rounding
	  $windowNumIndex = 0;
	  if ( $end <= $chromLength{$chrom} )
	  {
	    if ($start > $tagStarts[$#tagStarts] || $end < $tagStarts[0])
	    {
	      $totalNumTags[$windowNumIndex] += 0 ;
	      $startIndex = $#tagStarts+1 if ($start > $tagStarts[$#tagStarts]);
	      $tagDistribution[$windowNumIndex] .= "0" if ($tagDistribution[$windowNumIndex] eq "");
	      $tagDistribution[$windowNumIndex] .= " 0" if ($tagDistribution[$windowNumIndex] ne "");
	    }
	    else
	    {
  	      $startIndex = binarySearch(1, $start, \@tagStarts);	### 1 cos' window start
  	      $endIndex = binarySearch(0, $end, \@tagStarts);		### 0 cos' window end
	      $totalNumTags[$windowNumIndex] += ($endIndex-$startIndex+1);
	      $nTags = $endIndex-$startIndex+1;
	      $tagDistribution[$windowNumIndex] .= $nTags  if ($tagDistribution[$windowNumIndex] eq "");
	      $tagDistribution[$windowNumIndex] .= " ".$nTags if ($tagDistribution[$windowNumIndex] ne "");
	    }
	    $totalSeqLength[$windowNumIndex] += ($end - $start + 1);
	    #$dummy = $end - $start + 1; print "\t0 - $start\t$end\t$dummy\n";
	    for ($k = 1; $k <= $numWindows-1; $k++)
	    {
	      $end = $start - 1;
              $start = ($terms[2]+1) - int( (($k+1)*$numBasesInWindow) + 0.5 );
	      $windowNumIndex++;
	      if ($start > $tagStarts[$#tagStarts] || $end < $tagStarts[0])
	      {
	        $totalNumTags[$windowNumIndex] += 0; 
	        $startIndex = $#tagStarts+1 if ($start > $tagStarts[$#tagStarts]);
	        $tagDistribution[$windowNumIndex] .= "0" if ($tagDistribution[$windowNumIndex] eq "");
	        $tagDistribution[$windowNumIndex] .= " 0" if ($tagDistribution[$windowNumIndex] ne "");
	      }
	      else
	      {
	        $endIndex = $startIndex - 1;
  	        $startIndex = binarySearch(1, $start, \@tagStarts);	### 1 cos' window start
	        $totalNumTags[$windowNumIndex] += ($endIndex-$startIndex+1);
	        $nTags = $endIndex-$startIndex+1;
	        $tagDistribution[$windowNumIndex] .= $nTags  if ($tagDistribution[$windowNumIndex] eq "");
	        $tagDistribution[$windowNumIndex] .= " ".$nTags if ($tagDistribution[$windowNumIndex] ne "");
	      }
	      $totalSeqLength[$windowNumIndex] += ($end - $start + 1);
	      #$dummy = $end - $start + 1; print "\t$k - $start\t$end\t$dummy\n";
	    }
	  }
	}
      }
    }
  }
  return (\@totalNumTags, \@totalSeqLength,\@tagDistribution);
}



sub countTagsInWindow	### call: $numTags = countTagsInWindow($start, $end, @tagStarts)
{
  my($start, $end, @tagStarts);
  ($start, $end, @tagStarts) = @_;

  ### if window start falls to right of tagStarts or if window end falls to left of tagStarts
  return 0 if ($start > $tagStarts[$#tagStarts] || $end < $tagStarts[0]);

  my ($startIndex) = binarySearch(1, $start, \@tagStarts);	# 1 cos' window start
  my ($endIndex) = binarySearch(0, $end, \@tagStarts);		# 0 cos' window end
  return $endIndex-$startIndex+1;
}

sub binarySearch		### call: $returnValue = binarySearch($begin, $query, \@starts);
{
  my ($begin, $query, $starts);	### begin = 1 if query is a window start, 0 if window end 
  ($begin, $query, $starts) = @_;
  #my(@starts) = @$startsPtr;

  ### Handling cases when query is out of bounds
  return 0 if ($query < $starts->[0] && $begin == 1);
  return -1 if ($query < $starts->[0] && $begin == 0);
  return $#$starts if ($query > $starts->[$#$starts] && $begin == 0);
  return $#$starts+1 if ($query > $starts->[$#$starts] && $begin == 1);

  my ($left) = 0;
  my ($right) = $#$starts;
  my ($center);
  my ($prevCenter) = -1;
  while(1)
  {
    $center = int (($left + $right)/2 );
    if ($query == $starts->[$center])
    {
      if ($begin == 1)
      {
        while($center > 0 && $starts->[$center] == $starts->[$center-1]) 
	  { $center =  $center - 1; }
      }
      else
      {
        while($center < $#$starts && $starts->[$center] == $starts->[$center+1]) 
	  { $center =  $center + 1; }
      }
      return $center;
    }
    if ($center == $prevCenter)
    {
      return $right if ($begin == 1);
      return $right-1 if ($begin == 0);
    }
    $right = $center if ($query < $starts->[$center]);
    $left = $center if ($query > $starts->[$center]);
    $prevCenter = $center;
  }
}

