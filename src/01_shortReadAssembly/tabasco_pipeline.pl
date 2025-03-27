#!perl -w
use strict;
use warnings FATAL => 'all';
# Quantifies assembly completeness via alignment to 6000 Andropogoneae core genes
# Baoxing Song, modified by Charlie Hale (chale295@gmail.com)
# USAGE: tabasco_pipeline.pl <genome sequence file> <reference transcript seqs file> <threshold> <outputDir>

# get all the transcript id, since some transcript has no alignment match, so maybe not present in the sam file
my $genome_sequence_file=$ARGV[0];
my $basic_name=$genome_sequence_file;
$basic_name=~s/\.fasta$//;
$basic_name=~s/\.fa$//;

$basic_name=~s/.*\///g;
my $outputDir=$ARGV[3];
my $samFile = "$outputDir/$basic_name" . ".sam";
# Get path of transcript seqs
my $transcript_seqs=$ARGV[1];

system("/programs/minimap2-2.17/minimap2 -ax splice -a -uf -C1 -k 12 -t 12 $genome_sequence_file $transcript_seqs --cs  > $samFile");

my %counts;
my %fragedCounts;
my $seq="";
my $name="";
my %seqs;
my $maskthreadsHold = $ARGV[2];

open INPUT, $transcript_seqs;
while( my $line=<INPUT> ){
	if( $line=~/^>(\S+)/ ){
		$counts{$1}=0;
		if( length($name)>0 && length($seq)>0 ){
			$seq=~s/\s//g;
            $seqs{$name}=$seq;
		}
		$name=$1;
		$seq="";
	}else{
		$seq = $seq . $line;
	}
}
close INPUT;

if( length($name)>0 && length($seq)>0 ){
	$seq=~s/\s//g;
    $seqs{$name}=$seq;
}

print "fasta file reading done\n";

open INPUT, "$samFile";
while( my $line=<INPUT> ){
	if($line=~/^(?!@)/){
		if( $line=~/^(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t(.*?)\t.*NM:i:(\d+).*cs:(\S+)/ ){
			if( $3 ne "*" ){
				my $transcriptId=$1;
				my $maskLength=0;
				my $ciga = $6;
				my $noMatch = $7;
				my $cs = $8;
				while( $ciga=~/(\d+)[SH]/g ){
					$maskLength+=$1;
				}
				my $match = 0;
				while( $cs=~/:(\d+)/g ){
					$match+=$1;
				}
				if( exists $seqs{$transcriptId} ){
					if( ($match)/(length($seqs{$transcriptId})) > $maskthreadsHold ){
						if( exists $counts{$transcriptId} ){
							$counts{$transcriptId}=$counts{$transcriptId}+1;
						}else{
							$counts{$transcriptId}=1;
						}
					}else{
						if( exists $fragedCounts{$transcriptId} ){
							$fragedCounts{$transcriptId}=$fragedCounts{$transcriptId}+1;
						}else{
							$fragedCounts{$transcriptId}=1;
						}
					}
				}else{
					print "there is some wrong with $line\n";
				}
			}
		}
	}
}
close INPUT;

print "sam file reading done\n";

my %everOccrur;
while( my ($key, $value) = each %counts ){
	if( $value>0 ){
		if( $key=~/^(.*?)_(.*?)_(.*?)$/ ){
			$everOccrur{$1}=1;
		}
	}
}

my $occurSize = keys %everOccrur;
print "occurSize\t$occurSize\n";

my %fragedHash;
while( my ($key, $value) = each %fragedCounts ){
	if( $key=~/^(.*?)_(.*?)_(.*?)$/ ){
		if ( exists  $everOccrur{$1} ){

		}else{
			$fragedHash{$1}=1;
		}
	}
}

my $complete=0;
my $missing=0;
my $fraged=keys %fragedHash;
my $duplicated=0;

my %finalCounts;
while( my ($key, $value) = each %counts ){
	if( $key=~/^(.*?)_(.*?)_(.*?)$/ ){
		if( exists $finalCounts{$1} ){

		}else{
			$finalCounts{$1}=0;
		}
		if( exists $everOccrur{$1} ){
			if( 0 == $value ){
				$value = 1;
			}
		}
		if( $value > $finalCounts{$1} ){
			$finalCounts{$1} = $value;
		}
	}
}

open OUTPUT, ">$outputDir/$basic_name\.tabasco_summary";
while( my ($key, $value) = each %finalCounts ){
	if( $value==1 ){
		$complete += 1;
	}elsif( $value==0 ){
		$missing += 1;
	}else{
		$duplicated += 1;
	}
	print OUTPUT "$key\t$value\n";
}
close OUTPUT;

$missing=$missing-$fraged;
open OUTPUT, ">$outputDir/$basic_name\.count";
print OUTPUT "complete\tduplicated\tfraged\tmissing\n";
print OUTPUT "$complete\t$duplicated\t$fraged\t$missing\n";
close OUTPUT;
