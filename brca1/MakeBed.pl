#/usr/bin/perl -w 

use strict;
use warnings; 

open (my $data1, '<', 'start.txt')or die "Could not open file";
open (my $data2, '<', 'stop.txt')or die "Could not open file";
my $outfile1='NM_007294.bed';
open(FILE1,">> $outfile1") or die "problem opening file";

my @fields1;
my @fields2;

while(my $line= <$data1>)  
{ 
	chomp $line;
	@fields1=split "," , $line;	
}
while(my $line= <$data2>)  
{ 
	chomp $line;
	@fields2=split "," , $line;	
}

my $num=scalar @fields1;
for(my $i=0; $i<$num; $i++)
{
	print FILE1 "chr17\t$fields1[$i]\t$fields2[$i]\tNM_007294\tNA\t-\n";
}


