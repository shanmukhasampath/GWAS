#!/usr/bin/perl

my $filedir = $ARGV[0];
open (my $fh,"<",$filedir) or die "Can't open the file: $filedir";
my %snpdict;
my %posdict;
while (my $line = <$fh>){
    chomp ($line);
    my @values = split(" ", $line);
	$position = join("-",$values[0],$values[3]);
	$snpdict{$values[1]} = $position;
	$posdict{$position} = $values[1];
}
close $fh;

my @dir = split("/",$filedir);
my @filename;
if (index($dir[$#dir], '_') != -1) {
	@filename = split('_',$dir[$#dir]);
} else {
	@filename = split('\.',$dir[$#dir]);
}
my $direct = join('/', @dir[0..$#dir-1]);
$commonsnps = $direct.'/'.$filename[0].'_common_snps.txt';
$updatesnps = $direct.'/'.$filename[0].'_update_snps.txt';

my $filewithdir = $ARGV[1];
open (my $fh,"<",$filewithdir) or die "Can't open the file: $filewithdir";
open (my $cfh,">",$commonsnps) or die "Can't open the file $commonsnps: ";
open (my $ufh,">",$updatesnps) or die "Can't open the file $updatesnps: ";
while (my $line = <$fh>){
	chomp ($line);
	my @values = split(" ",$line);
	my $pos = join("-",$values[0],$values[3]);
	if( exists($snpdict{$values[1]}) and exists($posdict{$pos}) and length($values[4]) == 1 and length($values[5]) == 1 ){
		print $cfh $values[1],"\n";
	} 
	elsif( exists($posdict{$pos}) and length($values[4]) == 1 and length($values[5]) == 1 ){
		print $ufh $posdict{$pos},"\t",$values[1],"\n";
	}
}
close $fh;
close $cfh;
close $ufh;
