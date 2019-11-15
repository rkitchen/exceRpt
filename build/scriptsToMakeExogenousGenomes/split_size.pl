#!/usr/bin/perl
#####################################
# Program: split_size.pl  -  Date: Mon Feb  9 11:32:24 EST 2015
# Autor: Fabio C. P. Navarro - Ludwig
# Goal:
#
# Input:
#
# Output:
#
#####################################


use Getopt::Long;
use strict;
use warnings;

my ( $file, $file2, $usage ) = ( "", "", 0 );
GetOptions( "f|file=s"      => \$file,
#           "f2|file2=s"      => \$file2,
           "h|help|usage"  => \$usage )
    || die "Error while parsing comand line arguments";

usage() if( $usage );

if( $file eq "" ) {
    usage();
}

my $size = 7.5*1024*1024*1024;
my $suffix = 0;
my $currsize = 0;
my $temp = "";
my $OUT;

open ( $OUT, '>', "$file.$suffix");
open ( IN0, "<$file" );
while ( <IN0> ) {
	#chomp $_;
	#my @tokens = split ( /[ \t\cI]+/,$_ );
	if ( $_ =~ /^[\>]/ ) {
		$currsize += length($temp);
		if ( $size - $currsize < 0)  { 
			$suffix++;
			close ($OUT);
			open ( $OUT, '>', "$file.$suffix");
			$currsize = 0;
		}
		print $OUT $temp; 
		$temp = $_;
	}
	else {
		$temp .= $_;
	}
}
print $OUT $temp;
close($OUT);

#open ( IN1, "<$file2" );
#while ( <IN1> ) {
#	chomp $_;
#	my @tokens = split ( /[ \t\cI]+/,$_ );
#}



sub usage {
    die "Usage: $0 [options]
Available Options :
   -f  | --file         : Filename in.
   -h  | --help --usage : Print this message.
";

}
