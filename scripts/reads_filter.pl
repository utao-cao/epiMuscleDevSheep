#!/usr/bin/perl -w
use Getopt::Long;
use Cwd qw/abs_path/;
use List::Util qw/max min/;

# Version : v1.0
# Data    : 2020-03-13

my ( $in );
GetOptions(
	'i:s' => \$in,
	'o:s' => \$out,
	'len:i' => \$len,
);
=cut
if (defined $help || !defined $file || !defined $id || !defined $out) {
	print "#############################################################\n";
	print "reads filter\n";
	print "Usage:\n";
	print "perl reads_filter.pl -i ${sample}.sort.bed -o ${sample}.bed.filter -len 150\n";
	print "Options:\n";
	print "  -i   input bed filename\n";
	print "  -o   output bed filename\n";
	print "  -len insert size threshold\n";
	print "#############################################################\n";
	exit;
}
=cut


open (OUT,">$out");
open (IN,"<$in");
while (my $line1 = <IN>) {
	chomp $line1;
	my $line2 = <IN>;
	chomp $line2;
	my @a1 = split(/\t/,$line1);
	my @a2 = split(/\t/,$line2);
	my @b;
	push(@b, ($a1[1],$a1[2],$a2[1],$a2[2]) );
	my $max = max @b;
	my $min = min @b;
	my $insert = $max - $min;
	if ($insert <= $len) {
		print OUT "$line1\n$line2\n";
	}
}
close IN;
close OUT;


