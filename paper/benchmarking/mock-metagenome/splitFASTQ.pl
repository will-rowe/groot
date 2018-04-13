use strict;
use warnings;
my $file = $ARGV[0];
open(FILE, "<$file") || die "cannot open $file\n";
open(OUT1, ">$file\_1") || die "cannot open $file\_1\n";
open(OUT2, ">$file\_2") || die "cannot open $file\_2\n";
while(<FILE>){
    chomp;
    print OUT1 "$_\/1\n";
    print OUT2 "$_\/2\n";
    my $newline = <FILE>; chomp($newline);
    print OUT1 substr($newline, 0, length($newline)/2)."\n";
    print OUT2 substr($newline, length($newline)/2, length($newline)/2)."\n";
    $newline = <FILE>; chomp($newline);
    print OUT1 "$newline\/1\n";
    print OUT2 "$newline\/2\n";
    $newline = <FILE>; chomp($newline);
    print OUT1 substr($newline, 0, length($newline)/2)."\n";
    print OUT2 substr($newline, length($newline)/2, length($newline)/2)."\n";
}
close(FILE);
