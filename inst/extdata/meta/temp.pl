#!/usr/local/bin/perl 

@files = glob ("*.doc");

foreach $i (@files) {
$base = substr($i,0,index ($j,".pgm")-3);
system("mv $i $base.txt");
print "\n$base.txt"; 
}
