#!/usr/bin/perl
use mygff;
open infi,"<","info.txt" or die "open file open...\n";
	while(<infi>){
		chomp;
		my ($name,$gene)=getName2Gene($_);
		print "$name,$gene\n";
	}
		