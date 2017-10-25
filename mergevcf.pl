#!/usr/bin/perl
@files=@ARGV;
undef %index,%has,%snp;
print "building index......\n";
 foreach(@files){
	$strain=(split("_",$_))[0];
	open vcff, "<",$_ or die "open $_ failded...\n";
	while(<vcff>){
		if(!/^#/){
			$bak=$_;
			chomp;
			($chr,$pos,$ref,$alt,$genotype)=(split)[0,1,3,4,9];
			($gt,$gl,$score)=(split(":",$genotype))[0,1,2];
			$alt=(split /,/,$alt)[0];
			# if($ref=~/\b[ATCG]\b/ && ($alt=~/\b[ATCG]\b/ || $alt=~/\b[ATCG],[ATCG]\b/ || $alt=~/\b[ATCG],[ATCG],[ATCG]\b/)){$snp{$site}{'snp'}=1;}
			# else{$snp{$site}{'snp'}=0;}
			$site=$chr . "_" . $pos;
			$index{$site}=$ref;
			$has{$strain}{$site}=1;
			$snp{$site}{$strain}{'alt'}=$alt;
			$snp{$site}{$strain}{'ref'}=$ref;
			$snp{$site}{$strain}{'gt'}=$gt;
					}
				}
	close vcff;
	}
# print "testing hash...\n";
# $test="V_225841";$strain="F10";
# print join "*", $index{$site},$has{$strain}{$site},$snp{$site}{$strain}{'alt'},$snp{$site}{$strain}{'ref'},$snp{$site}{$strain}{'gt'},$snp{$site}{$strain}{'snp'};
# print "\n";

print "listing snps...\n";
print "*" x 100;
print "\n";
print "Chr\tPos\tF10\tF30\tF50\tF13\tF14\n";
foreach $site(keys %index){
	($chr,$pos)=(split(/_/,$site))[0,1];
	print "$chr\t$pos\t";
	 foreach $strain(qw(F10 F30 F50 F13 F14)){
		if($has{$strain}{$site} ==1){
			if($snp{$site}{$strain}{'gt'} eq "1/1"){
				print $snp{$site}{$strain}{'alt'} . "|" . $snp{$site}{$strain}{'alt'};
				print "\t";
				}
			elsif($snp{$site}{$strain}{'gt'} eq "0/1"){
				print $snp{$site}{$strain}{'ref'} . "|" . $snp{$site}{$strain}{'alt'};
				print "\t";
				}
			else{print "Bad_gt";print "\t";};
			print "\t";
		}
		else{print $index{$site} . "|" . $index{$site};print "\t";}
		}
		# print $snp{$site}{'snp'};
		print "\n";
	}
		