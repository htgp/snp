#!/usr/bin/perl
use snp_test;

$fa="sc61.fa";
$gff="sc61_rome.gff";
$snp="m1_uniq2_m2.txt";
# $snp="test.txt";
%gff=hash_gff($gff);
$gref=\%gff;

&loader;

sub loader{
open snp,"<","$snp" or die "open snp file failed...\n";
while(<snp>){#this snp file is like:I_106385_C_T->T/T:C/T;
	if(!/^#/){
		chomp;
		my ($key,$comp)=split('->',$_);
		my ($chr,$pos,$change)=split('_',$key);
		my @snps=($m1a,$m1b,$m2a,$m2b)=split('/|:',$comp);
		if(&is_snp(@snps)==1){
		my ($outaas,$gnm,$aa_ab_pos)=&snp2aa($chr,$pos,\@snps,$gref); 
		($m1a2aa,$m1b2aa,$m2a2aa,$m2b2aa)=@{$outaas};
		print "$chr:$pos->$m1a/$m1b:$m2a/$m2b,($m1a2aa/$m1b2aa:$m2a2aa/$m2b2aa) at site $aa_ab_pos of protein $gnm:";
		if($m1a2aa eq $m1b2aa && $m2a2aa eq $m1b2aa && $m2a2aa eq $m2b2aa){
			print "nonsense\n" ;
			}
		else{print "sense\n";
				}
			}
		else{#indel;
			my ($gnm,$type)=&nc_location($chr,$pos,$m1a,$gref);
			if($type eq "cds"){
			print "$chr:$pos->$m1a/$m1b:$m2a/$m2b,indel induced mutation at site of $pos of protein $gnm.sense\n";
			}
			elsif($type=~/RNA/){print "$chr:$pos->$m1a/$m1b:$m2a/$m2b,indel induced mutation at site of $pos of protein $gnm.\n";}
			else{print "$chr:$pos->$m1a/$m1b:$m2a/$m2b,indel induced mutation at site of $pos of $type $gnm that will be ignored.\n";}
			}
		}
		}
	close snp;
		
	}		