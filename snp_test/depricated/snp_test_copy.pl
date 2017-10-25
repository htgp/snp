#!/usr/bin/perl
BEGIN{
	unshift(@INC,".");
	}
use snp_test;

$fa="sc61.fa";
$gff="sc61_rome.gff";
$snp="tmp.txt";
%gff=hash_gff($gff);
$gref=\%gff;
# foreach my $key(keys %$gref){
	# print $key,"\n";
	# }
&loader;

sub loader{
open snp,"<","$snp" or die "open snp file failed...\n";
while(<snp>){#this snp file is like:I_106385_C_T->T/T:C/T;
	if(!/^#/){
		chomp;
		my ($key,$comp)=split('->',$_);
		$key=(split(' ',$key))[-1];
		my ($chr,$pos,$ref,$alt)=split('_',$key);
		# print join(' ',($chr,$pos,$ref,$alt)),"\n";
		# print $comp,"\n";
		my @snps=($m1a,$m1b,$m2a,$m2b)=split('_',$comp);
		# print join(',',@snps),"\n";
		if(&is_snp(@snps)==1){
				my ($outaas,$gnm,$aa_ab_pos,$type)=&snp2aa($chr,$pos,\@snps,$gref); 
				print "mark1:@$outaas,$gnm,$aa_ab_pos;\n";
			if($outaas ne "other" && $gnm ne "other" && $aa_ab_pos ne "other"){
				($m1a2aa,$m1b2aa,$m2a2aa,$m2b2aa)=@{$outaas};
				print "mark2:$chr:$pos->$m1a/$m1b:$m2a/$m2b,($m1a2aa/$m1b2aa:$m2a2aa/$m2b2aa) at site $aa_ab_pos of protein $gnm:";
				if($m1a2aa eq $m1b2aa && $m2a2aa eq $m1b2aa && $m2a2aa eq $m2b2aa){
					print "nonsense\n" ;
				}
				else{print "sense\n";
				}
			}
			else{print "mark3:$chr:$pos->$m1a/$m1b:$m2a/$m2b,snp induced mutation at site of $pos of a $type.ignored.\n";}
			}
		else{#indel;
			my ($gnm,$type)=&nc_location($chr,$pos,\@snps,$gref);
			if($gnm ne "other"){
			# print "type:$type\n";
			if($type eq "cds"){
			print "mark4:$chr:$pos->$m1a/$m1b:$m2a/$m2b,indel induced mutation at site of $pos of protein $gnm.sense\n";
			}
			elsif($type=~/RNA/){print "mark5:$chr:$pos->$m1a/$m1b:$m2a/$m2b,indel induced mutation at site of $pos of protein  $type $gnm.\n";}
			else{print "mark6:$chr:$pos->$m1a/$m1b:$m2a/$m2b,indel induced mutation at site of $pos of a $type that will be ignored.\n";}
			}
		}
		}#if
		}#while
	close snp;
		
	}		