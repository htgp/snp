#!/usr/bin/perl
#use warnings;
#use strict;
use Text::NSP::Measures::2D::Fisher::twotailed;
undef %selectsnp,%selectin;
our %selectsnp,%selectin;
&get_snp;
&get_indel;
#$snp_hash_ref=&get_snp;
#$indel_hash_ref=&get_indel;
#&print_hash($snp_hashi_ref);
#&print_hash($indel_ref_hash);


sub get_snp{
open snp,"<",$ARGV[0] or die "open snp file failed...\n";
open out_snp_tmp,">>out_snp_hi.txt" or die  "open out-snp file failed...\n";
while(<snp>){
	print "***************************************************************\n";
	if(!/^CHROM/){
	chomp;
	$_=~tr/\"//d;
	$_=~s/NA/0,0/g;
	my ($chr,$pos,$ref,$alt,$s1_1,$s1_2,$s2_1,$s2_2,$s3_1,$s3_2,$gid)=split;
	&ref_alt_alt;
#	print "$s1_1\t$s1_2\t$s2_1\t$s2_2\t$s3_1\t$s3_2\n";
#note here, for sample($s1_1,$s1_2|$s2_1,$s2_2|$s3_1,$s3_2), the original script analyze $s1,$s3 as mut,wild. so if you are going to analyze  other sample, i. e. s2, you had to use s2 to replace @hi1, @hi2;
#####################################################################################
#folowing is the original line which is commented.
	my @hi1=split(',',$s1_1);my @hi2=split(',',$s1_2);
#####################################################################################
	# my @hi1=split(',',$s2_1);my @hi2=split(',',$s2_2);
	my @wi1=split(',',$s3_1);my @wi2=split(',',$s3_2);
#	print join(',',@hi1,@hi2,@wi1,@wi2),"\n";
	my $hi_ref=&sum($hi1[0],$hi2[0]);
	my $hi_alt=&sum($hi1[1],$hi2[1]);
	my $wi_ref=&sum($wi1[0],$wil[0]);
	my $wi_alt=&sum($wi1[1],$wi2[1]);
	# print "$hi_ref,$hi_alt,$wi_ref,$wi_alt\n";
	@judge_wild_mut_sig=($hi_ref,$hi_ref+$hi_alt,$hi_ref+$wi_ref,$hi_ref+$hi_alt+$wi_ref+$wi_alt);
	@judge_heter=($hi_ref,$hi_ref+$hi_alt,3/2*$hi_ref+1/2*$hi_alt,2*($hi_ref+$hi_alt));
	@judge_ref_all=($hi_ref,$hi_ref+$hi_alt,2*$hi_ref+$hi_alt,2*($hi_ref+$hi_alt));
	@judge_alt_all=($hi_ref,$hi_ref+$hi_alt,$hi_ref,2*($hi_ref+$hi_alt));
	@judge_heter2=($wi_ref,$wi_ref+$wi_alt,3/2*$wi_ref+1/2*$wi_alt,2*($wi_ref+$wi_alt));
	@judge_ref_all2=($wi_ref,$wi_ref+$wi_alt,2*$wi_ref+$wi_alt,2*($wi_ref+$wi_alt));
	@judge_alt_all2=($wi_ref,$wi_ref+$wi_alt,$wi_ref,2*($wi_ref+$wi_alt));
	print $chr . "_" . $pos . "_" . $ref . "_" . $alt,":","$hi_ref,$hi_alt,$wi_ref,$wi_alt.\n";
#filt low coverage;
	if($hi_ref+$hi_alt<15 || $wi_alt+$wi_ref<15){
		print $chr . "_" . $pos . "_" . $ref . "_" . $alt, "->","too low num, skipped...\n";
			}
	else{
		undef $pm,$ph,$pa,$pr,$wi_type,$hi_type;
		print "begin analyze mut...\n";
		$pm=is_wild_mut_sig(@judge_wild_mut_sig);	
		print "Judge if mut-wild have same genotype. pvalue is $pm.";
		#if not wild mut same;
		if($pm<0.05){#null hypersis mut(ref:alt):wild(ref:alt) no difference;
			print "wild-mut different...\n";
			#if is not heterogus;null hypersis ref:alt=1:1;
			$ph=is_heter(@judge_heter);
			print "Judge if mut genotype are ref:alt. pvalue is $ph.";
			if($ph<0.05){
				print "Not heterogus\n";
				$pr=is_ref_all(@judge_ref_all);#null hypersis ref:ref(all ref);
				$pa=is_alt_all(@judge_alt_all);#null hypersis alt:alt(all alt);
				print "Judge if mut genotype are  alt:alt,pvalue is $pr; ref:ref, $pa. ";

				if($pr<0.05){
					print "alt:alt=",$alt . "_" . $alt,".\n";
					$hi_type=$alt . "_" . $alt;
					}
				elsif($pa<0.05){
					if($pa<$pr){
					print "ref/alt=",$ref . "_" . $ref,".\n";
					$hi_type=$ref . "_" . $ref;
					}
					}
				else{print "No properly catagorized...\n";}
						}
			else{
				print "heterogus ref:alt=", $ref . "_" . $alt,". \n";
				$hi_type=$ref . "_" . $alt;
				}
				}
		else{#pm>0.05,not significant;
			print "wild_mut same genotype.  ignored...\n";
			undef $hi_type,$pm,$ph,$pa,$pr;
			next;
			}
	print "begin analyze wild...\n";
		#if not wild mut same;
		undef $pa,$pr,$ph,$wi_type;
		if($pm<0.05){
			# print "wild-mut different...\n";
			#if not heterogus;
			$ph=is_heter(@judge_heter2);
			print "Judge if mut genotype are ref:alt. pvalue is $ph.";
			if($ph<0.05){
				print "Not heterogus\n";
				$pr=is_ref_all(@judge_ref_all2);
				$pa=is_alt_all(@judge_alt_all2);
				print "Judge if mut genotype are  alt:alt,pvalue is $pr; ref:ref, $pa. ";
				if($pr<0.05){
					print "alt:alt=",$alt . "_" . $alt,".\n";
					$wi_type=$alt . "_" . $alt;
					}
				elsif($pa<0.05){
					if($pa<$pr){
					print "ref/alt=",$ref . "_" . $ref,".\n";
					$wi_type=$ref . "_" . $ref;
					}
					}
				else{print "No properly catagorized...\n";}
						}
			else{
				print "heterogus ref:alt=",$ref . "_" . $alt,". \n";
				$wi_type=$ref . "_" . $alt;
				}
				}
		$selectsnp{$chr . "_" . $pos . "_" . $ref . "_" . $alt}=$hi_type . "_" . $wi_type if $fish<0.05 && $hi_type ne $wi_type;
		print "find real difference: ", $chr . "_" . $pos . "_" . $ref . "_" . $alt, "->", $hi_type . "_" . $wi_type,"\n" if $hi_type ne $wi_type;
		print out_snp_tmp "find real difference: ", $chr . "_" . $pos . "_" . $ref . "_" . $alt, "->", $hi_type . "_" . $wi_type,"\n" if $hi_type ne $wi_type;
		undef $alt,$ref,$wi_type,$hi_type,$pm,$ph,$pa,$pr;
		print "***************************************************************\n";
			}
		} 
	}
return \%selectsnp;
close snp;
}

sub get_indel{
open indel,"<",$ARGV[1] or die "open snp indel failed...\n";
open out_indel_tmp,">>","out_indel_hi.txt" or die "open out-indel failed...\n";
while(<indel>){
		print "***************************************************************\n";
	if(!/^CHROM/){
	chomp;
	$_=~tr/\"//d;
	$_=~s/NA/0,0/g;
	my ($chr,$pos,$ref,$alt,$s1_1,$s1_2,$s2_1,$s2_2,$s3_1,$s3_2,$gid)=split;
#	print "$s1_1\t$s1_2\t$s2_1\t$s2_2\t$s3_1\t$s3_2\n";
#note here, for sample($s1_1,$s1_2|$s2_1,$s2_2|$s3_1,$s3_2), the original script analyze $s1,$s3 as mut,wild. so if you are going to analyze  other sample, i. e. s2, you had to use s2 to replace @hi1, @hi2;
#####################################################################################
#folowing is the original line which is commented.
#	my @hi1=split(',',$s1_1);my @hi2=split(',',$s1_2);
#####################################################################################
	my @hi1=split(',',$s2_1);my @hi2=split(',',$s2_2);
	my @wi1=split(',',$s3_1);my @wi2=split(',',$s3_2);
#	print join(',',@hi1,@hi2,@wi1,@wi2),"\n";
	my $hi_ref=&sum($hi1[0],$hi2[0]);
	my $hi_alt=&sum($hi1[1],$hi2[1]);
	my $wi_ref=&sum($wi1[0],$wil[0]);
	my $wi_alt=&sum($wi1[1],$wi2[1]);
	# print "$hi_ref,$hi_alt,$wi_ref,$wi_alt\n";
	@judge_wild_mut_sig=($hi_ref,$hi_ref+$hi_alt,$hi_ref+$wi_ref,$hi_ref+$hi_alt+$wi_ref+$wi_alt);
	@judge_heter=($hi_ref,$hi_ref+$hi_alt,3/2*$hi_ref+1/2*$hi_alt,2*($hi_ref+$hi_alt));
	@judge_ref_all=($hi_ref,$hi_ref+$hi_alt,2*$hi_ref+$hi_alt,2*($hi_ref+$hi_alt));
	@judge_alt_all=($hi_ref,$hi_ref+$hi_alt,$hi_ref,2*($hi_ref+$hi_alt));
	@judge_heter2=($wi_ref,$wi_ref+$wi_alt,3/2*$wi_ref+1/2*$wi_alt,2*($wi_ref+$wi_alt));
	@judge_ref_all2=($wi_ref,$wi_ref+$wi_alt,2*$wi_ref+$wi_alt,2*($wi_ref+$wi_alt));
	@judge_alt_all2=($wi_ref,$wi_ref+$wi_alt,$wi_ref,2*($wi_ref+$wi_alt));
	print $chr . "_" . $pos . "_" . $ref . "_" . $alt,":","$hi_ref,$hi_alt,$wi_ref,$wi_alt.\n";
#filt low coverage;
	if($hi_ref+$hi_alt<15 || $wi_alt+$wi_ref<15){
		print $chr . "_" . $pos . "_" . $ref . "_" . $alt, "->","too low num, skipped...\n";
			}
	else{
		undef $pm,$ph,$pa,$pr,$wi_type,$hi_type;
		print "begin analyze mut...\n";
		$pm=is_wild_mut_sig(@judge_wild_mut_sig);	
		print "Judge if mut-wild have same genotype. pvalue is $pm.";
		#if not wild mut same;
		if($pm<0.05){#null hypersis mut(ref:alt):wild(ref:alt) no difference;
			print "wild-mut different...\n";
			#if is not heterogus;null hypersis ref:alt=1:1;
			$ph=is_heter(@judge_heter);
			print "Judge if mut genotype are ref:alt. pvalue is $ph.";
			if($ph<0.05){
				print "Not heterogus\n";
				$pr=is_ref_all(@judge_ref_all);#null hypersis ref:ref(all ref);
				$pa=is_alt_all(@judge_alt_all);#null hypersis alt:alt(all alt);
				print "Judge if mut genotype are  alt:alt,pvalue is $pr; ref:ref, $pa. ";

				if($pr<0.05){
					print "alt:alt=",$alt . "_" . $alt,".\n";
					$hi_type=$alt . "_" . $alt;
					}
				elsif($pa<0.05){
					if($pa<$pr){
					print "ref/alt=",$ref . "_" . $ref,".\n";
					$hi_type=$ref . "_" . $ref;
					}
					}
				else{print "No properly catagorized...\n";}
						}
			else{
				print "heterogus ref:alt=", $ref . "_" . $alt,". \n";
				$hi_type=$ref . "_" . $alt;
				}
				}
		else{#pm>0.05,not significant;
			print "wild_mut same genotype.  ignored...\n";
			undef $hi_type,$pm,$ph,$pa,$pr;
			next;
			}
	print "begin analyze wild...\n";
		#if not wild mut same;
		if($pm<0.05){
			# print "wild-mut different...\n";
			#if not heterogus;
			$ph=is_heter(@judge_heter2);
			print "Judge if mut genotype are ref:alt. pvalue is $ph.";
			if($ph<0.05){
				print "Not heterogus\n";
				$pr=is_ref_all(@judge_ref_all2);
				$pa=is_alt_all(@judge_alt_all2);
				print "Judge if mut genotype are  alt:alt,pvalue is $pr; ref:ref, $pa. ";
				if($pr<0.05){
					print "alt:alt=",$alt . "_" . $alt,".\n";
					$wi_type=$alt . "_" . $alt;
					}
				elsif($pa<0.05){
					if($pa<$pr){
					print "ref/alt=",$ref . "_" . $ref,".\n";
					$wi_type=$ref . "_" . $ref;
					}
					}
				else{print "No properly catagorized...\n";}
						}
			else{
				print "heterogus ref:alt=",$ref . "_" . $alt,". \n";
				$wi_type=$ref . "_" . $alt;
				}
				}
		$selectin{$chr . "_" . $pos . "_" . $ref . "_" . $alt}=$hi_type . "_" . $wi_type if $fish<0.05 && $hi_type ne $wi_type;
		print "find real difference: ", $chr . "_" . $pos . "_" . $ref . "_" . $alt, "->", $hi_type . "_" . $wi_type,"\n" if $hi_type ne $wi_type;
		print out_indel_tmp "find real difference: ", $chr . "_" . $pos . "_" . $ref . "_" . $alt, "->", $hi_type . "_" . $wi_type,"\n" if $hi_type ne $wi_type;
		undef $alt,$ref,$wi_type,$hi_type,$pm,$ph,$pa,$pr;
		print "***************************************************************\n";
			}
		}  
	}	
return \%selectin;
close indel;
}

sub get_gff(){
open gff,"<",$ARGV[2] or die "open gff failed...\n";
while(<gff>){
	if(!/^#/){
		chomp;
		split;	
}}
}

sub sum(){
	@ele= @_;
	$sum=0;
	for my $ele_member(@ele){
#	print $ele_member,"\n";
	our $sum += $ele_member;
#	print "now sum is: ",$sum,"\n";
	}
return $sum;
}

sub print_hash(){
	($ref)=@_;
	%hash=%{$ref};
	for $k(keys %hash){
		print "$k\t$hash{$k}\n";
	}
}
sub fisher_test{
	($n11,$n1p,$np1,$npp)=@_;
	print "in fisher_test process: $n11,$n1p,$np1,$npp.\n";
	$fisher=calculateStatistic(n11=>$n11,
                                        n1p=>$n1p,
                                        np1=>$np1,
                                        npp=>$npp,
						);
	return $fisher;
}
sub fisher_test_ref{
		($ref)=@_;
		($n11,$n1p,$np1,$npp)=@{$ref};
		print "arguments got: $n11,$n1p,$np1,$npp.\n";
		$fisher=calculateStatistic(n11=>$n11,
                                        n1p=>$n1p,
                                        np1=>$np1,
                                        npp=>$npp,
						);	
		return $fisher;
}
sub is_wild_mut_sig{
	@arr=@_;
	print "arguments1: ", join(',',@arr),"\n";
	$p=fisher_test(@arr);
	return $p;
}
sub is_heter{
	@arr=@_;
	print "arguments2: ", join(',',@arr),"\n";
	$p=fisher_test(@arr);
	return $p;
}
sub is_ref_all{
	@arr=@_;
	print "arguments3: ", join(',',@arr),"\n";
	$p=fisher_test(@arr);
	return $p;
}
sub is_alt_all{
	@arr=@_;
	print "arguments4: ", join(',',@arr),"\n";
	$p=fisher_test(@arr);
	return $p;
}
sub ref_alt_alt{
	if($alt=~/,/){
		print "this site has two alters: $ref/$alt.\n";
		(alt1,alt2)=(split(',',$alt))[0,1];
		($s1_1,$s1_2,$alt)=adapt($s1_1,$s1_2);
		($s2_1,$s2_2,$alt)=adapt($s2_1,$s2_2);
		($s3_1,$s3_2,$alt)=adapt($s3_1,$s3_2);
		}
		}
sub adapt{
	($t1,$t2)=@_;
	($t1a,$t1b,$t1c)=split(',',$t1);
	($t2a,$t2b,$t2c)=split(',',$t2);
	@judge_complex_ref_alt=($t1a,sum($t1a+$t1b),sum($t1a+$t1b+$t1c),sum($t1a+$t1b+$t1c+$t2a+$t2b+$t2c));
	@judge_complex_alt_alt=($t1a,sum($t1a+$t1b),$t1a,sum($t1a+$t1b+$t1c+$t2a+$t2b+$t2c));
	@judge_complex_ref_ref=($t1a,sum($t1a+$t1b),sum($t1b+$t1c+$t2b+$t2c),sum($t1a+$t1b+$t1c+$t2a+$t2b+$t2c));
	
	if(fisher_test(@judge_complex_ref_alt)>0.05){#ref:alt=1:1?
		return join(',',$t1a,sum($t1b,$t1c)),join(',',$t2a,sum($t2b,$t2c)),$alt1;
		}
	else{
		if(fisher_test(@judge_complex_alt_alt)>0.05){#not ref:alt=1:1;then alt:alt?
			return join(',',$1b,$1c),join(',',$2b,$2c),$alt2;
			$ref=$alt1;
		}
		if(fisher_test(@judge_complex_ref_ref)>0.05){
			if(fisher_test(@judge_complex_ref_ref)>fisher_test(@judge_complex_alt_alt)){
				if(sum($t1b,$t2b)>=sum($t1c,$t2c)){
					return join(',',$t1a,$t1b),join(',',$t2a,$t2b),$alt1;}
				else{return join(',',$t1a,$t1c),join(',',$t2a,$t2c),$alt2;}
			}
			}
	}
	}

		
		
		
		
		
		
		