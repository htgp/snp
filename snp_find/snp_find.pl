#!/usr/bin/perl
#use warnings;
#use strict;
use Text::NSP::Measures::2D::Fisher::twotailed;
use Text::NSP::Measures::2D::CHI::x2;
use Algorithm::MinMax;
use Statistics::Distributions;
undef %selectsnp,%selectin;
my %selectsnp,%selectin;
&get_snp;

sub get_snp{
	open snp,"<",$ARGV[0] or die "open snp file failed...\n";
	open out_snp,">>out_snp.txt" or die  "open out-snp file failed...\n";
	while(<snp>){
		print "***************************************************************\n";
		if(!/^CHROM/){
		chomp;
		$_=~tr/\"//d;
		$_=~s/NA/0,0/g;
		#sample1,2,3 and each with 2 replictes.
		# my ($chr,$pos,$ref,$alt,$s1_1,$s1_2,$s2_1,$s2_2,$s3_1,$s3_2,$gid)=split;
		#this is a tmp implementation for sample1,2 each with only 1 replicates;
		my ($chr,$pos,$ref,$alt,$s1_1,$s2_1,$gid)=split;
		#	print "$s1_1\t$s1_2\t$s2_1\t$s2_2\t$s3_1\t$s3_2\n";
	#note here, for sample($s1_1,$s1_2|$s2_1,$s2_2|$s3_1,$s3_2), the original script analyze $s1,$s3 as mut,wild. so if you are going to analyze  other sample, i. e. s2, you had to use s2 to replace @hi1, @hi2;
	#####################################################################################
	#folowing is the original line which is commented.
	#	my @hi1=split(',',$s1_1);my @hi2=split(',',$s1_2);
	#####################################################################################
		my @hi1=split(',',$s1_1);my @hi2=split(',',$s1_2);     #always remember mod here!
		my @wi1=split(',',$s2_1);my @wi2=split(',',$s2_2);     #key code...
	#####################################################################################
	#	print join(',',@hi1,@hi2,@wi1,@wi2),"\n";
		my $hi_ref=&sum($hi1[0],$hi2[0]);
		my $hi_alt=&sum($hi1[1],$hi2[1]);
		my $wi_ref=&sum($wi1[0],$wil[0]);
		my $wi_alt=&sum($wi1[1],$wi2[1]);
		print "$hi_ref,$hi_alt,$wi_ref,$wi_alt\n";
	
	          # word2   ~word2
  # word1    n11      n12 | n1p
 # ~word1    n21      n22 | n2p
           # --------------
           # np1      np2   npp
									  # 以下数组均包含以下四个参数：
									  # n11=>$n11,
                                      # n1p=>$n1p,
                                      # np1=>$np1,
                                      # npp=>$npp);
	@judge_wild_mut_sig=($hi_ref,$hi_ref+$hi_alt,$hi_ref+$wi_ref,$hi_ref+$hi_alt+$wi_ref+$wi_alt);
	@judge_heter=($hi_ref,$hi_ref+$hi_alt,3/2*$hi_ref+1/2*$hi_alt,2*($hi_ref+$hi_alt));
	@judge_ref_all=($hi_ref,$hi_ref+$hi_alt,2*$hi_ref+$hi_alt,2*($hi_ref+$hi_alt));
	@judge_alt_all=($hi_ref,$hi_ref+$hi_alt,$hi_ref,2*($hi_ref+$hi_alt));
	@judge_heter2=($wi_ref,$wi_ref+$wi_alt,3/2*$wi_ref+1/2*$wi_alt,2*($wi_ref+$wi_alt));
	@judge_ref_all2=($wi_ref,$wi_ref+$wi_alt,2*$wi_ref+$wi_alt,2*($wi_ref+$wi_alt));
	@judge_alt_all2=($wi_ref,$wi_ref+$wi_alt,$wi_ref,2*($wi_ref+$wi_alt));
	print $chr . "_" . $pos . "_" . $ref . "_" . $alt,":","$hi_ref,$hi_alt,$wi_ref,$wi_alt.\n";
#filt low coverage;
	if($hi_ref+$hi_alt<15 && $wi_alt+$wi_ref<15){
		print $chr . "_" . $pos . "_" . $ref . "_" . $alt, "->","too low num, skipped...\n";
			}
	else{
		undef $pm,$ph,$pa,$pr,$wi_type,$hi_type;
		print "begin analyze mut...\n";
		$pm=is_wild_mut_sig(@judge_wild_mut_sig);	
		print "Judge if mut-wild have same genotype. pvalue is: $pm.\n";
		#if not wild mut same;
		if($pm<0.05){#null hypersis mut(ref:alt):wild(ref:alt) no difference;
			print "wild-mut different...\n";
			#if is not heterogus;null hypersis ref:alt=1:1;
			$ph=is_heter(@judge_heter);
			print "Judge if mut genotype are ref:alt. pvalue is $ph\n.";
			if($ph<0.05){
				print "Not heterogus\n";
				$pr=is_ref_all(@judge_ref_all);#null hypersis ref:ref(all ref);
				$pa=is_alt_all(@judge_alt_all);#null hypersis alt:alt(all alt);
				print "Judge if mut genotype are  alt:alt,pvalue is $pr; ref:ref, $pa.\n";

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
		print out_snp "find real difference: ", $chr . "_" . $pos . "_" . $ref . "_" . $alt, "->", $hi_type . "_" . $wi_type,"\n" if $hi_type ne $wi_type;
		undef $alt,$ref,$wi_type,$hi_type,$pm,$ph,$pa,$pr;
		print "***************************************************************\n";
			}
		} 
	}
return \%selectsnp;
close snp;
}



sub sum(){
	my @ele= @_;
	my $sum=0;
	for my $ele_member(@ele){
#	print $ele_member,"\n";
	$sum += $ele_member;
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


sub is_wild_mut_sig{
	print "in is_wild_mut_sig...\n";
	my @arr=@_;
	my $p_mut;
	my $method=&which_method(@arr);
	print "the result of which_method: $method.\n";
	if( $method eq "chi"){
		print "arguments1: ", join(',',@arr),",using chi-square test.\n";
		$p_mut=&chi_test(@arr);
		if($p_mut<0.06 && $p_mut>0.04){
			print "The pValue is too close to 0.05. Better using fisher exact test.\n";
			$p_mut=fisher_test(@arr);
			}
		}
	else{
		$p_mut=fisher_test(@arr);
		print "method is fisher: $p_mut.\n";
		}
				return $p_mut;
}
sub is_heter{
	my @arr=@_;
	my $p;
	if( &which_method(@arr) eq "chi"){
	print "arguments2: ", join(',',@arr),",using chi-square test.\n";
	$p=&chi_test(@arr);
	if($p<=0.06 && $p>=0.04){
	print "using fisher exact test.\n";
	$p= &fisher_test(@arr)}
	}
	else{$p= &fisher_test(@arr)}
	return $p;
}
sub is_ref_all{
	my @arr=@_;
	my $p;
	if( &which_method(@arr) eq "chi"){
	print "arguments3: ", join(',',@arr),",using chi-square test.\n";
	$p=&chi_test(@arr);
	if($p<=0.06 && $p>=0.04){
	print "using fisher exact test.\n";
	$p= &fisher_test(@arr)}
	}
	else{$p= &fisher_test(@arr)}
	return $p;
}
sub is_alt_all{
	my @arr=@_;
	my $p;
	if( &which_method(@arr) eq "chi"){
	print "arguments4: ", join(',',@arr),",using chi-square test.\n";
	$p=&chi_test(@arr);
	if($p<=0.06 && $p>=0.04){
	print "using fisher exact test.\n";
	$p= &fisher_test(@arr)}
	}
	else{$p= &fisher_test(@arr)}
	return $p;
}
sub which_method{
	my @arr=@_;
	my $mn11=$arr[0];
	my $mn12=$arr[1]-$arr[0];
	my $mn21=$arr[2]-$arr[0];
	my $mn22=$arr[3]-$mn11-$mn12-$mn21;
	my $t11=$arr[1]*$arr[2]/$arr[3];
	my $t12=$arr[1]*($arr[3]-$arr[2])/$arr[3];
	my $t21=$arr[2]*($arr[3]-$arr[1])/$arr[3];
	my $t22=($arr[3]-$arr[2])*($arr[3]-$arr[1])/$arr[3];
	print "actual frequency: $mn11,$mn12,$mn21,$mn22.\n";
	print "thereical frequency: $t11,$t12,$t21,$t22.\n";
	my @a=($t11,$t12,$t21,$t22);
	my @r= Algorithm::MinMax->minmax( \@a );
	my $min=$r[0];
	print "min of TF: $min\n";
	if($min>=1 && $arr[3]>=40){
		$method="chi";
		}
	else{$method="fisher";}
	return $method;
	}
	
sub chi_test{
	my @arr=@_;
	print "data in are: ", join(',',@arr),".\n";
	my $x2_value=Text::NSP::Measures::2D::CHI::x2::calculateStatistic(n11=>$arr[0],n1p=>$arr[1],np1=>$arr[2],npp=>$arr[3]);
	my $pchi=Statistics::Distributions::chisqrprob(1,$x2_value);
	print "X2 and pCHI are $x2_value and $pchi.\n";
	return $pchi;
	}
sub fisher_test{
	my ($n11,$n1p,$np1,$npp)=@_;
	print "in fisher_test process: ",join(',',$n11,$n1p,$np1,$npp),".\n";
	$fisher=Text::NSP::Measures::2D::Fisher::twotailed::calculateStatistic(n11=>$n11,
                                        n1p=>$n1p,
                                        np1=>$np1,
                                        npp=>$npp,
						);
	if( $errorCode = Text::NSP::Measures::2D::Fisher::twotailed::getErrorCode()){
		print $errorCode." - ".getErrorMessage()."\n";
		}
		
	return $fisher;
	}	