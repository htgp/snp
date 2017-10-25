#!/usr/bin/perl
use Array::Compare;
%s1=%s2=%hash1=%hash2=%rev_s1=%rev_s2={};$conc=$aa1=$aa2=$bb1=$bb2=$cc1=$cc2=$dd1=$dd2=""; @a=@b=@mu1=@mu2=@wi1=@wi2=();
open s1, "<", $ARGV[0] or die "open file 1 failed...\n";
while(<s1>){
		chomp;
		($tip,$info)=split(': ',$_);
		($key,$nt)=split('->',$info);
		($m1,$m2,$w1,$w2)=split('_',$nt);
		# print "Sample1:$key,$m1,$m2,$w1,$w2\n";
		$s1{$key}=[$m1,$m2,$w1,$w2];
	}
	close s1;


	
open s2, "<", $ARGV[1] or die "open file 2 failed...\n";
while(<s2>){
		chomp;
		($tip,$info)=split(': ',$_);
		($key,$nt)=split('->',$info);
		($m1,$m2,$w1,$w2)=split('_',$nt);
		# print "Sample2: $key,$m1,$m2,$w1,$w2\n";
		$s2{$key}=[$m1,$m2,$w1,$w2];
	}
	close s2;

	
# $hr=\%s1;
# &test_hash($hr);
sub test_hash{
	($ref)=@_;
	print "accept arguments:$ref.\n";
	%h=%{$ref};
	for $key(keys %h){
			print "$key:";
			@arr=@{$h{$key}};
			print join(',',@arr);
			print "\n";	
	}
}


&m_uniq(\%s1,\%s2);
	
sub m_uniq{
	($hr1,$hr2)=@_;
	%hash1=%{$hr1};
	%hash2=%{$hr2};
	undef $diff;
	open m1u, ">>m1u.txt" or die "open m1u failed...\n";
	print m1u "#Site:(s1_alt1)_(s1_alt2):(s2_alt1):(s2_alt2)\n";
	for $key(keys %hash1){
		print "In sample1,";
		$diff=&equal4nt($hr1,$hr2,$key);
		# print "diff:$diff\n";
		if($diff ne ""){
			print  m1u "$key->$diff\n";
			}
		}
		close m1u;
	undef $diff;
	open m2u,">>m2u.txt" or die "open m2u failed...\n";
	print m2u "#Site:(s1_alt1)_(s1_alt2):(s2_alt1):(s2_alt2)\n";
	for $key(keys %hash2){
		print "In sample2,";
		$diff=&equal4nt($hr1,$hr2,$key);
		# print "diff:$diff\n";
		if($diff ne ""){
			print  m2u "$key->$diff\n";
			}
		}
	}
		close m2u;
sub equal4nt{
		($ref1,$ref2,$key)=@_;
		$conc=$aa1=$aa2=$bb1=$bb2=$cc1=$cc2=$dd1=$dd2=""; @mu1=@mu2=@wi1=@wi2=();
		$comp=Array::Compare->new;
		%rev_s1=%{$ref1};
		%rev_s2=%{$ref2};
		if(!exists $rev_s1{$key}){
			($aa2,$bb2,$cc2,$dd2)=@{$rev_s2{$key}};
			# print "aa2,bb2,cc2,dd2:$aa2,$bb2,$cc2,$dd2.\n";
			@mu2=($aa2,$bb2);
			@wi2=($cc2,$dd2);
			@mu1=($cc2,$dd2);
			@wi1=($cc2,$dd2);
			}
		if(!exists $rev_s2{$key}){
			($aa1,$bb1,$cc1,$dd1)=@{$rev_s1{$key}};
			@mu1=($aa1,$bb1);
			@wi1=($cc1,$dd1);
			@mu2=($cc1,$dd1);
			@wi2=($cc1,$dd1);
			}
		if(exists $rev_s2{$key} && exists $rev_s1{$key}){
				($aa2,$bb2,$cc2,$dd2)=@{$rev_s2{$key}};
				($aa1,$bb1,$cc1,$dd1)=@{$rev_s1{$key}};
				@mu1=($aa1,$bb1);
				@mu2=($aa2,$bb2);
				@wi1=($cc1,$dd1);
				@wi2=($cc2,$dd2);
			}
		# print "aa1,bb1,cc1,dd1:$aa1,$bb1,$cc1,$dd1\n";
		# print "array mu1:mu2:wi1:wi2--->", join("/",@mu1) . ":" . join("/",@mu2) .  ":" .  join("/",@wi1) . ":" . join("/",@wi2),"\n";
		if( ! $comp->compare(\@mu1,\@mu2)  && ! $comp->compare(\@mu1,&swap_arr_ele(@mu2))){
			$conc= join("/",@mu1) . ":" . join("/",@mu2);
			print "$key:$conc\n";
			return $conc;
			}
		else{
			print "same snp for $key.\n";
			$conc="";
			return $conc;
		}
		 $aa1=$aa2=$bb1=$bb2=$cc1=$cc2=$dd1=$dd2=""; @mu1=@mu2=@wi1=@wi2=();
		}
			
			
		
		
		
sub verify_same_arr{
	($mua,$mub,$muc,$mud)=@_;
	@a=($mua,$mub);
	@b=($muc,$mud);
	#return 1 if same;
	$comp=Array::Compare->new;
	if($comp->compare(\@a,\@b)){
		return 1;
		}
	else{return 0;}
}	

sub swap_arr_ele{
	($e1,$e2)=@_;
	@ret=($e2,$e1);
	return \@ret;
}
		