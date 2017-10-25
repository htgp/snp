#!/usr/bin/perl
use snp_test;

my $gff_file=$ARGV[0];
my $diff_snp_file=$ARGV[1];


%gff=hash_gff($gff_file);
# for $key(keys %gff){print "in main test hash...";print $key," ",$gff{$key},"\n";}

$gref=\%gff;



open snp,"<",$diff_snp_file or die "open snp file failed...\n";
#the diff_snp_file is a result of snp_find.pl;the format of the file is like:
#find real difference: I_37595_G_A->G_A_G_G
while(<snp>){
	if(!/^#/){
		chomp;
		my $info=(split(/ /,$_))[-1];
		my($info_pos_ref_alt,$alles)=(split(/->/,$_))[0,1];
		my ($chr,$pos)=(split(/_/,$info_pos_ref_alt))[0,1];
		my ($m1,$m2,$w1,$w2)=(split(/_/,$alles))[0,1,2,3];
		print "snpfile split result: $chr,$pos,$m1,$m2,$w1,$w2\n";
		$m12aa=&snp2aa($chr,$pos,$m1,$gref); 
		$m22aa=&snp2aa($chr,$pos,$m2,$gref); 
		$w12aa=&snp2aa($chr,$pos,$w1,$gref); 
		$w22aa=&snp2aa($chr,$pos,$w2,$gref); 
		print "$chr,$pos,$m12aa,$m22aa,$w12aa,$w22aa,";
		if($m12aa eq $m22aa && $m22aa eq $w12aa &&$w22aa eq $w12aa){
			print "nonsense\n" ;
			}
		else{print "sense\n";}
		print $m12aa,"\n";
		}
		}
	close snp;
		
			