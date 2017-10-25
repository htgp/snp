#!/usr/bin/perl
BEGIN{
	unshift(@INC,".");
	}
use snp_test;
if($ARGV[0] eq ""){print "Usage: perl snp_test.pl snp_find_output_file > output_file\n";exit;}
$fa="sc61.fa";
$gff="sc61_rome.gff";
$snp=$ARGV[0];
%gff=hash_gff($gff);
$gref=\%gff;

&loader;


sub loader{
		open snp,"<","$snp" or die "open snp file failed...\n";
		#this snp file is like:I_106385_C_T->T/T:C/T;
 while(<snp>){
	if(!/^#/){#if_1
		print;
		chomp;
		
		my ($key,$comp)=split('->',$_);
		$key=(split(' ',$key))[-1];
		my ($chr,$pos,$ref,$alt)=split('_',$key);
		my @snps=($m1a,$m1b,$m2a,$m2b)=split('_',$comp);
		my $sr=\@snps;
		# print "snps: @snps; its ref: $sr; after recovery: @{$sr}\n";
		if(&is_snp(@snps)==1){#if_snp;
			my $outaas=&snp2aa_multi($chr,$pos,$sr,$gref); 
			# print "returned arr ref: $outaas\n";
			my @list=@{$outaas};
			# print "outaas: $outaas. unpack it: @list\n";
			while(scalar @list>0){
				my($aas,$gnm,$aa_pos,$type)=(shift @list,shift @list,shift @list,shift @list);
				# print "**$aas,$gnm,$aa_pos,$type\n";
				if($type eq "CDS"){
					($m1a2aa,$m1b2aa,$m2a2aa,$m2b2aa)=@{$aas};
					# print "$m1a2aa,$m1b2aa,$m2a2aa,$m2b2aa\n";
					if($m1a2aa eq $m1b2aa && $m2a2aa eq $m1b2aa && $m2a2aa eq $m2b2aa){
						print "$chr:$pos->nt($m1a/$m1b:$m2a/$m2b),aa($m1a2aa/$m1b2aa:$m2a2aa/$m2b2aa) at site $aa_pos of protein $gnm:Nosense.\n";
					}#if_nosense;
					else{
						print "$chr:$pos->nt($m1a/$m1b:$m2a/$m2b),aa($m1a2aa/$m1b2aa:$m2a2aa/$m2b2aa) at site $aa_pos of protein $gnm:Sense.\n";
						}#else_sense;
					}#if_cds;
				elsif($type=~/RNA/){print "$chr:$pos->nt($m1a/$m1b:$m2a/$m2b)at site $aa_pos of $type $gnm.\n";}
				elsif($type eq "upstream"){print "$chr:$pos->nt($m1a/$m1b:$m2a/$m2b)at upstream $aa_pos of protein $gnm.\n";}
				else{print "$chr:$pos->nt($m1a/$m1b:$m2a/$m2b)at site $aa_pos of $type that will be ignored.\n";}
				$aas=$gnm=$aa_pos=$type="";
				}#while;
				$outaas="";undef @list;
				}#if_snp;
		else{#indel;
			my $outaas=&nc_location($chr,$pos,\@snps,$gref);
			print "returned arr ref: $outaas\n";
			my @list=@{$outaas};
			# print "outaas: $outaas. unpack it: @list\n";
			while(scalar @list>0){
				my($aas,$gnm,$aa_pos,$type)=(shift @list,shift @list,shift @list,shift @list);
				# print "**$aas,$gnm,$aa_pos,$type\n";
				if($type eq "CDS"){
						print "$chr:$pos->nt($m1a/$m1b:$m2a/$m2b),indel induced frame shift at site $aa_pos of protein $gnm:Sense.\n";
					}#if_cds;
				elsif($type=~/RNA/){print "$chr:$pos->nt($m1a/$m1b:$m2a/$m2b),indel variation at site $aa_pos of $type $gnm.\n";}
				elsif($type eq "upstream"){print "$chr:$pos->nt($m1a/$m1b:$m2a/$m2b),indel at upstream $aa_pos of protein $gnm.\n";}
				else{print "$chr:$pos->nt($m1a/$m1b:$m2a/$m2b),indel variation at site $aa_pos of $type that will be ignored.\n";}
				$aas=$gnm=$aa_pos=$type="";
				}
				undef @list;		
			}#else_indel;
				}#if_1;
			$key=$comp=	$chr=$pos=$change=$m1a=$m1b=$m2a=$m2b="";
			undef @snps;
				}#while
		close snp;
				}#sub;
			