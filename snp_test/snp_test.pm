#!/usr/bin/perl
use Bio::Seq;
use Bio::SeqIO;

my $upstream=100;

sub snp2aa{
	my($snp_chr,$snp_pos,$snps_ref,$gff_ref)=@_;
	print "in snp2aa, paras are : $snp_chr,$snp_pos,$snps_ref,$gff_ref.snps are @{$snps_ref}.\n";
	my($aa_ref, $gnm,$abpos,$type)=&snp_location_and_translation($snp_chr,$snp_pos,$snps_ref,$gff_ref);
	# print "gnm in snp2aa: ",$gnm,"\n";
	return $aa_ref, $gnm,$abpos,$type;
	}
	
sub snp2aa_multi{
	$snp_chr=$snp_pos=$snps_ref=$gff_ref=$aa_ref="";
	my($snp_chr,$snp_pos,$snps_ref,$gff_ref)=@_;
	# print "in snp2aa_multi, paras are : $snp_chr,$snp_pos,$snps_ref,$gff_ref.snps are @{$snps_ref}.\n";
	my $aa_ref=&snp_location_and_translation($snp_chr,$snp_pos,$snps_ref,$gff_ref);
	# print "in snp2aa_multi,return arr ref if: $aa_ref, unpacked is: @{$aa_ref}.\n";
	return $aa_ref;
	}
sub nt2aa{
{   
    my($codon) = @_;   
    # print "in nt2aa, para condon accepted is: $codon.\n";
    $codon = uc $codon;#uc=uppercase;lc=lowercase  
    my(%genetic_code) = (
    'TCA' => 'S',    # Serine   
    'TCC' => 'S',    # Serine   
    'TCG' => 'S',    # Serine   
    'TCT' => 'S',    # Serine   
    'TTC' => 'F',    # Phenylalanine   
    'TTT' => 'F',    # Phenylalanine   
    'TTA' => 'L',    # Leucine   
    'TTG' => 'L',    # Leucine   
    'TAC' => 'Y',    # Tyrosine    
    'TAT' => 'Y',    # Tyrosine   
    'TAA' => '*',    # Stop   
    'TAG' => '*',    # Stop   
    'TGC' => 'C',    # Cysteine   
    'TGT' => 'C',    # Cysteine   
    'TGA' => '*',    # Stop   
    'TGG' => 'W',    # Tryptophan   
    'CTA' => 'L',    # Leucine   
    'CTC' => 'L',    # Leucine   
    'CTG' => 'L',    # Leucine   
    'CTT' => 'L',    # Leucine   
    'CCA' => 'P',    # Proline   
    'CCC' => 'P',    # Proline   
    'CCG' => 'P',    # Proline   
    'CCT' => 'P',    # Proline   
    'CAC' => 'H',    # Histidine   
    'CAT' => 'H',    # Histidine   
    'CAA' => 'Q',    # Glutamine   
    'CAG' => 'Q',    # Glutamine   
    'CGA' => 'R',    # Arginine   
    'CGC' => 'R',    # Arginine   
    'CGG' => 'R',    # Arginine   
    'CGT' => 'R',    # Arginine   
    'ATA' => 'I',    # Isoleucine   
    'ATC' => 'I',    # Isoleucine   
    'ATT' => 'I',    # Isoleucine   
    'ATG' => 'M',    # Methionine   
    'ACA' => 'T',    # Threonine   
    'ACC' => 'T',    # Threonine   
    'ACG' => 'T',    # Threonine   
    'ACT' => 'T',    # Threonine   
    'AAC' => 'N',    # Asparagine   
    'AAT' => 'N',    # Asparagine   
    'AAA' => 'K',    # Lysine   
    'AAG' => 'K',    # Lysine   
    'AGC' => 'S',    # Serine   
    'AGT' => 'S',    # Serine   
    'AGA' => 'R',    # Arginine   
    'AGG' => 'R',    # Arginine   
    'GTA' => 'V',    # Valine   
    'GTC' => 'V',    # Valine   
    'GTG' => 'V',    # Valine   
    'GTT' => 'V',    # Valine   
    'GCA' => 'A',    # Alanine   
    'GCC' => 'A',    # Alanine   
    'GCG' => 'A',    # Alanine   
    'GCT' => 'A',    # Alanine       
    'GAC' => 'D',    # Aspartic Acid   
    'GAT' => 'D',    # Aspartic Acid   
    'GAA' => 'E',    # Glutamic Acid   
    'GAG' => 'E',    # Glutamic Acid   
    'GGA' => 'G',    # Glycine   
    'GGC' => 'G',    # Glycine   
    'GGG' => 'G',    # Glycine   
    'GGT' => 'G',    # Glycine   
    );   
   
    if(exists $genetic_code{$codon})   
    {   
       # print $genetic_code{$codon}," is the aa.\n";
       return $genetic_code{$codon};   
    }  
    else  
    {   
		return "?";
            # print STDERR "Bad codon \"$codon\"!!\n";   
            # exit;   
    }   
}  
	}
sub get_tripleton{
	#return a tripleton;
	my($snp_p,$gl,$gr)=@_;
	# print "in get_tripleton,paras are: $snp_p,$gl,$gr.\n";
	my $ogl,$ogr;
	my $modul=($snp_p-$gl+1)%3;
	# print "modul:$modul.\n";
	if($modul==0){$ogl=$snp_p-2;$ogr=$snp_p;}
	elsif($modul==2){$ogl=$snp_p-1;$ogr=$snp_p+1;}
	else{$ogl=$snp_p,$ogr=$snp_p+2;}
	return $modul,$ogl,$ogr;
	}
sub verify_snp{
	#is the snp provided identical to the fasta file at same position?
	}
sub process_cds{
	my($snp_pos,$pckey,$snp_ref,$g_hash)=@_;
	# print "process_cds: paras are $snp_pos,$pckey,$snp,$g_hash.\n";
	my @snp_arr=@{$snp_ref};
	# print "in process_cds, snp arr are: ",@snp_arr,@outaas,"\n";
	my $uniq_snp_ref=simplify_snp($snp_ref);
	my %uniq_snp=%{$uniq_snp_ref};
	my %gff=%{$g_hash};
	my $aref=&parse_line($pckey);
	my @tmparr=@{$aref};
	my $gene_left=$tmparr[1];my $gene_right=$tmparr[2];my $direct=$tmparr[4];
	my($mod,$ngl,$ngr)=&get_tripleton($snp_pos,$gene_left,$gene_right);
	# print "module,tripleton_left,tripleton_right are: $mod,$ngl,$ngr.\n";
	# print "in process_cds: paras for extract_sequence are $tmparr[0],$ngl,$ngr;\n";
	my $tripleton=&extract_sequence($tmparr[0],$ngl,$ngr,,);
	# print $tripleton,"\n";
	for my $key(keys %uniq_snp){
	if($mod==0){$tripleton=&replace($tripleton,3,$key);}#print "kicked in tripleton: $tripleton.\n";
	if($mod==1){$tripleton=&replace($tripleton,1,$key);}#print "kicked in tripleton: $tripleton.\n";
	if($mod==2){$tripleton=&replace($tripleton,2,$key);}#print "kicked in tripleton: $tripleton.\n";
	if($tmparr[4] eq "-"){$tripleton=&rc_short_seq($tripleton);}
	$uniq_snp{$key}=&nt2aa($tripleton);
	# print "\%uniq_snp: $key, $uniq_snp{$key}.\n";
		}
	@outaas=map{$uniq_snp{$_}} @snp_arr;
	# print "mapping... @snp_arr, @outaas.\n";
	return \@outaas;
	}
sub replace{
		my($ori,$site,$kicker)=@_;
		my @parse=split('',$ori);
		$parse[$site-1]=$kicker;
		return join('',@parse);
		}
sub process_RNAs{
	print "RNA change";
}
sub rc_short_seq{
	my($ss)=(@_);
	my $rc_out="";
	my @ss_arr=split('',$ss);
	# print "split tripleton in rc_short_seq, ", join(',',@ss_arr),"\n";
	for(my $i=scalar(@ss_arr);$i>0;$i--){
		if(uc $ss_arr[$i-1] eq "A"){$rc_out .=  "T";}
		if(uc $ss_arr[$i-1] eq "T"){$rc_out .=  "A";}
		if(uc $ss_arr[$i-1] eq "C"){$rc_out .=  "G";}
		if(uc $ss_arr[$i-1] eq "G"){$rc_out .=  "C";}
}
		return $rc_out;
}
sub extract_sequence{
	my($chr,$start,$end,$direction,$info)=@_;
	# print $chr,$start,$end,$direction,$info,"\n";
	my $outseq;
	my $seq;
	my $io=Bio::SeqIO->new(-file=>'sc61.fa',-format=>'Fasta');
	while($seq=$io->next_seq){
		my $gene_name=$seq->display_id();
		# print $gene_name,"\n";
		if($gene_name eq $chr){
			$outseq=$seq->subseq($start,$end);
			# print "$outseq\n";
		}
		}
				return $outseq;
	}
sub random_seq{
	my($len)=@_;
	my $outseq="";
	while(my $i<$len){
		my @nucleotides=qw/A T G C/;  
		my $newbase=$nucleotides[rand @nucleotides];  
		$outseq .= $newbase;
		++$i;
		} 
		return $outseq
}		
sub hash_gff{
	my ($file)=@_;
	my %gff;
	open fl,"<",$file or die open "open gff file failed...";
	while(<fl>){
		chomp;
		if(!/^#/){
			# {'chr'}=$arr[0];
			# {'start'}=$arr[1];
			# {'end'}=$arr[2];
			# {'type'}=$arr[3];
			# {'direction'}=$arr[4];
			# {'info'}=$arr[5];
			$gff{$_}=1;
		}
	}
	close f1;
	# print "chr:$arr[0]\nstart:$arr[1]\nend:$arr[2]\ntype:$arr[3]\ndirection:$arr[4]\ninfo:$arr[5]\n";
	return %gff;
	}
sub snp_location_and_translation{
	my($chr,$pos,$snp_ref,$href)=@_;
	# print  "in SLAT, para snps are: ",@{$snp_ref},"\n";
	my @total_out=();
	# print "in process snp_location_and_translation,paras are: $chr,$pos,$snp_ref,$href. snp_ref rec are: @{$snp_ref}.\n";
	for my $key(keys %$href){
		my $aref=&parse_line($key);
		my @tmparr=@{$aref};
		# print "Doing SLAT, parsed keys are: ",$tmparr[0]," ",$tmparr[1]," ",$tmparr[2]," ",$tmparr[3]," ",$tmparr[4],"...\n";
		if(($tmparr[0] eq $chr && $tmparr[1]-$upstream<=$pos && $pos<=$tmparr[1] && $temarr[4] eq "+" && $tmparr[3] eq "CDS") || ($tmparr[0] eq $chr && $tmparr[2]<=$pos && $pos<=$tmparr[2]+$upstream && $temarr[4] eq "-" && $tmparr[3] eq "CDS")) {#up;
			my $gnm=&get_gene_name($tmparr[5]);
			my $outaa_ref="other";
			my $up_pos;
			if($temarr[4] eq "+"){$up_pos=$pos-$tmparr[1]}
			else{$up_pos=$tmparr[2]-$pos}
			my $type="upstream";
			# print "$chr-$pos:-100 upstream of coding region which may have an impact at site $pos of $gnm.\n";
			push (@total_out,  ($outaa_ref,$gnm,$up_pos,$type));
			}#up
		if($tmparr[0] eq $chr && $tmparr[1]<=$pos && $pos<=$tmparr[2]){#hit
			# print "Doing SLAT, parsed keys are: ",$tmparr[0]," ",$tmparr[1]," ",$tmparr[2]," ",$tmparr[3]," ",$tmparr[4],"...\n";
			if($tmparr[3] eq 'CDS'){
				# print "CDS region.\n";
				# print "in SLAT#, para snp_ref: ",@{$snp_ref},"\n";
				my $outaa_ref=&process_cds($pos,$key,$snp_ref,$href);
				# print "in SLAT*: ", @{$outaa_ref}, " are the outaas\n";
				my $gnm=&get_gene_name($tmparr[5]);
				# print "gnm in SLAT is: ",$gnm,"\n";
				my $type=$tmparr[3];
				my $ab_aa_pos;
				if($tmparr[4] eq "+"){
					$ab_aa_pos=&aa_abs_pos($tmparr[1],$pos);
					}
				else{
					$ab_aa_pos=&aa_abs_pos($pos,$tmparr[2]);
					}
				push (@total_out,  ($outaa_ref,$gnm,$ab_aa_pos,$type));
					}
			elsif($tmparr[3] =~ /RNA/){
				# &process_RNAs($key);
				my $type=$tmparr[3];
				my $outaa_ref="";
				my $gnm=&get_gene_name($tmparr[5]);
				my $ab_aa_pos;
				if($tmparr[4] eq "+"){
					$ab_aa_pos=$pos-$tmparr[1];
					}
				else{
					$ab_aa_pos=$tmparr[2]-$pos;
					}
				push (@total_out,( $outaa_ref,$gnm,$ab_aa_pos,$type));		
				}
			# elsif($tmparr[3] =~/chromosome/){next;}
			else{
				my ($outaa_ref,$gnm,$ab_aa_pos,$type)=("other","other",$pos,$tmparr[3]);;
				push (@total_out,  ( $outaa_ref,$gnm,$ab_aa_pos,$type));
				}

			}#hit
			}
		return \@total_out;
			}
sub nc_location{
	my($chr,$pos,$snp_ref,$href)=@_;
	my @total_out;
	my @snps=@{$snp_ref};
	for my $key(keys %$href){
		my $aref=&parse_line($key);
		my @tmparr=@{$aref};
		if(($tmparr[0] eq $chr && $tmparr[1]-$upstream<=$pos && $pos<=$tmparr[1] && $temarr[4] eq "+" && $tmparr[3] eq "CDS") || ($tmparr[0] eq $chr && $tmparr[2]<=$pos && $pos<=$tmparr[2]+$upstream && $temarr[4] eq "-" && $tmparr[3] eq "CDS")) {
			my $gnm=&get_gene_name($tmparr[5]);
			my $outaa_ref="other";
			my $up_pos;
			if($temarr[4] eq "+"){$up_pos=$pos-$tmparr[1]}
			else{$up_pos=$tmparr[2]-$pos}
			my $type="upstream";
			# print "$chr-$pos:-100 upstream of coding region which may have an impact at site $pos of $gnm.\n";
			push (@total_out,  ($outaa_ref,$gnm,$up_pos,$type));
			}
		if($tmparr[0] eq $chr && $tmparr[1]<=$pos && $pos<=$tmparr[2]){
			# print "Doing SLAT, parsed keys are: ",$tmparr[0]," ",$tmparr[1]," ",$tmparr[2]," ",$tmparr[3]," ",$tmparr[4],"...\n";
			if($tmparr[3] eq 'CDS'){
				my $outaa_ref="other";
				my $gnm=&get_gene_name($tmparr[5]);
				# print "gnm in SLAT is: ",$gnm,"\n";
				my $type=$tmparr[3];
				my $ab_aa_pos;
				if($tmparr[4] eq "+"){
					$ab_aa_pos=&aa_abs_pos($tmparr[1],$pos);
					}
				else{
					$ab_aa_pos=&aa_abs_pos($pos,$tmparr[2]);
					}
				push (@total_out,  ( $outaa_ref,$gnm,$ab_aa_pos,$type));
					}
			elsif($tmparr[3] =~ /RNA/){
				# &process_RNAs($key);
				my $type=$tmparr[3];
				my $outaa_ref="other";
				my $gnm=&get_gene_name($tmparr[5]);
				my $ab_aa_pos;
				if($tmparr[4] eq "+"){
					$ab_aa_pos=$pos-$tmparr[1];
					}
				else{
					$ab_aa_pos=$tmparr[2]-$pos;
					}
				push (@total_out,( $outaa_ref,$gnm,$ab_aa_pos,$type));		
				}
			else{
				my ($outaa_ref,$gnm,$ab_aa_pos,$type)=("other","other",$pos,$tmparr[3]);;
				push (@total_out,  ( $outaa_ref,$gnm,$ab_aa_pos,$type));
				}

			}
			}
		return \@total_out;
			}
sub aa_abs_pos{
	my ($p1,$p2)=@_;
	my $out=int(($p2-$p1)/3+0.99);
	# print "aa_abs_pos, paras accepted are: $p1,$p2, result is: $out.\n";
	return $out;
}
sub get_gene_name{
	my ($info)=@_;
	if($info!~/;gene=.*;Alias=.*;.*/){#tRNA;ID=tK(UUU)O;Name=tK(UUU)O;Ontology_ter;
		$info=~/.*;Name=(.*);Ontology_term=.*/;
		# print "tRNA:$1. info now is $info\n";
		return $1;
		}
	else{
		$info=~/.*;Name=(.*);gene=.*;Alias=.*;Ontology_term=.*/;#cds,snRNA,rRNA;
		# print "cds or rRNA:$1; . info now is $info\n";#ID=NME1;Name=NME1;gene=NME1;Alias=NME1,RRP2;Ontology_term=
		return $1;
		}
}	
sub parse_line{
	(my $line)=@_;
	my ($chr,$db,$type,$gl,$gr,$u1,$dir,$u2,$info)=split('\t',$line);
	my @arr=($chr,$gl,$gr,$type,$dir,$info);
	return \@arr;
	}
sub is_snp{
	my @nts=@_;
	# print join('_',@nts),"\n";
	my $flag=1;
	for(@nts){
		# print "$_,$flag\n";
		if(length($_) != 1){
			$flag=0;
			}
		if($flag eq 0){
			return 0;
			last;
			}
		}
		if($flag eq 1){
			return 1;
			}
		}
sub simplify_snp{
		my($snp_ref)=@_;
		my @snps=@{$snp_ref};
		my %uniq_snp;
		for my $key(@snps){
			$uniq_snp{$key}=1;
			}
		return \%uniq_snp;
		}

			
	
1;