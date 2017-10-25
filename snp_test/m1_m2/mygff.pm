#!/usr/bin/perl
sub getName2Gene{
	my ($info)=@_;
	if($info=~/;gene=/){
		# $info=~/.*;Name=(t*[A-Z]*\(*[A-Za-z]{1,3}\)*[A-Z]*[0-9]{0,4}[A-Z]{0,1}-{0,1}[A-Za-z]{0,2});gene=(t*[A-Z]{0,2}\(*[A-Za-z]{1,8}\)*[0-9]{0,4}-*_*[A-Za-z0-9]{0,1});.*/;
		$info=~/.*;Name=([0-9]{0,3}_*-*t*.*\(*[A-Za-z]{0,8}\)*.*_*-*[0-9]{0-3}_*-*[A-Za-z]{0,8}_*-*[0-9]{0,3});gene=([0-9]{0,3}_*-*t*.*\(*[A-Za-z]{0,8}\)*.*_*-*[0-9]{0-3}_*-*[A-Za-z]{0,8}_*-*[0-9]{0,3});.*/;
		if($1 ne ""  && $2 ne ""){
			return $1,$2;
			}
		elsif($info=~/;Name=/){
			$info=~/.*;Name=([A-Za-z]{1,4}[0-9]{0,4}[A-Z]{0,1}-{0,1}[A-Z0-2]{0,2});Ontology.*/;
			if($1 ne ""){return $1,$1;}
			else{return $info,"err";}
		}
		else{}
	}
}	
1;