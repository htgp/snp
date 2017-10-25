#!/usr/bin/perl
sub getName2Gene{
	my ($info)=@_;
	if($info=~/;gene=/ && $info=~/;gene=/ ){
		$info=~/.*;Name=([^;]*);gene=([^;]*);.*/;
		if($1 ne ""  && $2 ne ""){
			return $1,$2;
			}
		else{return $infl,"err1";}
			}
	elsif($info=~/;Name=/){
			$info=~/.*;Name=([^;]*);.*/;
			if($1 ne ""){return $1,$1;}
			else{return $info,"err2";}
		}
		else{return $info,"err3";}
}	
1;