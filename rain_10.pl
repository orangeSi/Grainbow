#!/usr/bin/perl -w
use PerlIO::gzip;

if(@ARGV lt 3){
	print "
	perl $0 <ref.fa> <libxx.PE.gz list;insert size that's from big to short is best for view> <yes or no;whether filter some PE>\n
		detail: 
			color:
				red -> gap region;
				black -> genome sequence
			filter: will filter some PE if read1 of one pair is overlaped with read1 of another pair reads.
	author: sikaiwei\@genome.cn

";
	exit;
}
my $filterTAG=$ARGV[2];
my @lists;
for(`cat $ARGV[1]`){
	chomp;
	push @lists,$_;
}
my @colors=("green","blue","yellow");
my $ref="$ARGV[0]";
my %reflen;
my %scaffolds;
`rm *.svg 2>/dev/null`;
`fastalength $ref >$ref.length`;
open LEN,"$ref.length" or die "$!";
while(<LEN>){
	$_=~ /^(\S+)\s+(\S+)$/;
	$reflen{$2}=$1;
}
close LEN;

open IN,"$ref" or die "$!";
$/=">";
<IN>;
while(<IN>){
	chomp;
	my @arr=split(/\n/,$_,2);
	$arr[1]=~ s/\s+//g;
	my @bps=split(//,$arr[1]);
	my @tmp=$arr[1]=~ /N+/g;

#	print "$arr[0] has ",scalar(@tmp),"gaps\n";
	my $flag=0;
	my $i=0;
	my $gap=0;
	for(@bps){
		$i++;
		if($_ eq "N"){


			if(!$gap){$flag++;$scaffolds{$arr[0]}{'gaps'}{$flag}{'s'}=$i;#print "$arr[0] gap$flag start $i\n"
			}
			$gap=1;
			if($bps[$i] ne "N"){$scaffolds{$arr[0]}{'gaps'}{$flag}{'e'}=$i;#print "$arr[0] gap$flag end $i\n"
			}
		}else{
			$gap=0;
		}

	}
	my $numGaps=scalar(keys %{$scaffolds{$arr[0]}{'gaps'}});
#	print "$arr[0] has $numGaps gaps\n";
#	$scaffolds{$arr[0]}=$arr[1];

}
close IN;
$/="\n";


my $flag2;
## fetch  short reads  #############################################################
foreach my $k(@lists){
#	my $color=shift @colors;
	my %lib;       
#`gzip -dc $k|awk -v line=1 '{j=line%2;if(j==1){print length(\$2),\$8,\$9};line=line+1}' |gzip >$k.tmp.gz`;
	my $tmpfile="$k.tmp.gz";
	if( -f $tmpfile){$tmpfile=1}else{$tmpfile=0}
	if(!$tmpfile){
		open OUT,"|sort -k 1,1 -k 2n,2n -S 5000M |gzip >$k.tmp.gz" or die "$!";
		open IN,"<:gzip(autopop)","$k" or die "$!";

		while(<IN>){
			chomp;
			my @arr=split(/\s+/,$_);
			my $len=length $arr[1];
			my $scf=$arr[7];
			$scf =~ s/\|/_/g;
			my $pos_start1=$arr[8];
			my $line2=<IN>;
			@arr=split(/\s+/,$line2);
			my $pos_start2=$arr[8];
			if($pos_start1 > $pos_start2){my $end=$pos_start1;my $start=$pos_start2;$pos_start1=$start;$pos_start2=$end}
			if(!$tmpfile){
				print OUT "$scf\t$pos_start1\t$pos_start2\t$len\n";
			}


		}
		close IN;
		close OUT;
	}
	my $flag1;
	my $before1;
	my $beforeid;
#	my $before2;
#	my $scf;
	open IN,"gzip -dc $k.tmp.gz|" or die "$!";
	while(<IN>){
		chomp;
		$flag1+=1;
#		print "flag is $flag1\n";
		my @arr=split(/\s+/,$_);
		my $len=$arr[3];
		my $scf=$arr[0];
		$scf =~ s/\|/_/g;
		my $pos_start1=$arr[1];
		my $pos_start2=$arr[2];

### filter some PE 
		if($flag1!=1 && $beforeid eq $scf && $filterTAG eq "yes"){

			my $dis1=$pos_start1 - $before1;
#my $dis2=abs($pos_start2 - $before2);
			if($dis1 < $len*0.2){
				next;
			}	

		}
#print "$scf\t$pos_start1\n";
		my $insert=abs($pos_start2 - $pos_start1) + $len;
		my $color;
		if($insert < 800 ){
			$color="red"
		}elsif($insert < 1500 ){
			$color="green";
		}elsif($insert < 3000){
			$color="blue";
		}elsif($insert < 7000 ){
			$color="gray";
		}elsif($insert < 10000){
			$color="black";
		}else{
			## gt 10k 
			next;
			}
		$color=$color ? $color:"yellow";
#		print "flag1 is $flag1\n";	
		$lib{$scf}{$flag1}{'c'}=$color;## clolor
			$lib{$scf}{$flag1}{'l'}=$len;  ## reads length
			$lib{$scf}{$flag1}{'start1'}=$pos_start1; ## PE read1 start positon in ref
			$lib{$scf}{$flag1}{'s2'}=$pos_start2;  ## PE read2 start positon in ref
			$before1=$pos_start1;
		$beforeid=$scf;
#		$before2=$pos_start2;
	}
	close IN;

	foreach my $scf(keys %lib){
		open SVG,">>$scf.svg" or die "$!";
		my $scflen=$reflen{$scf};
		$radio=0.05;                        ## $radio piex per bp
			my $width1=$radio*$scflen + 100; 
		my $height1=$radio*1000*2;
## scffold sequence line background
		my $p1x=50;
		my $p1y=$height1*2/3;
		my $p2x=$p1x;
		my $p2y=$p1y+5;;
		my $p3x=$p1x+$radio*$scflen;
		my $p3y=$p2y;
		my $p4x=$p3x;
		my $p4y=$p1y;
		my $svg;
		if(-z "$scf.svg" ){
			my $tmp="<svg width=\"$width1\" height=\"$height1\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">\n";
### draw genome line
			$tmp.="<path d=\"M$p1x $p1y L$p2x $p2y L$p3x $p3y L$p4x $p4y Z\" />\n";
			print SVG $tmp;


			foreach my $kk(keys %{$scaffolds{$scf}{'gaps'}}){

				my $start=$scaffolds{$scf}{'gaps'}{$kk}{'s'};
				my $end=$scaffolds{$scf}{'gaps'}{$kk}{'e'};
				my $pp1x=$p1x+$start*$radio;	
				my $pp1y=$p1y+6;
				my $pp2x=$pp1x;
				my $pp2y=$p2y+6;
				my $pp3x=$p1x+$end*$radio;
				my $pp3y=$pp2y;
				my $pp4x=$pp3x;
				my $pp4y=$pp1y;
### draw gap region
				print SVG "<path d=\"M$pp1x $pp1y L$pp2x $pp2y L$pp3x $pp3y L$pp4x $pp4y Z\" style=\"fill:red;\" />\n";
			}



		}

### draw PE half circos		
		foreach my $flag(sort {$lib{$scf}{$a}{'start1'}<=>$lib{$scf}{$b}{'start1'}} keys %{$lib{$scf}}){
			my $s1=$lib{$scf}{$flag}{'start1'};
			my $length=$lib{$scf}{$flag}{'l'};

			my $width=$length * $radio /2 ;# for not that thick
				my $s2=$lib{$scf}{$flag}{'s2'};
			my $sx= $p1x + $s1*$radio+$width/2;
			my $sy= $p1y;
			my $ex= $p1x + $s2*$radio + $width/2;
			my $ey= $p1y;
			my $tag=($s1 < $s2)? 1:0;
			my $d=abs($s1 - $s2) +1 - $length;
			my $r=$d/2*$radio;
			my $r2=$r/1.2;
#		        if(!($svg && $scf && $s1 && $s2 && $sx && $sy && $r && $r2 && $tag && $ex && $ey && $lib{$scf}{$flag}{'c'} && $width )){
#				print "svg is $svg \n scf is $scf \n s1 $s1 \n s2 $s2 \n sx $sx \n xy $sy \n r $r \n r2 $r2 \n tag $tag \n ex $ex \n ey $ey \n col $lib{$scf}{$flag}{'c'} \n width $width\n";

#				}
			$svg.="<g style=\"fill: none\"><title>$scf;$s1;$s2</title>\n<path d=\"M$sx $sy A$r $r2 0 0 $tag  $ex $ey\"  style=\"stroke:$lib{$scf}{$flag}{'c'};stroke-width:$width\" />\n</g>\n";
		}
		print SVG $svg;
		close SVG;
	}


}
### while end



### 
foreach my $k (glob("*.svg")){
	`echo "</svg>" >>$k`;

}





