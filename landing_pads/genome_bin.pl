# genome_bin.pl -bin 2000 -index /data/borth/hdhiman/ref/picr_refseq/picr.fa.fai

$bin_size=$ARGV[1];
open(ip,$ARGV[3]) or die "Cant open ip";

open(out,">genome_binned-$bin_size\.bed");

foreach(<ip>)
{
	chomp $_;
	@info=split(/\t/,$_);
	
	for($start=1,$id=1;$start<=$info[1];$start+=$bin_size,$id++)
	{
		$stop=$start+$bin_size-1;
		if($stop<=$info[1])
		{print out "$info[0]\t$start\t$stop\t$id\n";}
		else
		{print out "$info[0]\t$start\t$info[1]\t$id\n";}
	}
}
