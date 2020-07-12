#!/usr/bin/perl
use List::Util qw( min max );
use Math::Round;

$src="/home/hdhiman/Data/work/hd/janssen/integration_sites/filter_regions/src2";
$results="/home/hdhiman/Data/work/hd/janssen/integration_sites/filter_regions/results2";

###################################################################################################
############################ Read source files ####################################################
###################################################################################################

open(VAR2, "$src/all_snp-2kb_mval_gt0.7.bedGraph") or die "Can't open VAR";
# print "Scanning high variability bins ...\n";
foreach(<VAR2>)
{
    chomp $_;
    @info_1=split(/\t/,$_);
    $var_start{$info_1[0]}.=$info_1[1]."~";
    $var_stop{$info_1[0]}.=$info_1[2]."~";
}

open(EXP, "$src/exp_gt20.bed") or die "Can't open EXP";
# print "#Scanning expression profile ...\n";
foreach(<EXP>)
{
    chomp $_;
    @info_2=split(/\t/,$_);
    $exp_start{$info_2[0]}.=$info_2[1]."~";
    $exp_stop{$info_2[0]}.=$info_2[2]."~";
}

open(REP, "$src/Tp5_11_segments_rep.bed") or die "Can't open REP";
# print "#Scanning repressed chromatin states ...\n";
foreach(<REP>)
{
    chomp $_;
    @info_3=split(/\t/,$_);
    $rep_start{$info_3[0]}.=$info_3[1]."~";
    $rep_stop{$info_3[0]}.=$info_3[2]."~";
}

open(QUI, "$src/Tp5_11_segments_qui.bed") or die "Can't open QUI";
# print "#Scanning quiescent state ...\n";
foreach(<QUI>)
{
    chomp $_;
    @info_4=split(/\t/,$_);
    $qui_start{$info_4[0]}.=$info_4[1]."~";
    $qui_stop{$info_4[0]}.=$info_4[2]."~";
}

open(ENH, "$src/Tp5_11_segments_enh.bed") or die "Can't open ENH";
# print "#Scanning enhancer states ...\n";
foreach(<ENH>)
{
    chomp $_;
    @info_5=split(/\t/,$_);
    $enh_start{$info_5[0]}.=$info_5[1]."~";
    $enh_stop{$info_5[0]}.=$info_5[2]."~";
}

open(PROM, "$src/Tp5_11_segments_prom.bed") or die "Can't open PROM";
# print "#Scanning promoter states ...\n";
foreach(<PROM>)
{
    chomp $_;
    @info_6=split(/\t/,$_);
    $prom_start{$info_6[0]}.=$info_6[1]."~";
    $prom_stop{$info_6[0]}.=$info_6[2]."~";
}

open(ACT, "$src/Tp5_11_segments_active.bed") or die "Can't open ACTIVE";
# print "#Scanning active transctription state ...\n";
foreach(<ACT>)
{
    chomp $_;
    @info_7=split(/\t/,$_);
    $act_start{$info_7[0]}.=$info_7[1]."~";
    $act_stop{$info_7[0]}.=$info_7[2]."~";
}

open(gff, "$src/GCF_003668045.1_CriGri-PICR_genomic.gff") or die "Can't open gff";
# print "#Scanning genes gff ...\n";

@gene_anno=<gff>;
foreach(@gene_anno)
{
    chomp $_;
    if($_ =~ /^NW/)
    {
	@info_8=split(/\t/,$_);
	$gene_nm=$info_8[8];

	if($info_8[2] =~ /^gene$/ || $info_8[2] =~ /pseudogene/)
	{
	    $chr_gene{$info_8[0]}.=$info_8[8]."~";
	    $gene_chr{$gene_nm}=$info_8[0];
	    $gene_start{$gene_nm}=$info_8[3];
	    $gene_stop{$gene_nm}=$info_8[4];
	    $gene_nm++;
	}

	if($info_8[2] =~ /exon/)
	{ 
	    $exon_start{$gene_nm}.=$info_8[3]."~";
	    $exon_stop{$gene_nm}.=$info_8[4]."~";
	}
    }
}


###################################################################################################
###################	Browse the genome 	###########################################################
###################################################################################################
$file_in=$ARGV[1].".txt"; 
$file_out=$ARGV[1]."_filtered.txt";

open(INDEX, "<$results/$file_in") or die "Can't open INDEX";
open(OUT, ">$results/$file_out") or die "Can't write OUT";

# print "#Reading the index ...\n";

# print "#Chr\tStart\tStop\tHighExpressionPeak\tHighVariabilityPeak\tQuiescent\tRepressed\tEnhancer\tPromoter\tActiveTranscription\tGeneRegion\n";

foreach(<INDEX>)
{
        chomp $_;
        @index=split(/\t/,$_);
        $start=$index[1];
        $stop=$index[2];

        $pos=($index[1]+$index[2])/2;
        @spl=split(/\./,$pos);
        $pos=$spl[0];
        
                $distance_exp=check($pos, $exp_start{$index[0]}, $exp_stop{$index[0]});
                $distance_var=check($pos, $var_start{$index[0]}, $var_stop{$index[0]});
                $distance_qui=check($pos, $qui_start{$index[0]}, $qui_stop{$index[0]});
                $distance_rep=check($pos, $rep_start{$index[0]}, $rep_stop{$index[0]});
                $distance_enh=check($pos, $enh_start{$index[0]}, $enh_stop{$index[0]});
                $distance_prom=check($pos, $prom_start{$index[0]}, $prom_stop{$index[0]});
                $distance_act=check($pos, $act_start{$index[0]}, $act_stop{$index[0]});

                $flag_gene=0;
                @check_gene=split(/~/,$chr_gene{$index[0]});

                for($c=1;$c<=@check_gene;$c++)
                {       
                        if($pos>=$gene_start{$check_gene[$c]} && $pos<=$gene_stop{$check_gene[$c]})
                        {
                                $flag_gene=1;
                                @exon_start_coord=split(/~/,$exon_start{$check_gene[$c]});
                                @exon_stop_coord=split(/~/,$exon_stop{$check_gene[$c]});
                                $flag_exon=0;
                                for($d=0; $d<@exon_start_coord;$d++)
                                {
                                        if($pos>=$exon_start_coord[$d] && $pos<=$exon_stop_coord[$d])
                                        {$flag_exon=1;last;}
                                }
                        last;
                        }
                }
                if($flag_gene==1 && $flag_exon==1){ $anno="exon";}
                if($flag_gene==1 && $flag_exon==0){ $anno="intron";}
                if($flag_gene==0){ $anno="intergenic";}
		print OUT "$index[0]\t$start\t$stop\t$distance_exp\t$distance_var\t$distance_qui\t$distance_rep\t$distance_enh\t$distance_prom\t$distance_act\t$anno\n";
                #print OUT "$index[0]\t$start\t$stop\t$distance_exp\t$distance_enh\t$distance_prom\t$distance_act\t$distance_var\t$distance_rep\t$distance_qui\t$anno\n";
}
close(INDEX);


###################################################################################################
##############################	Functions	#######################################################
###################################################################################################

sub check
{
    my @passed=@_; $min=0;

    @check_start=split(/~/, $passed[1]);
    @check_stop=split(/~/, $passed[2]);

    @diff_start=map { $_ - $passed[0] } @check_start;
    @diff_stop=map { $_ - $passed[0] } @check_stop;

    foreach(@diff_start){$_=abs($_);}
    foreach(@diff_stop){$_=abs($_);}

    $min_start = min(@diff_start);
    $min_stop = min(@diff_stop);

    if($min_start<= abs($min_stop)) {$min= $min_start;}
    else {$min=abs($min_stop);}

    for($i=0;$i<@check_start;$i++)
    {
	if($passed[0]>=$check_start[$i] && $passed[0]<=$check_stop[$i])
	{	$min=0;	}
    }
return($min);
}
