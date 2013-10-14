#
# Please cite:
# 1. Yulia A. Medvedeva, Marina V. Fridman, Nina J. Oparina, Dmitri B. Malko, Ekaterina O. Ermakova, Ivan V. Kulakovskiy, Vsevolod J. Makeev
# Evidence for transcriptional regulation by non-5' CpG islands in the human genome. BMC Genomics, 2008
# 2. Ina Y 
# New methods for estimating the numbers of synonymous and nonsynonymous substitutions. J Mol Evol, 1995
# 
# Three arguments required:
# 1. input file
# 2. output file
# 3. R (float-point number)
#
# Input file format (two codon sequences aligned, frame 1; sequence in one line required; '-' for gaps; 'N' recognized, but skipped; case doesn't matter): 
#<line 1> Some comments
#<line 2> >name1
#<line 3> atggggtga
#<line 4> *
#<line 5> >name2
#<line 6> atggg-ttt
#<line 7> *
#
use strict;
use MyModules::General;
use MyModules::Ina;
use MyModules::Identity;

open FPR, "<$ARGV[0]" or die "Cannot open input file ".$ARGV[0]."\n";
my $R=$ARGV[2];

#reading sequences
$_=<FPR>; $_=<FPR>;
my $seq1=<FPR>; 
$_=<FPR>; $_=<FPR>;
my $seq2=<FPR>; 
close FPR;
chomp ($seq1,$seq2);

#lowercase
$seq1=lc $seq1;
$seq2=lc $seq2;

print "Making appearance table...\n";
my @a=qw/a c t g/;
my @signif_codons=();
foreach my $i (@a)
{
	foreach my $j (@a)
	{
		foreach my $k (@a)
		{
			unless ($gct{$i.$j.$k} eq "Stop")
			{
				push @signif_codons, $i.$j.$k;
			}
		}
	}
}
my %pairs=(); my %s_table=(); my %n_table=(); my %sts_table=(); my %nts_table=(); my %stv_table=(); my %ntv_table=();
foreach my $c1 (@signif_codons)
{
	foreach my $c2 (@signif_codons)
	{
		$pairs{$c1.$c2}=0;
		$s_table{$c1.$c2}=0;
		$n_table{$c1.$c2}=0;
		$sts_table{$c1.$c2}=0;
		$nts_table{$c1.$c2}=0;
		$stv_table{$c1.$c2}=0;
		$ntv_table{$c1.$c2}=0;
	}
}
my $l=min(length($seq1),length($seq2));
for (my $j=0;$j+2<$l;$j+=3)
{
	my $c1=lc(substr($seq1,$j,3));
	my $c2=lc(substr($seq2,$j,3));
	if (defined $pairs{$c1.$c2})
	{
		$pairs{$c1.$c2}++;
	}
}

print "syntable...\n";
my %s=syntable($R);

print "s,n,sts,nts...\n";
foreach my $codon1 (@signif_codons)
{
	foreach my $codon2 (@signif_codons)
	{
		my $s1=$s{$codon1};
		my $s2=$s{$codon2};
		$s_table{$codon1.$codon2}=0.5*($s1+$s2);
		$n_table{$codon1.$codon2}=(3-0.5*($s1+$s2));
		my @c1=split "",$codon1;
		my @c2=split "",$codon2;
		my ($_sts,$_nts,$_stv,$_ntv)=(0,0,0,0);
		my $good_pathways=0;
		#Checking 6 pathways	
		#123
		my $pc1=$c2[0].$c1[1].$c1[2];
		my $pc2=$c2[0].$c2[1].$c1[2];
		if ( defined($s{$pc1}) and defined($s{$pc2}) )
		{
			$good_pathways++;
			$_sts+=( sts($codon1,$pc1) + sts($pc1,$pc2) + sts($pc2,$codon2) );
			$_nts+=( nts($codon1,$pc1) + nts($pc1,$pc2) + nts($pc2,$codon2) );
			$_stv+=( stv($codon1,$pc1) + stv($pc1,$pc2) + stv($pc2,$codon2) );
			$_ntv+=( ntv($codon1,$pc1) + ntv($pc1,$pc2) + ntv($pc2,$codon2) );
		}
		#132
		$pc1=$c2[0].$c1[1].$c1[2];
		$pc2=$c2[0].$c1[1].$c2[2];
		if ( defined($s{$pc1}) and defined($s{$pc2}) )
		{
			$good_pathways++;
			$_sts+=( sts($codon1,$pc1) + sts($pc1,$pc2) + sts($pc2,$codon2) );
			$_nts+=( nts($codon1,$pc1) + nts($pc1,$pc2) + nts($pc2,$codon2) );
			$_stv+=( stv($codon1,$pc1) + stv($pc1,$pc2) + stv($pc2,$codon2) );
			$_ntv+=( ntv($codon1,$pc1) + ntv($pc1,$pc2) + ntv($pc2,$codon2) );
		}
		#213
		$pc1=$c1[0].$c2[1].$c1[2];
		$pc2=$c2[0].$c2[1].$c1[2];
		if ( defined($s{$pc1}) and defined($s{$pc2}) )
		{
			$good_pathways++;
			$_sts+=( sts($codon1,$pc1) + sts($pc1,$pc2) + sts($pc2,$codon2) );
			$_nts+=( nts($codon1,$pc1) + nts($pc1,$pc2) + nts($pc2,$codon2) );
			$_stv+=( stv($codon1,$pc1) + stv($pc1,$pc2) + stv($pc2,$codon2) );
			$_ntv+=( ntv($codon1,$pc1) + ntv($pc1,$pc2) + ntv($pc2,$codon2) );
		}
		#231
		$pc1=$c1[0].$c2[1].$c1[2];
		$pc2=$c1[0].$c2[1].$c2[2];
		if ( defined($s{$pc1}) and defined($s{$pc2}) )
		{
			$good_pathways++;
			$_sts+=( sts($codon1,$pc1) + sts($pc1,$pc2) + sts($pc2,$codon2) );
			$_nts+=( nts($codon1,$pc1) + nts($pc1,$pc2) + nts($pc2,$codon2) );
			$_stv+=( stv($codon1,$pc1) + stv($pc1,$pc2) + stv($pc2,$codon2) );
			$_ntv+=( ntv($codon1,$pc1) + ntv($pc1,$pc2) + ntv($pc2,$codon2) );
		}
		#312
		$pc1=$c1[0].$c1[1].$c2[2];
		$pc2=$c2[0].$c1[1].$c2[2];
		if ( defined($s{$pc1}) and defined($s{$pc2}) )
		{
			$good_pathways++;
			$_sts+=( sts($codon1,$pc1) + sts($pc1,$pc2) + sts($pc2,$codon2) );
			$_nts+=( nts($codon1,$pc1) + nts($pc1,$pc2) + nts($pc2,$codon2) );
			$_stv+=( stv($codon1,$pc1) + stv($pc1,$pc2) + stv($pc2,$codon2) );
			$_ntv+=( ntv($codon1,$pc1) + ntv($pc1,$pc2) + ntv($pc2,$codon2) );
		}
		#321
		$pc1=$c1[0].$c1[1].$c2[2];
		$pc2=$c1[0].$c2[1].$c2[2];
		if ( defined($s{$pc1}) and defined($s{$pc2}) )
		{
			$good_pathways++;
			$_sts+=( sts($codon1,$pc1) + sts($pc1,$pc2) + sts($pc2,$codon2) );
			$_nts+=( nts($codon1,$pc1) + nts($pc1,$pc2) + nts($pc2,$codon2) );
			$_stv+=( stv($codon1,$pc1) + stv($pc1,$pc2) + stv($pc2,$codon2) );
			$_ntv+=( ntv($codon1,$pc1) + ntv($pc1,$pc2) + ntv($pc2,$codon2) );
		}
		if ($good_pathways)
		{
			$sts_table{$codon1.$codon2}=$_sts/$good_pathways;
			$nts_table{$codon1.$codon2}=$_nts/$good_pathways;
			$stv_table{$codon1.$codon2}=$_stv/$good_pathways;
			$ntv_table{$codon1.$codon2}=$_ntv/$good_pathways;
		}
	}
}
#
# S, N, STS, STV, NTS, NTV
#
print "S_arr...\n";
my @S_arr=(); my @N_arr=(); my @STS_arr=(); my @STV_arr=(); my @NTS_arr=(); my @NTV_arr=(); 
foreach my $key (keys %pairs)
{
	push @S_arr,$pairs{$key}*$s_table{$key};
	push @N_arr,$pairs{$key}*$n_table{$key};
	push @STS_arr,$pairs{$key}*$sts_table{$key};
	push @STV_arr,$pairs{$key}*$stv_table{$key};
	push @NTS_arr,$pairs{$key}*$nts_table{$key};
	push @NTV_arr,$pairs{$key}*$ntv_table{$key};
}
my ($S,$N,$STS,$STV,$NTS,$NTV)=(&sum_arr(@S_arr),&sum_arr(@N_arr),&sum_arr(@STS_arr),&sum_arr(@STV_arr),&sum_arr(@NTS_arr),&sum_arr(@NTV_arr)); 
#
# Pn Ps Qn Qs
#
print "Pn...";
my $Ps=$STS/$S;
my $Qs=$STV/$S;
my $Pn=$NTS/$N;
my $Qn=$NTV/$N;
#
# applicability
#
print "success...\n";
my $success="";
if ((1-2*$Ps-$Qs)<=0)
{
	$success="no1";
}
elsif ((1-2*$Qs)<=0)
{
	$success="no2";
}
elsif ((1-2*$Pn-$Qn)<=0)
{
	$success="no3";
}
elsif ((1-2*$Qn)<=0)
{
	$success="no4";
}
else
{
	$success="yes";
}

print "Writing results...\n";
open FPW, ">$ARGV[1]" or die "Cannot open output file ".$ARGV[1]."\n";
select FPW;
print "R\tsuccess\tIDn\tIDa\tdn\tds\tdn/ds\n";
print "$R\t$success\t";
print STDOUT "IDn...\n";
my $idenn=0;
my $lengthn=0;
for (my $j=0;$j<$l;$j++)
{
	my $a=lc(substr ($seq1,$j,1));
	my $b=lc(substr ($seq2,$j,1));
	if ($nuc{$a} and $nuc{$b})
	{
		$lengthn++;
		if ($a eq $b)
		{
			$idenn++;
		}
	}
}
my $idn=$idenn/$lengthn;
print STDOUT "IDa...\n";
my $idena=0;
my $lengtha=0;
for (my $j=0;$j+2<$l;$j+=3)
{
	my $codon1=lc(substr($seq1,$j,3)); 
	my $codon2=lc(substr($seq2,$j,3)); 
	if (defined($gct{$codon1}) and defined($gct{$codon2}))
	{
		$lengtha++;
		if ($gct{$codon1} eq $gct{$codon2})
		{
			$idena++;
		}
	}
}
my $ida=$idena/$lengtha;
print "$idn\t$ida\t";
my ($ds,$dn)=(0,0);
if ($success eq "yes")
{
	$ds=-0.5*log ( (1-2*$Ps-$Qs)*sqrt(1-2*$Qs) );
	$dn=-0.5*log ( (1-2*$Pn-$Qn)*sqrt(1-2*$Qn) );
	print "$dn\t$ds\t";
	if ($ds>0)
	{
		print $dn/$ds;
	}
	print "\n";
	
}
else
{
	print "\t\t\n";
}
print "STS $STS STV $STV NTS $NTS NTV $NTV S $S N $N LENGTHN $lengthn IDENN $idenn LENGTHA $lengtha IDENA $idena\n";
close FPW;

sub sum_arr (@)
{
	my @a=sort {$a <=> $b} @_;
	my $sum=0;
	for (my $j=0;$j<=$#a;$j++)
	{
		$sum+=$a[$j];
	}
	$sum;
}
