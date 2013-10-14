use strict;
use MyModules::General;
#
# Please cite:
# 1. Yulia A. Medvedeva, Marina V. Fridman, Nina J. Oparina, Dmitri B. Malko, Ekaterina O. Ermakova, Ivan V. Kulakovskiy, Vsevolod J. Makeev
# Evidence for transcriptional regulation by non-5' CpG islands in the human genome. BMC Genomics, 2008
# 2. Ina Y 
# New methods for estimating the numbers of synonymous and nonsynonymous substitutions. J Mol Evol, 1995
# 
# Two arguments required:
# 1. input file
# 2. output file
#
# Input file format (two codon sequences aligned, frame 1; sequence in one line required; '-' for gaps; 'N' recognized, but skipped; case doesn't matter): 
#<line 1> Some comments
#<line 2> >name1
#<line 3> atggggtga
#<line 4> *
#<line 5> >name2
#<line 6> atggg-ttt
#<line 7> *
open FPR, "<".$ARGV[0] or die "Cannot open input file ".$ARGV[0]."\n";
$_=<FPR>; $_=<FPR>;
my $seq1=<FPR>;
$_=<FPR>; $_=<FPR>;
my $seq2=<FPR>;
close FPR;
chomp ($seq1,$seq2);
open FPW, ">".$ARGV[1] or die "Cannot open output file ".$ARGV[1]."\n";
select FPW;
print "R=",&countR($seq1,$seq2),"\n";
close FPW;
#-------------------------------------------------------------------------------
sub countR ($$)
{
	my $R;
	my $seq1=$_[0];
	my $seq2=$_[1];
	my $l=&min(length($seq1),length($seq2));
	print "l, total = $l ";
	my $m=$l;
	my $P3=0;
	my $Q3=0;
	for (my $j=2;$j<$m;$j+=3)
	{
		my $codon1=substr($seq1,$j-2,3);
		my $codon2=substr($seq2,$j-2,3);
		my $c1=substr($seq1,$j,1);
		my $c2=substr($seq2,$j,1);
		if ($codon1=~/[-nN]+/ or $codon2=~/[-nN]+/)
		{
			$l-=3;
			next;
		}
		$P3+=&ts($c1,$c2);
		$Q3+=&tv($c1,$c2);
	}
	print "l, codons with deletions excluded = $l\n";
	print "$P3 $Q3\n";
	$P3/=($l/3);
	$Q3/=($l/3);
	print "$P3 $Q3\n";
	if ((1-2*$P3-$Q3)<=0)
	{
		print ("Cannot count R: ln(1-2P-Q) undefined\n");
		return;
	}
	if((1-2*$Q3)<=0)
	{
		print ("Cannot count R: ln(1-2Q) undefined\n");
		return;
	}
	if (log(1-2*$Q3)==0)
	{
		print ("Cannot count R: ln(1-2Q)=0\n");
		return;
	}
	$R=(2*log((1-2*$P3)-$Q3)/log(1-2*$Q3))-1;
	$R;
}
