package MyModules::Ina;
#
# Please cite:
# 1. Yulia A. Medvedeva, Marina V. Fridman, Nina J. Oparina, Dmitri B. Malko, Ekaterina O. Ermakova, Ivan V. Kulakovskiy, Vsevolod J. Makeev
# Evidence for transcriptional regulation by non-5' CpG islands in the human genome. BMC Genomics, 2008
# 2. Ina Y 
# New methods for estimating the numbers of synonymous and nonsynonymous substitutions. J Mol Evol, 1995
# 
use strict;
use vars qw(@ISA @EXPORT $VERSION);
use Exporter;
$VERSION=1.00;
@ISA=qw(Exporter);
@EXPORT=qw(syntable nts ntv sts stv);
use MyModules::General;

sub syntable ($) 
{
        my $R=$_[0];
        my %s=();
        
        $s{"ttt"}=$R/(1+$R);
        $s{"ttc"}=$R/(1+$R);
        $s{"tta"}=2*$R/(1+$R);
        $s{"ttg"}=2*$R/(1+$R);
        
        $s{"ctt"}=1;
        $s{"ctc"}=1;
        $s{"cta"}=(1+2*$R)/(1+$R);
        $s{"ctg"}=(1+2*$R)/(1+$R);
        
        $s{"att"}=(0.5+$R)/(1+$R);
        $s{"atc"}=(0.5+$R)/(1+$R);
        $s{"ata"}=1/(1+$R);
        $s{"atg"}=0;
        
        $s{"gtt"}=1;
        $s{"gtc"}=1;
        $s{"gta"}=1;
        $s{"gtg"}=1;
        
        $s{"tct"}=1;
        $s{"tcc"}=1;
        $s{"tca"}=1;
        $s{"tcg"}=1;
        
        $s{"cct"}=1;
        $s{"ccc"}=1;
        $s{"cca"}=1;
        $s{"ccg"}=1;
        
        $s{"act"}=1;
        $s{"acc"}=1;
        $s{"aca"}=1;
        $s{"acg"}=1;
        
        $s{"gct"}=1;
        $s{"gcc"}=1;
        $s{"gca"}=1;
        $s{"gcg"}=1;
        
        $s{"tat"}=1;
        $s{"tac"}=1;
        $s{"taa"}=undef;
        $s{"tag"}=undef;
        
        $s{"cat"}=$R/(1+$R);
        $s{"cac"}=$R/(1+$R);
        $s{"caa"}=$R/(1+$R);
        $s{"cag"}=$R/(1+$R);
        
        $s{"aat"}=$R/(1+$R);
        $s{"aac"}=$R/(1+$R);
        $s{"aaa"}=$R/(1+$R);
        $s{"aag"}=$R/(1+$R);
        
        $s{"gat"}=$R/(1+$R);
        $s{"gac"}=$R/(1+$R);
        $s{"gaa"}=$R/(1+$R);
        $s{"gag"}=$R/(1+$R);
        
        $s{"tgt"}=$R/(0.5+$R);
        $s{"tgc"}=$R/(0.5+$R);
        $s{"tga"}=undef;
        $s{"tgg"}=0;
                    
        $s{"cgt"}=1;
        $s{"cgc"}=1;
        $s{"cga"}=1.5;
        $s{"cgg"}=(1.5+$R)/(1+$R);
        
        $s{"agt"}=$R/(1+$R);
        $s{"agc"}=$R/(1+$R);
        $s{"aga"}=1/(1+2*$R)+$R/(1+$R);
        $s{"agg"}=(0.5+$R)/(1+$R);
        
        $s{"ggt"}=1;
        $s{"ggc"}=1;
        $s{"gga"}=1;
        $s{"ggg"}=1;
        return %s;
}
sub sts #synonimous transitions between two codons that differ at <=1 positions
{
	my $diff_n=codon_diff($_[0],$_[1]);
	unless (defined($diff_n)) {print "sts=undef\n"; return undef;}
	if ($diff_n==0) 
	{
		#print "sts=0\n"; 
		return 0;
	}
	if ($diff_n>1) 
	{
		#print "sts=undef\n"; 
		return undef;
	}
	if ($gct{$_[0]} cmp $gct{$_[1]}) 
	{
		#print "sts=0\n"; 
		return 0;
	}
	my @c1=split "",$_[0];
	my @c2=split "",$_[1];
	for (my $j=0;$j<3;$j++)
	{
		if ($c1[$j] cmp $c2[$j])
		{
			if (ts($c1[$j],$c2[$j])) 
			{
				#print "sts=1\n"; 
				return 1;
			}
			else 
			{
				#print "sts=0\n"; 
				return 0;
			}
		}
	}
}

sub nts #nonsynonimous transitions between two codons that differ at <=1 positions
{
	my $diff_n=codon_diff($_[0],$_[1]);
	unless (defined($diff_n)) 
	{
		#print "nts=undef\n"; 
		return undef;
	}
	if ($diff_n==0) 
	{
		#print "nts=0\n"; 
		return 0;
	}
	if ($diff_n>1) 
	{
		#print "nts=undef\n"; 
		return undef;
	}
	unless ($gct{$_[0]} cmp $gct{$_[1]}) 
	{
		#print"nts=0\n"; 
		return 0;
	}
	my @c1=split "",$_[0];
	my @c2=split "",$_[1];
	for (my $j=0;$j<3;$j++)
	{
		if ($c1[$j] cmp $c2[$j])
		{
			if (ts($c1[$j],$c2[$j])) 
			{
				#print "nts=0\n"; 
				return 1;
			}
			else 
			{
				#print "nts=0\n"; 
				return 0;
			}
		}
	}
}

sub stv #synonimous transvertions between two codons that differ at <=1 positions
{
	my $diff_n=codon_diff($_[0],$_[1]);
	unless (defined($diff_n)) 
	{
		#print "stv=undef\n"; 
		return undef;
	}
	if ($diff_n==0) 
	{
		#print "stv=0\n"; 
		return 0;
	}
	if ($diff_n>1) 
	{
		#print "stv=undef\n"; 
		return undef;
	}
	if ($gct{$_[0]} cmp $gct{$_[1]}) 
	{
		#print "stv=0\n"; 
		return 0;
	}
	my @c1=split "",$_[0];
	my @c2=split "",$_[1];
	for (my $j=0;$j<3;$j++)
	{
		if ($c1[$j] cmp $c2[$j])
		{
			if (tv($c1[$j],$c2[$j])) 
			{
				#print "stv=1\n"; 
				return 1;
			}
			else 
			{
				#print "stv=0\n"; 
				return 0;
			}
		}
	}
}

sub ntv #nonsynonimous transvertions between two codons that differ at <=1 positions
{
	my $diff_n=codon_diff($_[0],$_[1]);
	unless (defined($diff_n)) 
	{
		#print "ntv=undef\n"; 
		return undef;
	}
	if ($diff_n==0) 
	{
		#print "ntv=0\n"; 
		return 0;
	}
	if ($diff_n>1) 
	{
		#print "ntv=undef\n"; 
		return undef;
	}
	unless ($gct{$_[0]} cmp $gct{$_[1]}) 
	{
		#print "ntv=0\n"; 
		return 0;
	}
	my @c1=split "",$_[0];
	my @c2=split "",$_[1];
	for (my $j=0;$j<3;$j++)
	{
		if ($c1[$j] cmp $c2[$j])
		{
			if (tv($c1[$j],$c2[$j])) 
			{
				#print "ntv=1\n"; 
				return 1;
			}
			else 
			{
				#print "ntv=0\n"; 
				return 0;
			}
		}
	}
}


1;
