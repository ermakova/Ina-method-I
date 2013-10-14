package MyModules::Identity;
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
@ISA=qw(Exporter);
@EXPORT=qw(IDn IDa);
$VERSION=1.00;
use MyModules::General;


sub IDn ($$) 
{
	my $l=min(length($_[0]),length($_[1]));
	my $s=0;
	my $z=0;
	for (my $j=0;$j<$l;$j++)
	{
		my $a=lc(substr ($_[0],$j,1));
		my $b=lc(substr ($_[1],$j,1));
		if ($nuc{$a} and $nuc{$b})
		{
			$z++;
			if ($a eq $b)
			{
				$s++;
			}
		}
	}
	($s,$z);
}

sub IDa ($$) 
{
	my $l=min(length($_[0]),length($_[1]));
	my $s=0;
	my $z=0;
	for (my $j=0;$j+2<$l;$j+=3)
	{
		my $codon1=lc(substr($_[0],$j,3)); 
		my $codon2=lc(substr($_[0],$j,3)); 
		if (defined($gct{$codon1}) and defined($gct{$codon2}))
		{
			$z++;
			if ($gct{$codon1} eq $gct{$codon2})
			{
				$s++;
			}
		}
	}
	($s,$z);
}

1;
