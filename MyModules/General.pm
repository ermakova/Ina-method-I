package MyModules::General; 
#
# Please cite:
# 1. Yulia A. Medvedeva, Marina V. Fridman, Nina J. Oparina, Dmitri B. Malko, Ekaterina O. Ermakova, Ivan V. Kulakovskiy, Vsevolod J. Makeev
# Evidence for transcriptional regulation by non-5' CpG islands in the human genome. BMC Genomics, 2008
# 2. Ina Y 
# New methods for estimating the numbers of synonymous and nonsynonymous substitutions. J Mol Evol, 1995
# 
use strict;
use vars qw(@ISA @EXPORT $VERSION @all_aa %gct %vargct %nuc %pyr %pur %digit);
use Exporter;
@ISA=qw(Exporter);
@EXPORT=qw(min codon_diff ts tv @all_aa %gct %vargct %nuc %pyr %pur %digit);
$VERSION=1.01;


my @all_digit=qw/0 1 2 3 4 5 6 7 8 9/;
%digit=();
for (@all_digit) {$digit{$_}=1;}
my @all_nuc=qw/c t C T a g A G/;
%nuc=();
for (@all_nuc) {$nuc{$_}=1;}
my @all_pyr=qw/c t C T/;
%pyr=();
for (@all_pyr) {$pyr{$_}=1;}
my @all_pur=qw/a g A G/;
%pur=();
for (@all_pur) {$pur{$_}=1;}

my @all_aa=qw/Lys Arg His Asp Glu Gly Asn Gln Cys Ser Thr Tyr Ala Val Leu Ile Met Pro Phe Trp Stop/;

%gct=(); #Universal genetic code table, three-letter

$gct{"ttt"}="Phe";
$gct{"ttc"}="Phe";
$gct{"tta"}="Leu";
$gct{"ttg"}="Leu";

$gct{"ctt"}="Leu";
$gct{"ctc"}="Leu";
$gct{"cta"}="Leu";
$gct{"ctg"}="Leu";

$gct{"att"}="Ile";
$gct{"atc"}="Ile";
$gct{"ata"}="Ile";
$gct{"atg"}="Met";

$gct{"gtt"}="Val";
$gct{"gtc"}="Val";
$gct{"gta"}="Val";
$gct{"gtg"}="Val";

$gct{"tct"}="Ser";
$gct{"tcc"}="Ser";
$gct{"tca"}="Ser";
$gct{"tcg"}="Ser";

$gct{"cct"}="Pro";
$gct{"ccc"}="Pro";
$gct{"cca"}="Pro";
$gct{"ccg"}="Pro";

$gct{"act"}="Thr";
$gct{"acc"}="Thr";
$gct{"aca"}="Thr";
$gct{"acg"}="Thr";

$gct{"gct"}="Ala";
$gct{"gcc"}="Ala";
$gct{"gca"}="Ala";
$gct{"gcg"}="Ala";

$gct{"tat"}="Tyr";
$gct{"tac"}="Tyr";
$gct{"taa"}="Stop";
$gct{"tag"}="Stop";

$gct{"cat"}="His";
$gct{"cac"}="His";
$gct{"caa"}="Gln";
$gct{"cag"}="Gln";

$gct{"aat"}="Asn";
$gct{"aac"}="Asn";
$gct{"aaa"}="Lys";
$gct{"aag"}="Lys";

$gct{"gat"}="Asp";
$gct{"gac"}="Asp";
$gct{"gaa"}="Glu";
$gct{"gag"}="Glu";

$gct{"tgt"}="Cys";
$gct{"tgc"}="Cys";
$gct{"tga"}="Stop";
$gct{"tgg"}="Trp";
            
$gct{"cgt"}="Arg";
$gct{"cgc"}="Arg";
$gct{"cga"}="Arg";
$gct{"cgg"}="Arg";

$gct{"agt"}="Ser";
$gct{"agc"}="Ser";
$gct{"aga"}="Arg";
$gct{"agg"}="Arg";

$gct{"ggt"}="Gly";
$gct{"ggc"}="Gly";
$gct{"gga"}="Gly";
$gct{"ggg"}="Gly";

%vargct=();  #Universal genetic code table, one-letter

$vargct{"ttt"}="F";
$vargct{"ttc"}="F";
$vargct{"tta"}="L";
$vargct{"ttg"}="L";

$vargct{"ctt"}="L";
$vargct{"ctc"}="L";
$vargct{"cta"}="L";
$vargct{"ctg"}="L";

$vargct{"att"}="I";
$vargct{"atc"}="I";
$vargct{"ata"}="I";
$vargct{"atg"}="M";

$vargct{"gtt"}="V";
$vargct{"gtc"}="V";
$vargct{"gta"}="V";
$vargct{"gtg"}="V";

$vargct{"tct"}="S";
$vargct{"tcc"}="S";
$vargct{"tca"}="S";
$vargct{"tcg"}="S";

$vargct{"cct"}="P";
$vargct{"ccc"}="P";
$vargct{"cca"}="P";
$vargct{"ccg"}="P";

$vargct{"act"}="T";
$vargct{"acc"}="T";
$vargct{"aca"}="T";
$vargct{"acg"}="T";

$vargct{"gct"}="A";
$vargct{"gcc"}="A";
$vargct{"gca"}="A";
$vargct{"gcg"}="A";

$vargct{"tat"}="Y";
$vargct{"tac"}="Y";
$vargct{"taa"}="*";
$vargct{"tag"}="*";

$vargct{"cat"}="H";
$vargct{"cac"}="H";
$vargct{"caa"}="Q";
$vargct{"cag"}="Q";

$vargct{"aat"}="N";
$vargct{"aac"}="N";
$vargct{"aaa"}="K";
$vargct{"aag"}="K";

$vargct{"gat"}="D";
$vargct{"gac"}="D";
$vargct{"gaa"}="E";
$vargct{"gag"}="E";

$vargct{"tgt"}="C";
$vargct{"tgc"}="C";
$vargct{"tga"}="*";
$vargct{"tgg"}="W";
            
$vargct{"cgt"}="R";
$vargct{"cgc"}="R";
$vargct{"cga"}="R";
$vargct{"cgg"}="R";

$vargct{"agt"}="S";
$vargct{"agc"}="S";
$vargct{"aga"}="R";
$vargct{"agg"}="R";

$vargct{"ggt"}="G";
$vargct{"ggc"}="G";
$vargct{"gga"}="G";
$vargct{"ggg"}="G";

sub min
{
	unless (defined($_[0]) and defined($_[1])) {return undef;}
	if ($_[0]<=$_[1]) {return $_[0];}
	else {return $_[1];}
}

sub codon_diff 
#number of nucleotide differences between codons, 
#undef if at least one is stop or invalid
#two arguments: $codon1, $codon2
{
	unless (defined($gct{$_[0]}) and defined($gct{$_[1]}))
	{
		return undef;
	}
	unless (($gct{$_[0]} cmp "Stop") and ($gct{$_[1]} cmp "Stop"))
	{
		return undef;
	}
	my @c1=split "",$_[0];
	my @c2=split "",$_[1];
	my $n=0;
	if ($c1[0] cmp $c2[0]) {$n++;}
	if ($c1[1] cmp $c2[1]) {$n++;}
	if ($c1[2] cmp $c2[2]) {$n++;}
	$n;
}

sub ts #transition
{
	unless (defined($_[0]) and defined($_[1])) {return undef;};
	unless (defined($nuc{$_[0]}) and defined($nuc{$_[1]})) {return undef;};
	unless($_[0] cmp $_[1]) {return 0;}
	if ($pyr{$_[0]} and $pyr{$_[1]}){return 1;}
	if ($pur{$_[0]} and $pur{$_[1]}){return 1;}
	0;
}

sub tv #transversion
{
	unless (defined($_[0]) and defined($_[1])) {return undef;};
	unless (defined($nuc{$_[0]}) and defined($nuc{$_[1]})) {return undef;};
	if ($pur{$_[0]} and $pyr{$_[1]}){return 1;}
	if ($pyr{$_[0]} and $pur{$_[1]}){return 1;}
	0;
}
1;
