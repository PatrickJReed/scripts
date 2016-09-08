#!/usr/bin/perl -w
use strict;
use warnings;

my $snp138_coding	= "./dbSNP138_Coding.txt"; 
my $refmrna_fa		= "./refMrna.fa";
my $refFlat_file	= "/glusterfs/users/caseybrown/SEQanalysis/refFlat.hg18.txt";

#Input dbsnp_Coding_file and VCF_file

#load ref flat file

#load vcf file, only keep nonsense mutations

#load dbSNP file

# convert dbsnp to protein space

#foreach variant: check if its position has any SNP associated w/ it check if it shares any SNPs in protein Space.

# Build hash for heterozygote genotype conversion
my %iupac = (
	R => ["A", "G"],
	Y => ["C", "T"],
	M => ["A", "C"],
	K => ["T", "G"],
	W => ["A", "T"],
	S => ["C", "G"]
);

my %iupac_rev = (
	AG => "R",
	CT => "Y",
	AC => "M",
	GT => "K",
	AT => "W",
	CG => "S",
);

# Build hash for codon table
my %DNA_code = (
'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A', 'TTA' => 'L',
'TTG' => 'L', 'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L',
'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R', 'AGA' => 'R',
'AGG' => 'R', 'AAA' => 'K', 'AAG' => 'K', 'AAT' => 'N', 'AAC' => 'N',
'ATG' => 'M', 'GAT' => 'D', 'GAC' => 'D', 'TTT' => 'F', 'TTC' => 'F',
'TGT' => 'C', 'TGC' => 'C', 'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P',
'CCG' => 'P', 'CAA' => 'Q', 'CAG' => 'Q', 'TCT' => 'S', 'TCC' => 'S',
'TCA' => 'S', 'TCG' => 'S', 'AGT' => 'S', 'AGC' => 'S', 'GAA' => 'E',
'GAG' => 'E', 'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',
'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G', 'TGG' => 'W',
'CAT' => 'H', 'CAC' => 'H', 'TAT' => 'Y', 'TAC' => 'Y', 'ATT' => 'I',
'ATC' => 'I', 'ATA' => 'I', 'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V',
'GTG' => 'V', 'TAA' => '*', 'TGA' => '*', 'TAG' => '*'
);
