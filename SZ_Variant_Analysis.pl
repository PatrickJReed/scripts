#!/usr/bin/perl -w
use strict;
use warnings;

use Text::CSV;
use Array::Utils qw(:all);
use Statistics::Basic qw(:all);
use List::Util qw[min max];


my $SZ_Affected		= "./Affected.avinput.exome_summary.csv";
my $SZ_Control		= "./Control.avinput.exome_summary.csv";
my $dbSNP_135		= "./dbSNP_135.hg18.avinput.exome_summary.csv";
my $RefSeq			= "./hg18_Exons.txt";
my $Mendel_Genes	= "./genes_Mendelian.txt";
my $Pathogenic_snps = "./Pathogenic_SNPs.txt"; 
my @LineInput;
my %Affected_Data;
my %Control_Data;
my %db135_Data;
my %Mendel;
my %GeneData;
my %rsGERPs;
my @tmp;
#Load Annotation Files
my $csv = Text::CSV->new({ sep_char => ',' });

#Load Affected File
open (LIST, '<' , $SZ_Affected) || die "File not found\n";     
     while (my $line = <LIST>) {
     chomp $line;
     if ($csv->parse($line)){
     my @LineInput= $csv->fields();
     my $size = scalar(@LineInput)-1; 
     	if ($LineInput[22] =~ /^[\+-]*[0-9]*\.*[0-9]*$/ && $LineInput[22] !~ /^[\. ]*$/) {
     		if ($LineInput[2] ne "synonymous SNV") {
         $Affected_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"Location"} = $LineInput[0];
         $Affected_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"Class"} = $LineInput[2];
         $Affected_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"AA_Change"} = $LineInput[3];
         $Affected_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"ESP_Freq"} = $LineInput[6];
         @{$Affected_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"1KG_Freq"}} = @LineInput[7,8,9];
         $Affected_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"SNP_ID"} = $LineInput[10];
         $Affected_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"GERP"} = $LineInput[22];
         $Affected_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"SIFT"} = $LineInput[14];
         @{$Affected_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"Delteriousnes"}} = @LineInput[12,14,16,18,20];
         $Affected_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"Ref"} = $LineInput[26];
         $Affected_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"Obs"} = $LineInput[27];
         $Affected_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"Quality"} = $LineInput[28];
         @{$Affected_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"Genotypes"}} = @LineInput[38..$size];
                	}		
        		}
             }
         };
close(LIST);


#Load Control File
open (LIST, '<' , $SZ_Control) || die "File not found\n";     
     while (my $line = <LIST>) {
     chomp $line;
     if ($csv->parse($line)){
     my @LineInput= $csv->fields();
     my $size = scalar(@LineInput)-1;
     	if ($LineInput[22] =~ /^[\+-]*[0-9]*\.*[0-9]*$/ && $LineInput[22] !~ /^[\. ]*$/) {
     		if ($LineInput[2] ne "synonymous SNV") {
         $Control_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"Location"} = $LineInput[0];
         $Control_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"Class"} = $LineInput[2];
         $Control_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"AA_Change"} = $LineInput[3];
         $Control_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"ESP_Freq"} = $LineInput[6];
         @{$Control_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"1KG_Freq"}} = @LineInput[7,8,9];
         $Control_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"SNP_ID"} = $LineInput[10];
         $Control_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"GERP"} = $LineInput[22];
         $Control_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"SIFT"} = $LineInput[14];
         @{$Control_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"Delteriousnes"}} = @LineInput[12,14,16,18,20];
         $Control_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"Ref"} = $LineInput[26];
         $Control_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"Obs"} = $LineInput[27];
         $Control_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"Quality"} = $LineInput[34];
         @{$Control_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"Genotypes"}} = @LineInput[38..$size];
     			}
     		}
     	}
     };
close(LIST);


#Load dbSNP 135
open (LIST, '<' , $dbSNP_135) || die "File not found\n";     
     while (my $line = <LIST>) {
     chomp $line;
     if ($csv->parse($line)){
     my @LineInput= $csv->fields();
     my $size = scalar(@LineInput)-1;
     	if ($LineInput[22] =~ /^[\+-]*[0-9]*\.*[0-9]*$/ && $LineInput[22] !~ /^[\. ]*$/) {
     	if ($LineInput[2] ne "synonymous SNV") {
         $db135_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"Location"} = $LineInput[0];
         $db135_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"Class"} = $LineInput[2];
        $db135_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"AA_Change"} = $LineInput[3];
         $db135_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"ESP_Freq"} = $LineInput[6];
         @{$db135_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"1KG_Freq"}} = @LineInput[7,8,9];
         $db135_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"SNP_ID"} = $LineInput[10];
         $db135_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"GERP"} = $LineInput[22];
         $db135_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"SIFT"} = $LineInput[14];
        @{$db135_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"Delteriousnes"}} = @LineInput[12,14,16,18,20];
         $db135_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"Ref"} = $LineInput[26];
         $db135_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"Obs"} = $LineInput[27];
         $db135_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"Quality"} = $LineInput[34];
         @{$db135_Data{$LineInput[1]}{"$LineInput[23]\t$LineInput[24]"}{"Genotypes"}} = @LineInput[38..$size];        
     		$rsGERPs{$LineInput[10]}{"Gerp"}=$LineInput[22];
     		$rsGERPs{$LineInput[10]}{"Gene"}=$LineInput[1];
     		$rsGERPs{$LineInput[10]}{"Freq"}=$LineInput[6];
     			}
     		}
     	}
     };
close(LIST);

##################### Define het frequency of each variants in cohorts ###################
my @affected_vars;
my @control_vars;
foreach (keys %Affected_Data) {
my $position;
	foreach $position (keys %{$Affected_Data{$_}}) {
		push(@affected_vars, "$_\t$position");
	};
};		
foreach (keys %Control_Data) {
my $position;
	foreach $position (keys %{$Control_Data{$_}}) {
		push(@control_vars, "$_\t$position");
	};
};
foreach (@affected_vars) {
my @tmp = split(/\t/,$_);
my $gts;
my $hetc = 0;
my $homc = 0;
my $norm = 0;	
	foreach $gts (@{$Affected_Data{"$tmp[0]"}{"$tmp[1]\t$tmp[2]"}{"Genotypes"}}) {
		my @gt = split(":", $gts);
			if ($gt[0] eq "0/1") { $hetc++; }
			elsif ($gt[0] eq "1/1") { $homc++; }
			elsif ($gt[0] eq "0/0") { $norm++; }
			}	
		$Affected_Data{"$tmp[0]"}{"$tmp[1]\t$tmp[2]"}{"Het_Freq"} = $hetc;
		$Affected_Data{"$tmp[0]"}{"$tmp[1]\t$tmp[2]"}{"Hom_Freq"} = $homc;
		$Affected_Data{"$tmp[0]"}{"$tmp[1]\t$tmp[2]"}{"Nrm_Freq"} = $norm; 	
};	
foreach (@control_vars) {
my @tmp = split(/\t/,$_);
my $gts;
my $hetc = 0;
my $homc = 0;	
my $norm = 0;	
	foreach $gts (@{$Control_Data{"$tmp[0]"}{"$tmp[1]\t$tmp[2]"}{"Genotypes"}}) {		
		my @gt = split(":", $gts);	
			if ($gt[0] eq "0/1") { $hetc++; }
			elsif ($gt[0] eq "1/1") { $homc++; }
			elsif ($gt[0] eq "0/0") { $norm++; }
			}	
		$Control_Data{"$tmp[0]"}{"$tmp[1]\t$tmp[2]"}{"Het_Freq"} = $hetc; 	
		$Control_Data{"$tmp[0]"}{"$tmp[1]\t$tmp[2]"}{"Hom_Freq"} = $homc;
		$Control_Data{"$tmp[0]"}{"$tmp[1]\t$tmp[2]"}{"Nrm_Freq"} = $norm;	
};

#Load Refseq Data
#open (LIST, $RefSeq) || die "File not found\n";     
#     while (<LIST>) {
#     chomp;
#	 my @Data = split(/\t/, $_);
#	 my @exonStarts = split(/,/, $Data[1]);
#	 my @exonEnds = split(/,/, $Data[2]);
#	 my $codinglength = 0;
#	 my $i = 0;
#	 	foreach (@exonStarts) {
#	 	if (($exonStarts[$i] + $exonEnds[$i]) =~ /^[\+-]*[0-9]*\.*[0-9]*$/ && ($exonStarts[$i] + $exonEnds[$i]) !~ /^[\. ]*$/) {
#	 		$codinglength = $codinglength + abs($exonEnds[$i]-$exonStarts[$i]);
#	 		$i++; 
#	 		};
#	 	};
#	 $GeneData{$Data[3]}{"Coding_Length"} = $codinglength;
#	 	};
#close(LIST);	 
	 
	 
#Load Andrey Gene list
open (LIST, $Mendel_Genes) || die "File not found\n";     
     while (<LIST>) {
     chomp;
     my @meta = split(/\|/, $_);
     my @genes = split(/,/, $meta[2]);
     	foreach my $gene (@genes) {
     	#print $gene."\n"; 	
       	$Mendel{$gene}{$meta[1]}="true";
			}
     };        
close(LIST);

#Load Pathogenic rsIDs
my @pathogenic_snps;
open (LIST, $Pathogenic_snps) || die "File not found\n";     
     while (<LIST>) {
     chomp;
     my @snps = split(' ', $_);
     foreach (@snps) {
     #print "$_"."\n";
     push(@pathogenic_snps, $_);
     		}    
     	};        
close(LIST);

#Print list of pathogenic variants, genes and gerps
my $Pathogenic_Gerps = "./Gerp_Pathogenic_dbSNP.txt";
open(MYOUTPUTFILE, ">$Pathogenic_Gerps");
print MYOUTPUTFILE "Gene_Name"."\t"."SNP_Gerp"."\t"."SNP_Freq"."\n";
foreach  (@pathogenic_snps) {
		 #print MYOUTPUTFILE $_."\t";	
		 if ($rsGERPs{$_}){
		 print MYOUTPUTFILE $rsGERPs{$_}{"Gene"}."\t";
		 print MYOUTPUTFILE $rsGERPs{$_}{"Gerp"}."\t";
		 if ($rsGERPs{$_}{"Freq"}){
		 print MYOUTPUTFILE $rsGERPs{$_}{"Freq"}."\n";
		 }
		 else { print MYOUTPUTFILE "0\n"; }
         };
};
close(MYOUTPUTFILE);

my %non_patho = %rsGERPs;
foreach (@pathogenic_snps) {
	if ($non_patho{$_}) { delete($non_patho{$_}); }
};

my @non_patho = keys(%non_patho);
my @all_snps  = keys(%rsGERPs);
print scalar(@pathogenic_snps)."\n";
print scalar(@all_snps)."\n";
print scalar(@non_patho)."\n";

my $nonPatho_Gerps = "./Gerp_nonPatho_dbSNP.txt";
open(MYOUTPUTFILE, ">$nonPatho_Gerps");
print MYOUTPUTFILE "Gene_Name"."\t"."SNP_Gerp"."\t"."SNP_Freq"."\n";
foreach  (@non_patho) {
		 #print MYOUTPUTFILE $_."\t";	
		 if ($rsGERPs{$_}){
		 print MYOUTPUTFILE $rsGERPs{$_}{"Gene"}."\t";
		 print MYOUTPUTFILE $rsGERPs{$_}{"Gerp"}."\t";
		 if ($rsGERPs{$_}{"Freq"}){
		 print MYOUTPUTFILE $rsGERPs{$_}{"Freq"}."\n";
		 }
		 else { print MYOUTPUTFILE "0\n"; }
         };
};
close(MYOUTPUTFILE);



#Analyze Genes in mendelian list, determine degree of variation for each gene in affected, control, and snp files. Sum gerp/coding length and n=number of variants

my @mendel_genes = keys(%Mendel);
my @affected_genes = keys(%Affected_Data);
my @control_genes = keys(%Control_Data);
my @dbsnp_genes = keys(%db135_Data);

#Analyze intersect of Mendel and affected genes
my @mendel_affected = intersect(@mendel_genes, @affected_genes);
print scalar(@mendel_affected)."\n";
my %mendel_affected_results;
foreach (@mendel_affected) {
	my @variants = keys %{$Affected_Data{$_}};
	my $Gene_Variant_Count = scalar(@variants);
	my @affected_gerps;
	my @affected_hets;
	my @affected_homs;
	my @affected_nrms;
	foreach my $variant (@variants) {
		push(@affected_gerps, $Affected_Data{$_}{$variant}{"GERP"});
		if ($Affected_Data{$_}{$variant}{"Het_Freq"}) {
		push(@affected_hets, $Affected_Data{$_}{$variant}{"Het_Freq"});
		push(@affected_homs, $Affected_Data{$_}{$variant}{"Hom_Freq"});
		push(@affected_nrms, $Affected_Data{$_}{$variant}{"Nrm_Freq"});
		}
	};
	my $cumulative_gerp;
	my $het_count;
	my $hom_count;
	my $nrm_count;
	$cumulative_gerp += $_ for @affected_gerps;
	$het_count += $_ for @affected_hets;
	$hom_count += $_ for @affected_homs;
	$nrm_count += $_ for @affected_nrms;
	$mendel_affected_results{$_}{"cumulative_gerp"} = $cumulative_gerp;
	$mendel_affected_results{$_}{"variant_count"} = $Gene_Variant_Count;
	$mendel_affected_results{$_}{"mean"} = mean(@affected_gerps);
	$mendel_affected_results{$_}{"stdev"} = stddev(@affected_gerps);
	$mendel_affected_results{$_}{"min"} = min(@affected_gerps);
	$mendel_affected_results{$_}{"max"} = max(@affected_gerps);	
	$mendel_affected_results{$_}{"het_count"} = $het_count;
	$mendel_affected_results{$_}{"hom_count"} = $hom_count;
	$mendel_affected_results{$_}{"nrm_count"} = $nrm_count;
};


#Analyze intersect of Mendel and control genes
my @mendel_control = intersect(@mendel_genes, @control_genes);
print scalar(@mendel_control)."\n";
my %mendel_control_results;
foreach (@mendel_control) {
	my @variants = keys %{$Control_Data{$_}};
	my $Gene_Variant_Count = scalar(@variants);
	my @control_gerps;
	my @control_hets;
	my @control_homs;
	my @control_nrms;
	foreach my $variant (@variants) {
		push(@control_gerps, $Control_Data{$_}{$variant}{"GERP"});
		if ($Control_Data{$_}{$variant}{"Het_Freq"}){
		push(@control_hets, $Control_Data{$_}{$variant}{"Het_Freq"});
		push(@control_homs, $Control_Data{$_}{$variant}{"Hom_Freq"});
		push(@control_nrms, $Control_Data{$_}{$variant}{"Nrm_Freq"});
		}
	};
	my $cumulative_gerp;
	my $het_count;
	my $hom_count;
	my $nrm_count;
	$cumulative_gerp += $_ for @control_gerps;
	$het_count += $_ for @control_hets;
	$hom_count += $_ for @control_homs;
	$nrm_count += $_ for @control_nrms;
	$mendel_control_results{$_}{"cumulative_gerp"} = $cumulative_gerp;
	$mendel_control_results{$_}{"variant_count"} = $Gene_Variant_Count;
	$mendel_control_results{$_}{"mean"} = mean(@control_gerps);
	$mendel_control_results{$_}{"stdev"} = stddev(@control_gerps);
	$mendel_control_results{$_}{"min"} = min(@control_gerps);
	$mendel_control_results{$_}{"max"} = max(@control_gerps);
	$mendel_control_results{$_}{"het_count"} = $het_count; 
	$mendel_control_results{$_}{"hom_count"} = $hom_count; 
	$mendel_control_results{$_}{"nrm_count"} = $nrm_count;
};


#Analyze intersect of Mendel and dbsnp genes
my @mendel_dbsnp = intersect(@mendel_genes, @dbsnp_genes);
print scalar(@mendel_dbsnp)."\n";
my %mendel_dbsnp_results;
foreach (@mendel_dbsnp) {
	my @variants = keys %{$db135_Data{$_}};
	my $Gene_Variant_Count = scalar(@variants);
	my @dbsnp_gerps;
	foreach my $variant (@variants) {
		push(@dbsnp_gerps, $db135_Data{$_}{$variant}{"GERP"});
	};
	my $cumulative_gerp;
	$cumulative_gerp += $_ for @dbsnp_gerps;
	$mendel_dbsnp_results{$_}{"cumulative_gerp"} = $cumulative_gerp;
	$mendel_dbsnp_results{$_}{"variant_count"} = $Gene_Variant_Count;
	$mendel_dbsnp_results{$_}{"mean"} = mean(@dbsnp_gerps);
	$mendel_dbsnp_results{$_}{"stdev"} = stddev(@dbsnp_gerps);
	$mendel_dbsnp_results{$_}{"min"} = min(@dbsnp_gerps);
	$mendel_dbsnp_results{$_}{"max"} = max(@dbsnp_gerps);
};


#################### Print Gerp Output analysis ############################

my $Gene_Gerp_Analysis = "./Gene_Gerp_Analysis.txt";
my $var;
open(MYOUTPUTFILE, ">$Gene_Gerp_Analysis");
print MYOUTPUTFILE "Gene_Name"."\t"."Affected_Cumulative_Gerp"."\t"."Affected_Var_Count"."\t"."Affected_Mean_Gerp"."\t"."Affected_Stdev"."\t"."Affected_min_Gerp"."\t"."Affected_max_Gerp"."\t"."Affected_HET"."\t"."Affected_HOM"."\t"."Affected_NRM"."\t"."Control_Cumulative_Gerp"."\t"."Control_Var_Count"."\t"."Control_Mean_Gerp"."\t"."Control_Stdev"."\t"."Control_min_Gerp"."\t"."Control_max_Gerp"."\t"."Control_HET"."\t"."Control_HOM"."\t"."Control_NRM"."\t"."dbSNP135_Cumulative_Gerp"."\t"."dbSNP135_Var_Count"."\t"."dbSNP135_Mean_Gerp"."\t"."dbSNP135_Stdev"."\t"."dbSNP135_min_Gerp"."\t"."dbSNP135_max_Gerp"."\n";
foreach $var (keys %Mendel) {
		 print MYOUTPUTFILE $var."\t";
		 if ($mendel_affected_results{$var}{"het_count"} > $mendel_affected_results{$var}{"hom_count"}){
		 print MYOUTPUTFILE $mendel_affected_results{$var}{"cumulative_gerp"}."\t";
         print MYOUTPUTFILE $mendel_affected_results{$var}{"variant_count"}."\t";
         print MYOUTPUTFILE $mendel_affected_results{$var}{"mean"}."\t";
		 print MYOUTPUTFILE $mendel_affected_results{$var}{"stdev"}."\t";
		 print MYOUTPUTFILE $mendel_affected_results{$var}{"min"}."\t";
		 print MYOUTPUTFILE $mendel_affected_results{$var}{"max"}."\t";
         print MYOUTPUTFILE $mendel_affected_results{$var}{"het_count"}."\t";
         print MYOUTPUTFILE $mendel_affected_results{$var}{"hom_count"}."\t";
         print MYOUTPUTFILE $mendel_affected_results{$var}{"nrm_count"}."\t";
         }
         else {
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
		 };
		 if ($mendel_control_results{$var}{"het_count"} > $mendel_control_results{$var}{"hom_count"}){
		 print MYOUTPUTFILE $mendel_control_results{$var}{"cumulative_gerp"}."\t";
         print MYOUTPUTFILE $mendel_control_results{$var}{"variant_count"}."\t";
         print MYOUTPUTFILE $mendel_control_results{$var}{"mean"}."\t";
	     print MYOUTPUTFILE $mendel_control_results{$var}{"stdev"}."\t";
	     print MYOUTPUTFILE $mendel_control_results{$var}{"min"}."\t";
	     print MYOUTPUTFILE $mendel_control_results{$var}{"max"}."\t";
	     print MYOUTPUTFILE $mendel_control_results{$var}{"het_count"}."\t";
         print MYOUTPUTFILE $mendel_control_results{$var}{"hom_count"}."\t";
         print MYOUTPUTFILE $mendel_control_results{$var}{"nrm_count"}."\t";
         }
         else {
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
		 };
		 if ($mendel_dbsnp_results{$var} and $mendel_control_results{$var} and $mendel_affected_results{$var}){
		 print MYOUTPUTFILE $mendel_dbsnp_results{$var}{"cumulative_gerp"}."\t";
         print MYOUTPUTFILE $mendel_dbsnp_results{$var}{"variant_count"}."\t";
         print MYOUTPUTFILE $mendel_dbsnp_results{$var}{"mean"}."\t";
	     print MYOUTPUTFILE $mendel_dbsnp_results{$var}{"stdev"}."\t";
	     print MYOUTPUTFILE $mendel_dbsnp_results{$var}{"min"}."\t";
	     print MYOUTPUTFILE $mendel_dbsnp_results{$var}{"max"}."\n";

         }
         else {
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\n";
		 };
	

	};
close(MYOUTPUTFILE);

#Analyze affected genes
print scalar(@affected_genes)."\n";
my %affected_results;
foreach (@affected_genes) {
	my @variants = keys %{$Affected_Data{$_}};
	my $Gene_Variant_Count = scalar(@variants);
	#print "$Gene_Variant_Count\n";
	my @affected_gerps;
	my @affected_hets;
	my @affected_homs;
	my @affected_nrms;
	foreach my $variant (@variants) {
		push(@affected_gerps, $Affected_Data{$_}{$variant}{"GERP"});
		if ($Affected_Data{$_}{$variant}{"Het_Freq"}) {
		push(@affected_hets, $Affected_Data{$_}{$variant}{"Het_Freq"});
		push(@affected_homs, $Affected_Data{$_}{$variant}{"Hom_Freq"});
		push(@affected_nrms, $Affected_Data{$_}{$variant}{"Nrm_Freq"});
		}
	};
	my $cumulative_gerp;
	my $het_count;
	my $hom_count;
	my $nrm_count;
	$cumulative_gerp += $_ for @affected_gerps;
	$het_count += $_ for @affected_hets;
	$hom_count += $_ for @affected_homs;
	$nrm_count += $_ for @affected_nrms;
	$affected_results{$_}{"cumulative_gerp"} = $cumulative_gerp;
	$affected_results{$_}{"variant_count"} = $Gene_Variant_Count;
	$affected_results{$_}{"mean"} = mean(@affected_gerps);
	$affected_results{$_}{"stdev"} = stddev(@affected_gerps);
	$affected_results{$_}{"min"} = min(@affected_gerps);
	$affected_results{$_}{"max"} = max(@affected_gerps);
	$affected_results{$_}{"het_count"} = $het_count;
	$affected_results{$_}{"hom_count"} = $hom_count;
	$affected_results{$_}{"nrm_count"} = $nrm_count;
};

#Analyze control genes
print scalar(@control_genes)."\n";
my %control_results;
foreach (@control_genes) {
	my @variants = keys %{$Control_Data{$_}};
	my $Gene_Variant_Count = scalar(@variants);
	my @control_gerps;
	my @control_hets;
	my @control_homs;
	my @control_nrms;
	foreach my $variant (@variants) {
		push(@control_gerps, $Control_Data{$_}{$variant}{"GERP"});
		if ($Control_Data{$_}{$variant}{"Het_Freq"}){
		push(@control_hets, $Control_Data{$_}{$variant}{"Het_Freq"});
		push(@control_homs, $Control_Data{$_}{$variant}{"Hom_Freq"});
		push(@control_nrms, $Control_Data{$_}{$variant}{"Nrm_Freq"});
		}
	};
	my $cumulative_gerp;
	my $het_count;
	my $hom_count;
	my $nrm_count;
	$cumulative_gerp += $_ for @control_gerps;
	$het_count += $_ for @control_hets;
	$hom_count += $_ for @control_homs;
	$nrm_count += $_ for @control_nrms;
	$control_results{$_}{"cumulative_gerp"} = $cumulative_gerp;
	$control_results{$_}{"variant_count"} = $Gene_Variant_Count;
	$control_results{$_}{"mean"} = mean(@control_gerps);
	$control_results{$_}{"stdev"} = stddev(@control_gerps);
	$control_results{$_}{"min"} = min(@control_gerps);
	$control_results{$_}{"max"} = max(@control_gerps);
	$control_results{$_}{"het_count"} = $het_count;
	$control_results{$_}{"hom_count"} = $hom_count;
	$control_results{$_}{"nrm_count"} = $nrm_count;
};

#Analyze dbsnp genes
print scalar(@dbsnp_genes)."\n";
my %dbsnp_results;
foreach (@dbsnp_genes) {
	my @variants = keys %{$db135_Data{$_}};
	my $Gene_Variant_Count = scalar(@variants);
	my @dbsnp_gerps;
	foreach my $variant (@variants) {
		push(@dbsnp_gerps, $db135_Data{$_}{$variant}{"GERP"});
	};
	my $cumulative_gerp;
	$cumulative_gerp += $_ for @dbsnp_gerps;
	$dbsnp_results{$_}{"cumulative_gerp"} = $cumulative_gerp;
	$dbsnp_results{$_}{"variant_count"} = $Gene_Variant_Count;
	$dbsnp_results{$_}{"mean"} = mean(@dbsnp_gerps);
	$dbsnp_results{$_}{"stdev"} = stddev(@dbsnp_gerps);
	$dbsnp_results{$_}{"min"} = min(@dbsnp_gerps);
	$dbsnp_results{$_}{"max"} = max(@dbsnp_gerps);

};


my $Gene_Gerp_Affected = "./Gene_Gerp_Affected.txt";
open(MYOUTPUTFILE, ">$Gene_Gerp_Affected");
print MYOUTPUTFILE "Gene_Name"."\t"."Affected_Cumulative_Gerp"."\t"."Affected_Var_Count"."\t"."Affected_Mean_Gerp"."\t"."Affected_Stdev"."\t"."Affected_min_Gerp"."\t"."Affected_max_Gerp"."\t"."Affected_HET"."\t"."Affected_HOM"."\t"."Affected_NRM"."\n";
foreach $var (keys %Affected_Data) {
		 print MYOUTPUTFILE $var."\t";
		 if ($affected_results{$var}{"het_count"}){
		 print MYOUTPUTFILE $affected_results{$var}{"cumulative_gerp"}."\t";
         print MYOUTPUTFILE $affected_results{$var}{"variant_count"}."\t";
         print MYOUTPUTFILE $affected_results{$var}{"mean"}."\t";
		 print MYOUTPUTFILE $affected_results{$var}{"stdev"}."\t";
         print MYOUTPUTFILE $affected_results{$var}{"min"}."\t";
         print MYOUTPUTFILE $affected_results{$var}{"max"}."\t";
         print MYOUTPUTFILE $affected_results{$var}{"het_count"}."\t";
         print MYOUTPUTFILE $affected_results{$var}{"hom_count"}."\t";
         print MYOUTPUTFILE $affected_results{$var}{"nrm_count"}."\n";
         }
         else {
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\n";
		 };
};
close(MYOUTPUTFILE);


my $Gene_Gerp_Control = "./Gene_Gerp_Control.txt";
open(MYOUTPUTFILE, ">$Gene_Gerp_Control");
print MYOUTPUTFILE "Gene_Name"."\t"."Control_Cumulative_Gerp"."\t"."Control_Var_Count"."\t"."Control_Mean_Gerp"."\t"."Control_Stdev"."\t"."Control_min_Gerp"."\t"."Control_max_Gerp"."\t"."Control_HET"."\t"."Control_HOM"."\t"."Control_NRM"."\n";
foreach $var (keys %Control_Data) {
		 print MYOUTPUTFILE $var."\t";
		 if ($control_results{$var}{"het_count"}){
		 print MYOUTPUTFILE $control_results{$var}{"cumulative_gerp"}."\t";
         print MYOUTPUTFILE $control_results{$var}{"variant_count"}."\t";
         print MYOUTPUTFILE $control_results{$var}{"mean"}."\t";
	     print MYOUTPUTFILE $control_results{$var}{"stdev"}."\t";
	     print MYOUTPUTFILE $control_results{$var}{"min"}."\t";
         print MYOUTPUTFILE $control_results{$var}{"max"}."\t";
         print MYOUTPUTFILE $control_results{$var}{"het_count"}."\t";
         print MYOUTPUTFILE $control_results{$var}{"hom_count"}."\t";
         print MYOUTPUTFILE $control_results{$var}{"nrm_count"}."\n";
         }
         else {
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\n";
		 };
};
close(MYOUTPUTFILE);


my $Gene_Gerp_dbsnp = "./Gene_Gerp_dbsnp.txt";
open(MYOUTPUTFILE, ">$Gene_Gerp_dbsnp");
print MYOUTPUTFILE "Gene_Name"."\t"."dbSNP_Cumulative_Gerp"."\t"."dbSNP_Var_Count"."\t"."dbSNP_Mean_Gerp"."\t"."dbSNP_Stdev"."\n";
foreach $var (keys %db135_Data) {
		 print MYOUTPUTFILE $var."\t";
		 if ($dbsnp_results{$var}){
		 print MYOUTPUTFILE $dbsnp_results{$var}{"cumulative_gerp"}."\t";
         print MYOUTPUTFILE $dbsnp_results{$var}{"variant_count"}."\t";
         print MYOUTPUTFILE $dbsnp_results{$var}{"mean"}."\t";
	     print MYOUTPUTFILE $dbsnp_results{$var}{"stdev"}."\n";
         }
         else {
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\t";
         print MYOUTPUTFILE "NA"."\n";
		 };
};
close(MYOUTPUTFILE);


#Analyze intersect of affected and control genes
my @affected_control = intersect(@control_genes, @affected_genes);
print scalar(@affected_control)."\n";
my $Affected_vs_Control = "./Affected_vs_Control_GerpAnalysis2.txt";
open(MYOUTPUTFILE, ">$Affected_vs_Control");
print MYOUTPUTFILE "Gene_Name"."\t"."Affected_Cumulative_Gerp"."\t"."Affected_Var_Count"."\t"."Affected_Mean_Gerp"."\t"."Affected_Stdev"."\t"."Affected_min_Gerp"."\t"."Affected_max_Gerp"."\t"."Affected_HET"."\t"."Affected_HOM"."\t"."Affected_NRM"."\t"."Control_Cumulative_Gerp"."\t"."Control_Var_Count"."\t"."Control_Mean_Gerp"."\t"."Control_Stdev"."\t"."Control_min_Gerp"."\t"."Control_max_Gerp"."\t"."Control_HET"."\t"."Control_HOM"."\t"."Control_NRM"."\n";
foreach (@affected_control) {
		 if ($affected_results{$_}{"het_count"} and $control_results{$_}{"het_count"}){
		 print MYOUTPUTFILE $_."\t";
		 print MYOUTPUTFILE $affected_results{$_}{"cumulative_gerp"}."\t";
         print MYOUTPUTFILE $affected_results{$_}{"variant_count"}."\t";
         print MYOUTPUTFILE $affected_results{$_}{"mean"}."\t";
		 print MYOUTPUTFILE $affected_results{$_}{"stdev"}."\t";
		 print MYOUTPUTFILE $affected_results{$_}{"min"}."\t";
         print MYOUTPUTFILE $affected_results{$_}{"max"}."\t";
         print MYOUTPUTFILE $affected_results{$_}{"het_count"}."\t";
         print MYOUTPUTFILE $affected_results{$_}{"hom_count"}."\t";
         print MYOUTPUTFILE $affected_results{$_}{"nrm_count"}."\t";
		 print MYOUTPUTFILE $control_results{$_}{"cumulative_gerp"}."\t";
         print MYOUTPUTFILE $control_results{$_}{"variant_count"}."\t";
         print MYOUTPUTFILE $control_results{$_}{"mean"}."\t";
	     print MYOUTPUTFILE $control_results{$_}{"stdev"}."\t";
	     print MYOUTPUTFILE $control_results{$_}{"min"}."\t";
         print MYOUTPUTFILE $control_results{$_}{"max"}."\t";
         print MYOUTPUTFILE $control_results{$_}{"het_count"}."\t";
         print MYOUTPUTFILE $control_results{$_}{"hom_count"}."\t";
         print MYOUTPUTFILE $control_results{$_}{"nrm_count"}."\n";
         }
      };
         
close(MYOUTPUTFILE);        








