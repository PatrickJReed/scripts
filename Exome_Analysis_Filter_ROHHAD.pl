#!/usr/bin/perl -w
use strict;
use warnings;
my $Input_File1		= "./Haplotypes_K87_GVCFs.recal.BOTH.annovar.PASS.genic.txt";
my $Output_File		= "./Haplotypes_K87_GVCFs.recal.BOTH.annovar.PASS.genic.Annotated.txt";
my (@tmp1, @tmp2, @IDs);
my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, $S1, $S2, $S3, $S4, $S5, $S6, $S7, %DATA1, @varID, $var, $PL, $GT, $DP, $GQ, $pPL, $pGT, $pDP, $pGQ, $aGT, $aPL, $aDP, $aGQ, $uGT, $uPL, $uDP, $uGQ, $dum, $depth, %Genotypes, $Ale1, $Ale2, $pAle1, $pAle2, $aAle1, $aAle2);


my @UID = ("450","452","925","926","927","929");

 
 
my %sample_info = (
 	
 	$UID[0], { "Sex", "M", "Mother", $UID[4], "Father", "NS", "Generation", "2", "Affected", "N" },
 	$UID[1], { "Sex", "F", "Mother", $UID[2], "Father", "NS", "Generation", "3", "Affected", "Y" },
 	$UID[2], { "Sex", "F", "Mother", $UID[4], "Father", "NS", "Generation", "2", "Affected", "Y" },
 	$UID[3], { "Sex", "F", "Mother", $UID[4], "Father", "NS", "Generation", "2", "Affected", "Y" },
 	$UID[4], { "Sex", "F", "Mother", "NS", "Father", "NS", "Generation", "1", "Affected", "Y" },
 	$UID[5], { "Sex", "F", "Mother", $UID[4], "Father", "NS", "Generation", "2", "Affected", "N" }
		);
 

open (LIST, $Input_File1) || die "File not found\n";     
     while (<LIST>) {
     chomp;
         @tmp1 = split(/\t/, $_);
         if ($tmp1[1] ne "Start"){
		 $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Chr"} = $tmp1[0];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Start"} = $tmp1[1];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"End"} = $tmp1[2];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Ref"} = $tmp1[3];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Alt"} = $tmp1[4];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Func.refGene"} = $tmp1[5];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Gene.refGene"} = $tmp1[6];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"GeneDetail.refGene"} = $tmp1[7];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"ExonicFunc.refGene"} = $tmp1[8];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"AAChange.refGene"} = $tmp1[9];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"exac02"} = $tmp1[10];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"CG69"} = $tmp1[11];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"ESP6500si_ALL"} = $tmp1[12];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"1000g2014sep_all"} = $tmp1[13]; 
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"snp141"} = $tmp1[14];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"snp138NonFlagged"} = $tmp1[15];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"clinvar_20140211"} = $tmp1[16];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"SIFT_score"} = $tmp1[17];      
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"SIFT_pred"} = $tmp1[18];       
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Polyphen2_HDIV_score"} = $tmp1[19];    
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Polyphen2_HDIV_pred"} = $tmp1[20];     
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Polyphen2_HVAR_score"} = $tmp1[21];    
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Polyphen2_HVAR_pred"} = $tmp1[22];     
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"LRT_score"} = $tmp1[23];       
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"LRT_pred"} = $tmp1[24];    
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"MutationTaster_score"} = $tmp1[25];    
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"MutationTaster_pred"} = $tmp1[26];     
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"MutationAssessor_score"} = $tmp1[27];  
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"MutationAssessor_pred"} = $tmp1[28];   
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"FATHMM_score"} = $tmp1[29];    
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"FATHMM_pred"} = $tmp1[30];     
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"RadialSVM_score"} = $tmp1[31]; 
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"RadialSVM_pred"} = $tmp1[32];  
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"LR_score"} = $tmp1[33];        
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"LR_pred"} = $tmp1[34]; 
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"VEST3_score"} = $tmp1[35];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"CADD_raw"} = $tmp1[36];        
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"CADD_pred"} = $tmp1[37];    
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"GERP++_RS"} = $tmp1[38];       
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"phyloP46way_placental"} = $tmp1[39];   
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"phyloP100way_vertebrate"} = $tmp1[40]; 
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"SiPhy_29way_logOdds"} = $tmp1[41];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"OMIM_Regions"} = $tmp1[42];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"pubsBlat"} = $tmp1[43];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"tfbs"} = $tmp1[44];         
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"gerp++gt2"} = $tmp1[45];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"gwasCatalog"} = $tmp1[46];
		 $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"wgRna"} = $tmp1[47];   
		 $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"targetScanS"} = $tmp1[48];         
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"genomicSuperDups"} = $tmp1[49];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"CHROM"} = $tmp1[50];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"POS"} = $tmp1[51];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"ID"} = $tmp1[52];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"REF"} = $tmp1[53];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"ALT"} = $tmp1[54];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"QUAL"} = $tmp1[55];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"FILTER"} = $tmp1[56];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"INFO"} = $tmp1[57];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"FORMAT"} = $tmp1[58];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Mother"} = $tmp1[59];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Father"} = $tmp1[60];
         $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"Child"} = $tmp1[61];




      	} 
     
	if ($DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"1000g2014sep_all"}) 
		{if ($DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"1000g2014sep_all"} > .0001) {delete $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"};}}
	elsif ($DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"ESP6500si_ALL"}) 
		{if ($DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"ESP6500si_ALL"} > .0001) {delete $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"};}}
	elsif ($DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"CG69"}) 
		{if ($DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"CG69"} > .0001) {delete $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"};}}
	elsif ($DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"exac02"}) 
		{if ($DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"}{"exac02"} > .0001) {delete $DATA1{"$tmp1[0]\t$tmp1[1]\t$tmp1[2]"};}}		
};
close(LIST);

@varID =keys(%DATA1);
print scalar(@varID);
print "\n";


# Filter variants by presence of homozygosity among affected Individuals 

if (scalar(@Affected) > 0) {
			foreach (@Affected) {
			@varID =keys(%DATA1);
			foreach $var (@varID) {
				if ($DATA1{$var}{$_}) {
				my @INFO = split(':', $DATA1{$var}{$_});
					my $GT = $INFO[0];
					my $AD = $INFO[1];
					my $DP = $INFO[2];
					my $GQ = $INFO[3];
					if (($GT ne "./.") and ($AD ne ".")) {
					my ($Ale1, $Ale2) = split('/', $GT);
						if (($GT ne "./.")) {
									if ($GQ >= 20) {
										if ($Ale1 eq $Ale2) {
											delete $DATA1{$var};
										}
									}			
								}	
							}	          
						}
					}
				}
			}				

@varID =keys(%DATA1);
print scalar(@varID);
print "\n";


#Filter variants by presence of variant
if (scalar(@Unaffected) > 0) {
			foreach (@Unaffected) {
			@varID =keys(%DATA1);
			foreach $var (@varID) {
				if ($DATA1{$var}{$_}) {
				my @INFO = split(':', $DATA1{$var}{$_});
					my $GT = $INFO[0];
					my $AD = $INFO[1];
					my $DP = $INFO[2];
					my $GQ = $INFO[3];
					if (($GT ne "./.")) { 			
						if ($GQ >= 20) { 
							if (($GT eq ("0/1" or "1/1"))){
									delete $DATA1{$var};
								}
							}			
						}	
					}
				}
			}			          			
		}
@varID =keys(%DATA1);
print scalar(@varID);
print "\n";

# Filter variants by pairwise comparison of each affected to unaffected subjects
foreach my $aftd (@Affected) {
	foreach my $unaftd (@Unaffected) {
		@varID =keys(%DATA1);
		foreach $var (@varID) {
			if ($DATA1{$var}{$aftd} and $DATA1{$var}{$unaftd}) {
					my @aINFO = split(':', $DATA1{$var}{$aftd});
					my $aGT = $aINFO[0];
					my $aAD = $aINFO[1];
					my $aDP = $aINFO[2];
					my $aGQ = $aINFO[3];
					my @uINFO = split(':', $DATA1{$var}{$unaftd});
					my $uGT = $uINFO[0];
					my $uAD = $uINFO[1];
					my $uDP = $uINFO[2];
					my $uGQ = $uINFO[3];
				if (($aGT ne "./.") and ($uGT ne "./.")) {
				my ($aAle1, $aAle2) = split('/', $aGT);
				my ($uAle1, $uAle2) = split('/', $uGT);
								if ($aGQ >= 20 and $uGQ >= 20) {
									if (($aGT eq $uGT)) {
										delete $DATA1{$var};
									}
								}			
							}				
						}
					}
				}
			}
		
					



@varID =keys(%DATA1);
print scalar(@varID);
print "\n";




open(MYOUTPUTFILE, ">$Output_File");
@IDs =keys(%DATA1);
print MYOUTPUTFILE "Chr\tStart\tEnd\tRef\tAlt\tFunc.refGene\tGene.refGene\tGeneDetail.refGene\tExonicFunc.refGene\tAAChange.refGene\tExAC\tCG69\tESP6500si_ALL\t1000g2014sep_all\tsnp141\tsnp138NonFlagged\tclinvar_20140929\tOMIM_Regions\tpubsBlat\tgerp++gt2\tgwasCatalog\tgenomicSuperDups\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t450\t452\t925\t926\t927\t929\n"; 
foreach (@IDs) {
	if ($DATA1{$_}{"FORMAT"}) { 
print MYOUTPUTFILE $DATA1{$_}{"Chr"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"Start"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"End"}."\t"; 
print MYOUTPUTFILE $DATA1{$_}{"Ref"}."\t"; 
print MYOUTPUTFILE $DATA1{$_}{"Alt"}."\t"; 
print MYOUTPUTFILE $DATA1{$_}{"Func.refGene"}."\t"; 
print MYOUTPUTFILE $DATA1{$_}{"Gene.refGene"}."\t"; 
print MYOUTPUTFILE $DATA1{$_}{"GeneDetail.refGene"}."\t"; 
print MYOUTPUTFILE $DATA1{$_}{"ExonicFunc.refGene"}."\t"; 
print MYOUTPUTFILE $DATA1{$_}{"AAChange.refGene"}."\t"; 
print MYOUTPUTFILE $DATA1{$_}{"exac02"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"CG69"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"ESP6500si_ALL"}."\t";  
print MYOUTPUTFILE $DATA1{$_}{"1000g2014sep_all"}."\t"; 
print MYOUTPUTFILE $DATA1{$_}{"snp141"}."\t"; 
print MYOUTPUTFILE $DATA1{$_}{"snp138NonFlagged"}."\t"; 
print MYOUTPUTFILE $DATA1{$_}{"clinvar_20140929"}."\t"; 
print MYOUTPUTFILE $DATA1{$_}{"OMIM_Regions"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"pubsBlat"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"gerp++gt2"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"gwasCatalog"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"genomicSuperDups"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"CHROM"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"POS"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"ID"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"REF"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"ALT"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"QUAL"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"FILTER"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"INFO"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"FORMAT"}."\t"; 
print MYOUTPUTFILE $DATA1{$_}{"450"}."\t"; 
print MYOUTPUTFILE $DATA1{$_}{"452"}."\t"; 
print MYOUTPUTFILE $DATA1{$_}{"925"}."\t"; 
print MYOUTPUTFILE $DATA1{$_}{"926"}."\t"; 
print MYOUTPUTFILE $DATA1{$_}{"927"}."\t";
print MYOUTPUTFILE $DATA1{$_}{"929"}."\n"; 
  		}
  	};		
close(MYOUTPUTFILE);

 
