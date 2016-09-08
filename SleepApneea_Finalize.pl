#!/usr/bin/perl -w
use strict;
use warnings;

 my $annovar_data	= "./SleepApnea.Exome_hg19.vcf.Filtered.avinput.genome_summary.txt";
 my $cg69 	= "./SleepApnea.Exome_hg19.vcf.Filtered.avinput.hg19_cg69_dropped";
 my esp6500si	= "./SleepApnea.Exome_hg19.vcf.Filtered.avinput.hg19_esp6500si_all_dropped";
 my $gwas 	= "./SleepApnea.Exome_hg19.vcf.Filtered.avinput.hg19_gwasCatalog";
 my $tfbs	= "./SleepApnea.Exome_hg19.vcf.Filtered.avinput.hg19_tfbsConsSites";
 my $phast46	= "./SleepApnea.Exome_hg19.vcf.Filtered.avinput.hg19_phastCons46way";
 my $kggseq_Data	= "./SleepApnea.kggseq.txt";
 my ($CHROM, $POS, $ID, %Annotation, @varID, $var, );

open (LIST, $annovar_data) || die "File not found\n";     
     while (<LIST>) {
         my @annovar = split(/\t/, $_);
         my $CHROM
         my $POS
         
       
         
         $Annotation{ "$CHROM\t$POS" }{ $UID[0] } = $S1;
         $Annotation{ "$CHROM\t$POS" }{ $UID[1] } = $S2;
         $Annotation{ "$CHROM\t$POS" }{ $UID[2] } = $S3;
         $Annotation{ "$CHROM\t$POS" }{ $UID[3] } = $S4;
         $Annotation{ "$CHROM\t$POS" }{ $UID[4] } = $S5;
         $Annotation{ "$CHROM\t$POS" }{ $UID[5] } = $S6;
         $Annotation{ "$CHROM\t$POS" }{ $UID[6] } = $S7;
         $Annotation{ "$CHROM\t$POS" }{"Depth"} = $depth; #Cumulative depth accross all samples
         $Annotation{ "$CHROM\t$POS" }{"Quality"} = $QUAL;
         $Annotation{ "$CHROM\t$POS" }{"EntireLine"} = "$CHROM\t$POS\t$ID\t$REF\t$ALT\t$QUAL\t$FILTER\t$INFO\t$FORMAT\t$S1\t$S2\t$S3\t$S4\t$S5\t$S6\t$S7";
     };
close(LIST);

`rm $tmpdata`;
@varID =keys(%VarData);
print scalar(@varID);
print "\n";


@varID =keys(%VarData);
foreach $var (@varID) {
	if ($VarData{$var}{"Depth"}<16){
		delete $VarData{$var};
	}
	elsif ($VarData{$var}{"Quality"}<30){
		delete $VarData{$var};
	}
	else {
	my ($GT1, $PL1, $DP1, $GQ1) = split(':', $VarData{$var}{$UID[0]});
	my ($GT2, $PL2, $DP2, $GQ2) = split(':', $VarData{$var}{$UID[1]});
	my ($GT3, $PL3, $DP3, $GQ3) = split(':', $VarData{$var}{$UID[2]});
	my ($GT4, $PL4, $DP4, $GQ4) = split(':', $VarData{$var}{$UID[3]});
	my ($GT5, $PL5, $DP5, $GQ5) = split(':', $VarData{$var}{$UID[4]});
	my ($GT6, $PL6, $DP6, $GQ6) = split(':', $VarData{$var}{$UID[5]});
	my ($GT7, $PL7, $DP7, $GQ7) = split(':', $VarData{$var}{$UID[6]});
	
	if ($GQ1 < 30 and $GQ2 < 30 and $GQ3 < 30 and $GQ4 < 30 and $GQ5 < 30 and $GQ6 < 30 and $GQ7 < 30) {
		delete $VarData{$var};
	}
	if ($DP1 < 8 and $DP2 < 8 and $DP3 < 8 and $DP4 < 8 and $DP5 < 8 and $DP6 < 8 and $DP6 < 7) {
		delete $VarData{$var};
	}
}
}		



# Filter variants by comparison to parents, top down (gen1 then gen2.......)
#sort UID array by value of Key Generation in Hash sample_info



my @Parents =("Mother", "Father");
foreach my $ID (@UID) {
	foreach my $Parent (@Parents) {
		if ($sample_info{$ID}{$Parent} ne "NS") {
			my $pID = $sample_info{$ID}{$Parent};
			if($sample_info{$ID}{"Affected"} eq "Y" and $sample_info{$pID}{"Affected"} eq "Y"){	
			#redefine @varID  
				@varID =keys(%VarData);	
				foreach $var (@varID) {
					my ($GT, $PL, $DP, $GQ) = split(':', $VarData{$var}{$ID});
					my ($pGT, $pPL, $pDP, $pGQ) = split(':', $VarData{$var}{$pID});
					my ($Ale1, $Ale2) = split('/', $GT);
					my ($pAle1, $pAle2) = split('/', $pGT);
					if ($DP >= 8 and $pDP >= 8) {
						if (($Ale1 eq $Ale2) or ($pAle1 eq $pAle2)) {
							delete $VarData{$var};
							}				
						}	
					}
				}	
					
			elsif($sample_info{$ID}{"Affected"} eq "Y" and $sample_info{$pID}{"Affected"} eq "N"){	
			#redefine @varID
				@varID =keys(%VarData);
				foreach $var (@varID) {
					my ($GT, $PL, $DP, $GQ) = split(':', $VarData{$var}{$ID});
					my ($pGT, $pPL, $pDP, $pGQ) = split(':', $VarData{$var}{$pID});
					my ($Ale1, $Ale2) = split('/', $GT);
					my ($pAle1, $pAle2) = split('/', $pGT);
					if ($DP >= 8 and $pDP >= 8) {
						if (($GT eq $pGT) or ($Ale1 eq $Ale2)){
							delete $VarData{$var};
							}
						}	
					}
				}	
		 	
		 	elsif($sample_info{$ID}{"Affected"} eq "N" and $sample_info{$pID}{"Affected"} eq "Y"){	
			#redefine @varID
				@varID =keys(%VarData);
				foreach $var (@varID) {
					my ($GT, $PL, $DP, $GQ) = split(':', $VarData{$var}{$ID});
					my ($pGT, $pPL, $pDP, $pGQ) = split(':', $VarData{$var}{$pID});
					my ($Ale1, $Ale2) = split('/', $GT);
					my ($pAle1, $pAle2) = split('/', $pGT);
					if ($DP >= 8 and $pDP >= 8) {
						if (($GT eq $pGT) or ($pAle1 eq $pAle2)){
							delete $VarData{$var};

						}	
					}		
				}			
			}
		}
	}
}

@varID =keys(%VarData);
print scalar(@varID);
print "\n";

#Define subset of unaffected Individuals
my @Unaffected;
foreach (@UID) {
	if ($sample_info{$_}{"Affected"} eq "N") {
		push(@Unaffected, $_);	
	}
}
print "@Unaffected\n";
# Define subset of affected Individuals
my @Affected;
foreach (@UID) {
	if ($sample_info{$_}{"Affected"} eq "Y") {
		push(@Affected, $_);	
	}
}
print "@Affected\n";

# Filter variants by presence of homozygosity among affected Individuals 

if (scalar(@Affected) > 0) {
			foreach (@Affected) {
			@varID =keys(%VarData);
			foreach $var (@varID) {
			my ($GT, $PL, $DP, $GQ) = split(':', $VarData{$var}{$_});
			my ($Ale1, $Ale2) = split('/', $GT);
			if ($DP >= 8) {
				if ($Ale1 == $Ale2) {
					delete $VarData{$var};
					}
				}			
			}	
		}	          
	}
@varID =keys(%VarData);
print scalar(@varID);
print "\n";

#Filter variants by presence of variant
if (scalar(@Unaffected) > 0) {
			foreach (@Unaffected) {
			@varID =keys(%VarData);
			foreach $var (@varID) {
			my ($GT, $PL, $DP, $GQ) = split(':', $VarData{$var}{$_});
			if ((($VarData{$var}{"Depth"}>=8) and ($GT eq "0/1")) or (($DP>=8) and ($GT eq "1/1"))){
					delete $VarData{$var};
					}
				}			
			}	
		}	          
	
@varID =keys(%VarData);
print scalar(@varID);
print "\n";



# Filter variants by pairwise comparison of each affected to unaffected subjects
foreach my $aftd (@Affected) {
	foreach my $unaftd (@Unaffected) {
		@varID =keys(%VarData);
		foreach $var (@varID) {
			my ($aGT, $aPL, $aDP, $aGQ) = split(':', $VarData{$var}{$aftd});
			my ($uGT, $uPL, $uDP, $uGQ) = split(':', $VarData{$var}{$unaftd});
			if ($aDP >= 8 and $uDP >= 8) {
			if (($aGT eq $uGT)) {
			delete $VarData{$var};
						}
					}			
				}				
			}
		}
		
		
my @FinalVariantPositions =keys(%VarData);
print scalar(@FinalVariantPositions);
print "\n";


open(MYOUTPUTFILE, ">$InputVCF.tmp.vcf");
foreach (@FinalVariantPositions) { 
	print MYOUTPUTFILE $VarData{$_}{"EntireLine"};
	#print MYOUTPUTFILE "\n";
  }		
close(MYOUTPUTFILE);

`cat $Header $InputVCF.tmp.vcf > $InputVCF.KindredFilteredFinal.vcf`;
`rm $InputVCF.tmp.vcf`; 
