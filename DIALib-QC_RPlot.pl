#!/usr/bin/perl -w

###############################################################################
# $Id: DIALib-QC_RPlot.pl 5692 2020-04-24 22:00:00Z cmidha $
#
# SBEAMS is Copyright (C) 2000-2020 Institute for Systems Biology
# This program is governed by the terms of the GNU General Public License (GPL)
# version 2 as published by the Free Software Foundation.  It is provided
# WITHOUT ANY WARRANTY.  See the full description of GPL terms in the
# LICENSE file distributed with this software.
###############################################################################

###############################################################################
# This script is used to generate plots of spectral ion library characterstics
# after assessed by DIALib-QC.1.0 tool.
###############################################################################


use strict;
#die "Usage: array_it.pl library.fullstats" unless $ARGV[0];
#die "Usage: DIALib-QC_RPlot.pl library.fullstats library.RT" unless @ARGV;
die "Usage: DIALib-QC_RPlot.pl library.fullstats library.RT" || "Usage: DIALib-QC_RPlot.pl library.fullstats" unless (scalar@ARGV !=0);

system("rm -rf DIALibQC.log");
#Loading File
my $filename = $ARGV[0]."pdf";
#my $filename = $ARGV[0]."png"; ### uncomment if want png file
$filename =~ s/fullstats//;
print ($filename," is generating\n");

open(LOG, ">>", "DIALibQC.log");

print LOG $filename," is generating\n";

#Installing R packages to run R code

use Statistics::R;
my $R = Statistics::R->new();
$R->run( q`library(ggplot2)` );
$R->run( q`library(scales)` );
$R->run( q`install.packages("ggpubr")`);
$R->run( q`library(ggpubr)`);

$R->set( 'filename', $filename );
$R->run(q`pdf(file=filename, width=16, height=16)`);

#$R->run(q`png(file=filename, width = 465, height =465 ,  units='mm', res = 300)`); ### uncomment if want png file


##3 Plot C if RT file provided

if (-e $ARGV[1]){
	
	open FIL2, $ARGV[1];
	my @RTfile= <FIL2>;
	my @rt2_value = ();
	my @rt3_value = (); 

	foreach $_(@RTfile){
		chomp($_);
		print LOG ($_);
		my ( $rt2, $rt3 ) = split( '\t', $_ );
		chomp($rt3);
		print LOG  "RT $rt2\n";
		push @rt2_value, $rt2;
		push @rt3_value, $rt3;
	}

	foreach $_(@rt2_value){print LOG $_;}
	foreach $_(@rt3_value){print LOG $_;}
	
	my $rt_array_size = scalar(@rt2_value);
	print LOG "RT ARRAY SIZE $rt_array_size\n";
	$R->set( 'rt2_value', \@rt2_value );
	$R->set( 'rt3_value', \@rt3_value );
	$R->run( q`dfrt<- data.frame(RT2value = c(rt2_value) , RT3value = c(rt3_value))` );
	$R->run( q`correlation = cor(dfrt$RT2value , dfrt$RT3value, method = "pearson")`);
	$R->run( q`correlation = sprintf("%1f",correlation)`);
	my $correlationval = $R->get('correlation');
	print LOG "CORRELATIONVAL",$correlationval;
	$R->run(q`N=NROW(dfrt)`);
	
	my $nrowcount = $R->get('N');
 	print LOG "NROWCOUNT ",$nrowcount;
	
	$R->run(q`title = paste0("Library +2/+3 pair RT correlation, n=",N, ", R2=",correlation)`);
	
	my $title = $R->get('title');
	print LOG "TITLE ",$title;
	$R->run(q`p3<-ggplot(dfrt, aes(x = RT2value, y = RT3value))+ 
	geom_point(shape = 1, size =4, color = '#008080')+
	theme_bw()+ labs(title =title, x = "+2RT (mins)", y = "+3RT (mins)", size = 40)+
	geom_smooth( method = "lm", color = 'black',alpha=0, size = .3  )`);
	close FIL2;
	}

open FIL, $ARGV[0] || die "library.fullstats is minimum file required to generate plots";
while ( my $line=<FIL> ) {
	chomp $line;
	my ( $key, $data ) = split( ':', $line );
	#open OUT, ">$key.csv"; ##UNCOMMENT ME TO GENERATE INDIVIDUAL CSVS
	my $i = 0;  
	my $j = 0; 
	my @item_array = ();
	my @val_array = (); 
	for my $itemvalue ( split( ',', $data ) ) {
		my ( $item, $val ) = split( '=', $itemvalue );
		push @item_array, $item; 
		push @val_array, $val;
		#print OUT join( ',', $item, $val ) . "\n";  ##UNCOMMENT ME TO GENERATE INDIVIDUAL CSVS

	}
  
	$R->set( 'Items', \@item_array );
	$R->set( 'Values', \@val_array );

	##p1 A PLOT
	if ($key =~ m/PrecusorCounts/){
		foreach $_(@val_array) {print LOG $_," precursor count\n";}
		$R->run( q`totalfragions<-sum(Values)` );
		$R->run(q` Values = round(Values*100/totalfragions, 2)`);
		$R->run(q`df<-data.frame(Parameters_type=c(Items),Frequency=c(Values))`);
		$R->run( q`minval<- min(Values)`);
		$R->run( q`maxval<- max(Values)`);
		$R->run(q`p1<-ggplot(df, aes(Parameters_type, Frequency))+ xlim(350,1300)+ 
		theme_bw()+ geom_col(fill="#008080")+ xlab("Precusor m/z")+
		ylab("Frequency in percent (%)") + ggtitle("Precusor m/z")`); 
	}
	
	##p2 B POLT
	if ($key =~ m/^pre_z/){

		$R->run( q`totalfragions<-sum(Values)` );
		$R->run(q` Values = round(Values*100/totalfragions,2)`);
	        $R->run(q`Items_char= as.character(Items)`);
		$R->run(q`df<-data.frame(Parameters_type=c(Items_char),Frequency=c(Values))`);
		$R->run(q`p2<-ggplot(df, aes(Parameters_type, Frequency))+ theme_bw()+
		geom_col(fill="#008080")+ xlab("Precursor charge")+
		geom_text(aes(label = Frequency), size = 6, hjust = 0.5, vjust = 0, position = "stack", colour = 'black')+
		ylab("Precursors in percent (%)")+ggtitle("Precursor charge")`);
		@item_array = @val_array = ();	
	}


	##p4 C/D PLOT
	if ($key =~ m/stripped_len/){      
		$R->run( q`totalfragions<-sum(Values)` );
		$R->run(q` Values = round(Values*100/totalfragions, 2)`);
		$R->run(q`df<-data.frame(Parameters_type=c(Items),Frequency=c(Values))`);
		$R->run(q`p4<-ggplot(df, aes(Parameters_type, Frequency))+ theme_bw()+
		geom_col(fill="#008080")+ xlab("Peptide length") + ylab("Peptides in percent (%)")+
		ggtitle("Peptide length")`); 
	}
	
	##p5 E PLOT
	if ($key =~ m/mods/){
		my $null_index = 0;
		for (my $i = 0; $i <@item_array; $i++){
			if ($item_array[$i] =~ m/none/){
				$null_index = $i;
			}
		}

                my $any_index = 0;
         	for (my $i = 0; $i <@item_array; $i++){
                        if ($item_array[$i] =~ m/any/){
                                $any_index = $i;
                        }
                }

		print LOG $null_index,"\t",$val_array[$null_index],"\n" ;
		splice @item_array, $null_index, 1;
		splice @val_array, $null_index, 1;


                splice @item_array, $any_index, 1;
                splice @val_array, $any_index, 1;

		$R->set( 'Items', \@item_array );
		$R->set( 'Values', \@val_array );

		$R->run( q`df<-data.frame(Parameters_type=c(Items),Frequency=c(Values))`);
		$R->run(q`p5<-ggplot(df, aes(Parameters_type, Frequency))+theme_bw()+
		geom_col(fill="#008080")+ xlab("Modification type")+
		geom_text(aes(label = Frequency), size = 6, hjust = 0.5, vjust = 0, position = "stack", colour = 'black')+
		ylab("Number of modifications")+ggtitle("Modifications")`); 
	}
	
	##p6 F PLOT
	if ($key =~ m/PepsPerProtein/){
		my @first4value = ();
		my @first4item = ();
		my $val4abovesum = 0;
		if (scalar(@item_array) <5){
			$R->run(q`Items_char= as.character(Items)`);
			$R->run(q`df<-data.frame(Parameters_type=c(Items_char),Frequency=c(Values))`);
			$R->run(q`p6<-ggplot(df, aes(Parameters_type, Frequency))+ 
			theme_bw()+ scale_x_discrete("Peptides per protein", limits =c(Items_char))+
			geom_text(aes(label = Frequency), size = 6, hjust = 0.5, vjust = 0, position = "stack", colour = 'black')+
			geom_col(fill="#008080") + ylab("Number of protiens")+
			ggtitle("Peptides per protein") `);
		}else{
			for (my $i= 0; $i<4; $i++){
				push @first4value, $val_array[$i];
				push @first4item, $item_array[$i];
			}

			for (my $i = 4; $i< scalar@val_array; $i++){
				$val4abovesum = $val4abovesum + $val_array[$i];
			}
			
			print LOG ($val4abovesum);
			$first4value[4]=$val4abovesum;
			$R->set( 'Items', \@first4item );
			$R->run(q`Items_char= as.character(Items)`);
			$R->run(q`Items_char= c(Items_char, ">5")`);
			$R->set( 'Values', \@first4value );
			$R->run(q`df<-data.frame(Parameters_type=c(Items_char),Frequency=c(Values))`);
			$R->run(q`p6<-ggplot(df, aes(Parameters_type, Frequency))+ 
			theme_bw()+ scale_x_discrete("Peptides per protein", limits =c(Items_char))+
			geom_text(aes(label = Frequency), size = 6, hjust = 0.5, vjust = 0, position = "stack", colour = 'black')+
			geom_col(fill="#008080") + ylab("Number of proteins")+
			ggtitle("Peptides per protein")`);
		}
	}	
	
	
	##p7 G POLT
	if ($key =~ m/FragmentsPerPrecursor/){
		$R->run( q`totalfragions<-sum(Values)` );
		$R->run(q` Values = round(Values*100/totalfragions, 2)`);
		$R->run(q`Items_char= as.character(Items)`);
		$R->run(q`df<-data.frame(Parameters_type=c(Items_char),Frequency=c(Values))`);
	
		print LOG scalar(@item_array),"\n";
		$R->run(q`p7<-ggplot(df, aes(Parameters_type, Frequency))+ theme_bw()+
		geom_col(fill="#008080")+ scale_x_discrete("Fragments per precursor", limits =Items_char)+
		geom_text(aes(label = Frequency), size = 3, hjust = 0.5, vjust = 0, position = "stack", colour = 'black')+	
		ylab("Frequency in percent (%)") + ggtitle("Fragments per precursor")`);

	}
	
	##p8 H PLOT 
	if ($key =~ m/frg_s/){
		$R->run( q`totalfragions<-sum(Values)` );
		$R->run(q` Values = signif(Values*100/totalfragions,digits = 1)`);
		$R->run( q`df<- data.frame(Parameters_type = c(Items) , Frequency = c(Values))` );
		$R->run(q`p8<-ggplot(df, aes(Parameters_type, Frequency))+ theme_bw()+
		geom_col(fill="#008080")+ xlab("Fragment ion type") + ylab("Frequency in percent (%)")+
		geom_text(aes(label = Frequency), size = 6, hjust = 0.5, vjust = 0, position = "stack", colour = 'black')+
		ggtitle("Fragment ion")`); 
	}
	
	
	##p9 I PLOT
	if ($key =~ m/frg_z/){

		foreach $_(@val_array){print LOG $_;}
		$R->run( q`totalfragions<-sum(Values)` );
		my $totalfragions = $R->get('totalfragions');
		print LOG "totalfragions, $totalfragions";
		$R->run(q`Values = round(Values*100/totalfragions,2)`);
		$R->run(q`Items_char= as.character(Items)`);
		$R->run(q`df<-data.frame(Parameters_type=c(Items_char),Frequency=c(Values))`);
		$R->run(q`p9<- ggplot(df, aes(Parameters_type, Frequency))+
		xlab("Fragment ion charge")+
		theme_bw()+geom_col(fill="#008080")+
		geom_text(aes(label = Frequency), size = 6, hjust = 0.5, vjust = 0, position = "stack", colour = 'black')+
		scale_x_discrete("Fragment ion charge", limits =c(Items_char))+
		ylab("Fragment ions in percent (%)") + ggtitle("Fragment ion charge") `); 
	}
	###p10
	# if ($key =~ m/frg_n/){
		# $R->run( q`totalfragions<-sum(Values)` );
		# $R->run(q` Values = round(Values*100/totalfragions,2)`);
		# $R->run(q`df<-data.frame(Parameters_type=c(Items),Frequency=c(Values))`);
		# $R->run(q`p10<-ggplot(df, aes(Parameters_type, Frequency))+ theme_bw()+
		# geom_col(fill="#008080")+ xlab("Fragment ion") + ylab("Precursor in percent(%)")+
		# ggtitle("Number of fragment ions per precursor")`); 
	# }
	###p11
	# if ($key =~ m/Distinct_pre_z/){
		# $R->run( q`totalfragions<-sum(Values)` );
		# #$R->run(q` Values = signif(Values*100/totalfragions, digits =1)`);
		# $R->run(q` Values = round((Values*100/totalfragions),0)`);
		# #$R->run(q` Values = round(Values/1000,0)`);
		# $R->run(q`df<-data.frame(Parameters_type=c(Items),Frequency=c(Values))`);
		# $R->run(q`p11<-ggplot(df, aes(Parameters_type, Frequency))+ theme_bw()+
		# geom_col(fill="#008080")+ xlab("Precursor charge")+
		# geom_text(aes(label = Frequency), size = 6, hjust = 0.5, vjust = 0, position = "stack", colour = 'black')+
		# ylab("Precursors in percent(%)")+ggtitle("Precursor charger")`);
	# }
		
	}
	
if (-e $ARGV[1]){
	$R->run(q`ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9 , labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"), ncol = 3, nrow = 3)`);
}else {
	$R->run(q`ggarrange(p1, p2, p4, p5, p6, p7, p8, p9 , labels = c("A", "B", "C", "D", "E", "F", "G", "H"), ncol = 3, nrow = 3)`);
}

$R->run(q`dev.off()`);
	
#close OUT;  ##UNCOMMENT ME TO GENERATE INDIVIDUAL CSVS
close LOG; 
