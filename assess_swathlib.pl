#!/usr/bin/perl
#
#  script to assess ion libraries (DIALib-QC)
#
=for comment

Copyright (C) 2020 Moritz Lab, Institute for Systems Biology

DIALib-QC.1.2 program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.


You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
=cut


# Import libraries
use strict;
use Data::Dumper; 
use File::Basename;
use Getopt::Long;
use Text::ParseWords;

my $t0 = time();
print STDERR "Starting run\n";

# Test for presence of optional (non-core) regression module
my $has_regression = 1;
eval "use Statistics::Regression; 1" or $has_regression = 0;
if ( !$has_regression ) {
   print STDERR "Unable to load module Statistics::Regression, skipping RT correlation analysis\n";
}

# Global vars
my %opts = read_options();  # Command-line options
my %problem_assays;
##+
## Hash of known modifications
##-
my %known_mods = (
		'C[CAM]' => {'mz' => 57.0215, 'aa' => 'c'},
		'C(UniMod:4)' => {'mz' => 57.0215, 'aa' => 'c'},
		'C[+57]' => {'mz' => 57.0215, 'aa' => 'c'},
		'C[+57.0]' => {'mz' => 57.0215, 'aa' => 'c'},
		'C[57]' => {'mz' => 57.0215, 'aa' => 'c'},
		'C[Carbamidomethyl (C)]' => {'mz' => 57.0215, 'aa' => 'c'},
		'C[+C2+H3+N+O]' => {'mz' => 57.0215, 'aa' => 'c'},
		'C[160]' => {'mz' => 57.0215, 'aa' => 'c'},
		'C(UniMod:39)' => {'mz' => 45.987721 , 'aa' => 'c'}, 
		'C[149]' => {'mz' => 45.9877721, 'aa' => 'c'}, 
		'C[46]' => {'mz' => 45.987721, 'aa' => 'c'}, 
		'C[+46.0]' => {'mz' => 45.987721, 'aa' => 'c'}, 
		'C[Methylthio]' => {'mz' => 45.987721 , 'aa' => 'c'}, 
                'C[MSH]' => {'mz' => 45.987721 , 'aa' => 'c'},
		'C[PCm]' => {'mz' => 39.994915, 'aa' => 'c'}, 
		'C[143]' => {'mz' => 39.994915, 'aa' => 'c'}, 
		'C[40]' => {'mz' => 39.994915, 'aa' => 'c'}, 
		'C[+40]' => {'mz' => 39.994915, 'aa' => 'c'}, 
		'C[+40.0]' => {'mz' => 39.994915, 'aa' => 'c'},
		'C(UniMod:35)' => {'mz' => 15.994915, 'aa' => 'c'}, 
		'L(UniMod:35)' => {'mz' => 15.994915, 'aa' => 'l'}, 
		'A(UniMod:35)' => {'mz' => 15.994915, 'aa' => 'a'},
		'D(UniMod:35)' => {'mz' => 15.994915, 'aa' => 'd'},
		'V(UniMod:35)' => {'mz' => 15.994915, 'aa' => 'v'},
		'F(UniMod:35)' => {'mz' => 15.994915, 'aa' => 'f'},
		'Y(UniMod:35)' => {'mz' => 15.994915, 'aa' => 'y'},
		'T(UniMod:35)' => {'mz' => 15.994915, 'aa' => 't'},
		'G(UniMod:35)' => {'mz' => 15.994915, 'aa' => 'g'},
		'Q(UniMod:35)' => {'mz' => 15.994915, 'aa' => 'q'},
		'N(UniMod:35)' => {'mz' => 15.994915, 'aa' => 'n'},
		'I(UniMod:35)' => {'mz' => 15.994915, 'aa' => 'i'},
		'E(UniMod:35)' => {'mz' => 15.994915, 'aa' => 'e'},
		'S(UniMod:35)' => {'mz' => 15.994915, 'aa' => 's'}, 
		'C(UniMod:26)' => {'mz' => 39.994915, 'aa' => 'c'}, 
		'C[Pyro-Carbamidomethyl (Any N-term)]' => {'mz' => -17.026549, 'aa' => 'c'}, 
		'M[CAM]' => {'mz' => 57.0215, 'aa' => 'm'},
		'M(UniMod:4)' => {'mz' => 57.0215, 'aa' => 'm'},
		'M[+57]' => {'mz' => 57.0215, 'aa' => 'm'},
		'M[+57.0]' => {'mz' => 57.0215, 'aa' => 'm'},
		'M[57]' => {'mz' => 57.0215, 'aa' => 'm'},
		'M[Oxi]' => {'mz' => 15.994915, 'aa' => 'm'},
		'M[Oxidation (M)]' => {'mz' => 15.994915, 'aa' => 'm'},
		'M(UniMod:35)' => {'mz' => 15.994915, 'aa' => 'm'},
		'M[147]' => {'mz' => 15.994915, 'aa' => 'm'},
		'M[+16]' => {'mz' => 15.994915, 'aa' => 'm'},
		'M[+O]' => {'mz' => 15.994915, 'aa' => 'm'},
		'M[+16.0]' => {'mz' => 15.994915, 'aa' => 'm'},
		'M[16]' => {'mz' => 15.994915, 'aa' => 'm'},
		'P[Oxi]' => {'mz' => 15.994915, 'aa' => 'p'},
		'P[Oxidation (P)]' => {'mz' => 15.994915, 'aa' => 'p'},
		'W[Oxi]' => {'mz' => 15.994915, 'aa' => 'w'},
		'W[Oxidation (W)]' => {'mz' => 15.994915, 'aa' => 'w'},
		'W(UniMod:35)' => {'mz' => 15.994915, 'aa' => 'w'},
		'W[+16]' => {'mz' => 15.994915, 'aa' => 'w'},
		'W[+16.0]' => {'mz' => 15.994915, 'aa' => 'w'},
		'W[16]' => {'mz' => 15.994915, 'aa' => 'w'},
		'W[220]' => {'mz' => 15.9949, 'aa' => 'w'},
		'Q[Dea]' => {'mz' => 0.984016, 'aa' => 'q'},
		'R[Dea]' => {'mz' => 0.984016, 'aa' => 'r'},
		'T[Dhy]' => {'mz' => -18.010565, 'aa' => 't'},
		'S[Pho]' => {'mz' => 79.966331, 'aa' => 's'},
		'[-1R]' => {'mz' => -156.101111, 'aa' => 'x'},
		'N[dAm]' => {'mz' => -17.026549, 'aa' => 'n'},
		'[CRM]-' => {'mz' => 43.005814, 'aa' => 'x'},
		'D[Sud]' => {'mz' => 31.972071, 'aa' => 'd'},
		'S[UGG]' => {'mz' => 114.042927, 'aa' => 's'},
		'Y[Amn]' => {'mz' => 15.010899, 'aa' => 'y'},
		'D[Oxi]' => {'mz' => 15.994915, 'aa' => 'd'},
		'[-1K]' => {'mz' => -128.094963, 'aa' => 'x'},
		'S[PPE]' => {'mz' => 340.085794, 'aa' => 's'},
		'S[HNc]' => {'mz' => 203.079373, 'aa' => 's'},
		'[AAS]-' => {'mz' => 26.01565, 'aa' => 'x'},
		'M[CRM]' => {'mz' => 43.005814, 'aa' => 'm'},
		'F[1Br]' => {'mz' => 77.910511, 'aa' => 'f'},
		'S[1Me]' => {'mz' => 14.01565, 'aa' => 's'},
		'N[1NF]' => {'mz' => 349.137281, 'aa' => 'n'},
		'K[HNE]' => {'mz' => 156.11503, 'aa' => 'k'},
		'D[Cox]' => {'mz' => 43.989829, 'aa' => 'd'},
		'K[Oxi]' => {'mz' => 15.994915, 'aa' => 'k'},
		'T[HNc]' => {'mz' => 203.079373, 'aa' => 't'},
		'E[GPE]' => {'mz' => 197.04531, 'aa' => 'e'},
		'N[1K2]' => {'mz' => 730.264392, 'aa' => 'n'},
		'S[-2H]' => {'mz' => -2.01565, 'aa' => 's'},
		'[CuX]' => {'mz' => 61.921774, 'aa' => 'x'},
		'N[Oxi]' => {'mz' => 15.994915, 'aa' => 'n'},
		'R[Oxi]' => {'mz' => 15.994915, 'aa' => 'r'},
		'Y[Oxi]' => {'mz' => 15.994915, 'aa' => 'y'},
		'F[Oxi]' => {'mz' => 15.994915, 'aa' => 'f'},
		'D[CuX]' => {'mz' => 61.921774, 'aa' => 'd'},
		'K[3Me]' => {'mz' => 42.04695, 'aa' => 'k'},
		'D[1Me]' => {'mz' => 14.01565, 'aa' => 'd'},
		'L[1Me]' => {'mz' => 14.01565, 'aa' => 'l'},
		'R[1Me]' => {'mz' => 14.01565, 'aa' => 'r'},
		'E[1Me]' => {'mz' => 14.01565, 'aa' => 'e'},
		'H[1Me]' => {'mz' => 14.01565, 'aa' => 'h'},
		'Q[1Me]' => {'mz' => 14.01565, 'aa' => 'q'},
		'R[2Me]' => {'mz' => 28.0313, 'aa' => 'r'},
		'P[2Ox]' => {'mz' => 31.989829, 'aa' => 'p'},
		'Y[2Ox]' => {'mz' => 31.989829, 'aa' => 'y'},
		'H[CAM]' => {'mz' => 57.0215, 'aa' => 'h'},
		'Y[Dhy]' => {'mz' => -18.010565, 'aa' => 'y'},
		'N[K2F]' => {'mz' => 876.322301, 'aa' => 'n'},
		'Y[-2H]' => {'mz' => -2.01565, 'aa' => 'y'},
		'D[MSH]' => {'mz' => 45.987721, 'aa' => 'd'},
		'T[NHN]' => {'mz' => 947.323029, 'aa' => 't'},
		'S[GCn]' => {'mz' => 176.032088, 'aa' => 's'},
		'K[OH2]' => {'mz' => 340.100562, 'aa' => 'k'},
		'N[G2d]' => {'mz' => 1216.422863, 'aa' => 'n'},
		'T[1Ac]' => {'mz' => 42.010565, 'aa' => 't'},
		'S[1Ac]' => {'mz' => 42.010565, 'aa' => 's'},
		'[Myr]-' => {'mz' => 210.198366, 'aa' => 'x'},
		# 'Y[Dhy]' => {'mz' => -18.010565, 'aa' => 'y'},#EXAMPLE LINE TO ADD MODIFICATION
		'Q(UniMod:7)' => {'mz' => 0.984016, 'aa' => 'q'}, 
		'H(UniMod:35)' => {'mz' => 15.994915, 'aa' => 'h'},
		# 'U' => {'mz' => 150.95363, 'aa' => 'u'}, #U = selenocysteine 
		'[CAM]-' => {'mz' => 57.0215, 'aa' => 'x'},
		'D[CAM]' => {'mz' => 57.0215, 'aa' => 'd'},
		'S[CAM]' => {'mz' => 57.0215, 'aa' => 's'},
		'E[CAM]' => {'mz' => 57.0215, 'aa' => 'e'},
		'T[CAM]' => {'mz' => 57.0215, 'aa' => 't'},
		'Y[CAM]' => {'mz' => 57.0215, 'aa' => 'y'},
		'K[CAM]' => {'mz' => 57.0215, 'aa' => 'k'},
		'D[NaX]' => {'mz' => 21.981943, 'aa' => 'd'},
		'E[NaX]' => {'mz' => 21.981943, 'aa' => 'e'},
		'E[KXX]' => {'mz' => 37.955882, 'aa' => 'e'},
		'[Ami]' => {'mz' => -0.984016, 'aa' => 'x'},
		'K[Frm]' => {'mz' => 27.994915, 'aa' => 'k'},
		'[Frm]-' => {'mz' => 27.994915, 'aa' => 'x'},
		'F[2Ox]' => {'mz' => 31.989829, 'aa' => 'f'},
		'M[2Ox]' => {'mz' => 31.989829, 'aa' => 'm'},
		'W[2Ox]' => {'mz' => 31.989829, 'aa' => 'w'},
		'S[Dhy]' => {'mz' => -18.010565, 'aa' => 's'},
		'[+1R]-' => {'mz' => 156.101111, 'aa' => 'x'},
		'T[-2H]' => {'mz' => -2.01565, 'aa' => 't'},
		'[-2H]' => {'mz' => -2.01565, 'aa' => 'x'},
		'D[KXX]' => {'mz' => 37.955882, 'aa' => 'd'},
		'D[Dhy]' => {'mz' => -18.010565, 'aa' => 'd'},
		'E[Dhy]' => {'mz' => -18.010565, 'aa' => 'e'},
		'K[AAR]' => {'mz' => 28.0313, 'aa' => 'k'},
		'W[Kyn]' => {'mz' => 3.994915, 'aa' => 'w'},
		'E[CuX]' => {'mz' => 61.921774, 'aa' => 'e'},
		'[+1K]-' => {'mz' => 128.094963, 'aa' => 'x'},
		'M[DTM]' => {'mz' => -48.003371, 'aa' => 'm'},
		'R[CRM]' => {'mz' => 43.005814, 'aa' => 'r'},
		'H[Oxi]' => {'mz' => 15.994915, 'aa' => 'h'},
		'K[AAS]' => {'mz' => 26.01565, 'aa' => 'k'}, 
		'H[AAS]' => {'mz' => 26.01565, 'aa' => 'h'}, 
		'N[Dea]' => {'mz' => 0.984016, 'aa' => 'n'},
		'W[HKy]' => {'mz' => 19.989829, 'aa' => 'w'},
		'D[2CM]' => {'mz' => 114.0429274, 'aa' => 'd'},
		'R[AGA]' => {'mz' => -43.053433, 'aa' => 'r'},
		'R[Orn]' => {'mz' => -42.021798, 'aa' => 'r'},
		'[NaX]' => {'mz' => 21.981943, 'aa' => 'x'},
		'P[PYD]' => {'mz' => -30.010565, 'aa' => 'p'},
		'N[MDe]' => {'mz' => 14.999666, 'aa' => 'n'},
		'N[Deamidation (N)]' => {'mz' => 0.984016, 'aa' => 'n'},
		'N[Deamidation (NQ)]' => {'mz' => 0.984016, 'aa' => 'n'},
		'N(UniMod:7)' => {'mz' => 0.984016, 'aa' => 'n'},
		'N[115]' => {'mz' => 0.984016, 'aa' => 'n'},
		'N[+1.0]' => {'mz' => 0.984016, 'aa' => 'n'},
		'N[1]' => {'mz' => 0.984016, 'aa' => 'n'}, 
		'E[PGE]' => {'mz' => -18.010565, 'aa' => 'e'},
		'[PGE]-' => {'mz' => -18.010565, 'aa' => 'x'},
		'P[PGP]' => {'mz' => 13.979265, 'aa' => 'p'},
		'[AAR]-' => {'mz' => 28.0313, 'aa' => 'x'},
		'K[CRM]' => {'mz' => 43.005814, 'aa' => 'k'},
		'E[111]' => {'mz' => -18.010565, 'aa' => 'e'},
		'E(UniMod:27)' => {'mz' => -18.010565, 'aa' => 'e'}, 
		'E[-18]' => {'mz' => -18.010565, 'aa' => 'e'},
		'E[-18.0]' => {'mz' => -18.010565, 'aa' => 'e'},
		'E[Glu->pyro-Glu (Any N-term)]' => {'mz' => -18.010565, 'aa' => 'e'},
		'E[Glu->pyro-Glu]' => {'mz' => -18.010565, 'aa' => 'e'},
		'Q[PGQ]' => {'mz' => -17.026549, 'aa' => 'q'},
		'[PGQ]-' => {'mz' => -17.026549, 'aa' => 'x'},
		'Q[111]' => {'mz' => -17.026549, 'aa' => 'q'},
		'Q(UniMod:28)' => {'mz' => -17.026549, 'aa' => 'q'}, 
		'(UniMod:28)' => {'mz' => -17.026549, 'aa' => 'x'}, 
		'Q[-17]' => {'mz' => -17.026549, 'aa' => 'q'},
		'Q[-17.0]' => {'mz' => -17.026549, 'aa' => 'q'},
		'Q[Glu pyro-Glu (Any N-term)]' => {'mz' => -17.026549, 'aa' => 'q'}, 
		'Q[Glu pyro-Glu]' => {'mz' => -17.026549, 'aa' => 'q'}, 
		'Q[Dea]' => {'mz' => 0.984016, 'aa' => 'q'}, 
		'Q[Deamidation (NQ)]' => {'mz' => 0.984016, 'aa' => 'q'},
		'[42]' => {'mz' => 42.010565, 'aa' => 'a'},
		'[42.0]' => {'mz' => 42.010565, 'aa' => 'a'},
		'[1Ac]-' => {'mz' => 42.010565, 'aa' => 'x'},
		'[Acetyl +(N-term)]' => {'mz' => 42.010565, 'aa' => 'a'},
		'[Acetyl +(Any N-term)]' => {'mz' => 42.010565, 'aa' => 'a'},
		'[Acetyl +(Protein N-term)]' => {'mz' => 42.010565, 'aa' => 'a'},
		'M[Acetyl +(Protein N-term)]' => {'mz' => 42.010565, 'aa' => 'm'},
		'K[Acetyl (K)]' => {'mz' => 42.010565, 'aa' => 'k'},

		'A[Acetyl (Protein N-term)]' => {'mz' => 42.010565, 'aa' => 'a'},
		'S[Acetyl (Protein N-term)]' => {'mz' => 42.010565, 'aa' => 's'},
		'D[Acetyl (Protein N-term)]' => {'mz' => 42.010565, 'aa' => 'd'},
		'M[Acetyl (Protein N-term)]' => {'mz' => 42.010565, 'aa' => 'm'},
		'T[Acetyl (Protein N-term)]' => {'mz' => 42.010565, 'aa' => 't'},
		'V[Acetyl (Protein N-term)]' => {'mz' => 42.010565, 'aa' => 'v'},
		'E[Acetyl (Protein N-term)]' => {'mz' => 42.010565, 'aa' => 'e'},
		'G[Acetyl (Protein N-term)]' => {'mz' => 42.010565, 'aa' => 'g'},
		'P[Acetyl (Protein N-term)]' => {'mz' => 42.010565, 'aa' => 'p'},
		'K[Acetyl (Protein N-term)]' => {'mz' => 42.010565, 'aa' => 'k'},
              );
my %mods_current;
my $nterm_mod;
my $cterm_mod;
my $LTOS;	#LossTypeOpenSwath;

# Variables for SWATH bin analysis
my %swaths;
my %swadefs;
my @swaths;

# Constants
use constant WATER_MASS => 18.010565;
use constant H_MASS => 1.00727646688;

if ( $opts{swath_file} ) {
  open SWA, $opts{swath_file};
  while ( my $line = <SWA> ) {
    chomp $line;
    my @line = split( /\s+/, $line );
    next unless $line[0] =~ /\d/; 
    next unless $line[1] =~ /\d/; 
    $swadefs{$line[0]} = $line[1];
  }
  close SWA;
}
if ( $opts{swath_file} ) {
  @swaths = ( qw( swa_defined swa_missing swa_conflict swa_ok swa_conflict_assay swa_5 swa_25 ) );
}

my %stats = ( mc => 0 ); # Hash of output values
# Initialize stats hash
for my $key ( qw( totpeps libpeps speps mc libprots sprots mprots precursor_ok precursor_bad fragment_ok fragment_bad fragment_na ), @swaths ) {
  $stats{$key} = 0;
}

# ppeps is a peptide/protein mapping file from reference DB; pr are proteotypic
my %ppeps = ( prots => {}, peps => {}, pr => {} );

# read in opional peptide digest file
if ( $opts{peptide_file} ) {
  print STDERR "File $opts{peptide_file} does not exist\n" if ! -e $opts{peptide_file};

  open PEPS, $opts{peptide_file}|| die;
  while ( my $line = <PEPS> ) {
    chomp $line;
    my @line = split( /\t/, $line );
    $ppeps{peps}->{$line[1]}++;
    if ( $line[0] !~ /,/ ) {
      $ppeps{pr}->{$line[0]}++;
    }
    for my $prot ( split( /,/, $line[0] ) ) {
      $ppeps{prots}->{$prot}++;
    }
  }
  close PEPS;
}
$stats{libpeps} = scalar( keys( %{$ppeps{peps}} ) );
$stats{libprots} = scalar( keys( %{$ppeps{prots}} ) );
$stats{libprprots} = scalar( keys( %{$ppeps{pr}} ) ); # proteins with at least one proteotypic peptide

my $libname = basename( $opts{ion_library} );
$libname =~ s/\.tsv$//;
$libname =~ s/\.csv$//;

# Open and read ion library file
open ILIB, $opts{ion_library}|| die;
my $head = 1;
my $type;

my $dta_dir = $libname . '_dta';
$dta_dir =~ s/\./_/g;
my $write_data = $opts{write} || 0;
if ( $write_data ) {
  if ( -e "$libname.mgf" ) {
    print "MGF file exists, exiting\n";
    exit;
  }
  open MGF, ">$libname.mgf" || die;
}

my $curr;
my %dta = ( frg => {}, pre => {z => '', mz => ''} );
my $pepkey;

# hash of all accessions recorded in ion library
my %prots;

my %peps = ( sing => {}, 
             mult => {},
             pre_z => {}, # precursor charge
             frg_z => {}, # fragment charge
             frg_n => {}, # fragment number
             frg_s => {}, # tally fragment series
             frg_i => { cnt => 0, inten => 0 }, # Summed intensity of top 5 
             ); 

my %rt = ( mseqs => {}, iRT => {} );
my $irt_cnt = 0;
my %is_irt;
for my $irt ( qw( ADVTPADFSEWSK DGLDAASYYAPVR GAGSSEPVTGLDAK GTFIIDPAAVIR GTFIIDPGGVIR LFLQFGAQGSPFLK LGGNEQVTR TPVISGGPYEYR TPVITGAPYEYR VEATFGVDESNAK YILAGVENSK ) ) {
  $is_irt{$irt}++;
}

# Hash of recognized neutral losses to their masses.
my %losses = ( 'H2O' => 18.010565,
               '1(+C+H4+O+S)' => 64.107989,	
               'NH3' => 17.026549,
               'noloss' => 0,
			   '-17' => 17.026549,
               '-18' => 18.010565,
	           '-64' => 64.107989
			); 

my %decoy = ( mixed => 0, decoy => 0, fwd => 0);
my %extrema = ( precursor_min => 10000, fragment_min => 10000, precursor_max => 0, fragment_max => 0 );

if ( $opts{print_mass_error} ) {
  my $me = $opts{output_file} . ".mass_error";
  if ( -e $me ) { 
    print STDERR "mass error file $me exists, cannot overwrite\n";
    exit;
  }
  open ME, ">$me";
}


my %cmap;            # Map of column headings to number (nth column)
my %idx;             # Map of target information to column index
my %massmap;  
my %precursorstats;  
my %ion_cnt;
my %tmax;
my $osw_altmod = 0;
my $is_CSV;
my $quote_char = '';
my @imarray;
my @imsort;
while ( my $line = <ILIB> ) {
  chomp $line;
  my @line = parse_line( $line );
  if ( !defined $is_CSV ) {
    if ( scalar( @line ) < 5 ) {
      @line = parse_csv( $line );
      if ( scalar( @line ) < 5 ) {
        print STDERR "Unknown file format, must be either TSV or CSV\n";
        exit 256;
      } else {
        $is_CSV++;
        print STDERR "CSV format detected\n";
        if ( $line =~ /^(["']).*["']\s*$/ ) {
          $quote_char = $1;
          print STDERR "Quote character is $quote_char\n";
        }
      }
    } else {
      $is_CSV = 0;
    }
  }

  if ( $head ) {
    my $idx = 0;
    for my $col ( @line ) {
      $cmap{$col} = $idx++;
    }
    if ( !$line[1] ) {
      print STDERR "Error: $ARGV[0] format invalid \n";
    } elsif ( $line =~ /[Q1 q1]/ && ($line =~ /modification_sequence/ || $line =~ /peptideModSeq/)) {
      $type = 'pv';
      print STDERR "Peakview library detected\n";
    } elsif ( $line =~ /PrecursorMz/ ) {
      if ( $line =~ /transition_group_id/ ) {
        print STDERR "OpenSWATH library detected\n";
        $type = 'os';
      } elsif ( $line =~ /iRT Value/ || $line =~ /ReferenceRun/ ) {
        print STDERR "Spectronaut library detected\n";
        $type = 'sn';
      } elsif ( $line =~ /FragmentLossType/ ) {
        print STDERR "Prosit/Fusion library format detected\n";
        $type = 'pf'; # 'Generic'
      } else {
       print STDERR "Error: $ARGV[0] format unknown\n";
      }
      if ( $type eq 'os' ) {
        if ( !defined $cmap{FullUniModPeptideName} ) {
          $cmap{FullUniModPeptideName} = $cmap{FullPeptideName};
        }
        if ( !defined $cmap{PrecursorCharge} ) {
          $cmap{PrecursorCharge} = $cmap{Charge};
        }
        if ( !defined $cmap{Tr_recalibrated} ) {
          $cmap{Tr_recalibrated} = $cmap{RetentionTime};
        }
        die "Can't find peptide name" if !defined $cmap{FullUniModPeptideName};
      }
    } else {
      print STDERR "Error: $ARGV[0] format invalid \n";
    }
    exit unless $type;
    $head = 0;

    # Modifed to fetch indexes because of perl default array indexing could hide issues.
    # find_index gets up to 3 args: arrayref of possible headings, checked in order, common
    # name of item being sought (for debug messages), and (opt) whether it is OK to have null
    $idx{frg_z} = find_index( [qw(frg_z FragmentCharge ProductCharge)], 'Fragment ion charge', 1); 
    $idx{ion_s} = find_index( ['frg_type','FragmentType','FragmentIonType'], 'Ion series', 1 ); 
    $idx{ion_n} = find_index( ['frg_nr','FragmentNumber','FragmentSeriesNumber'], 'Ion series number', 1); 
    $idx{mseq} = find_index( ['modification_sequence','FullUniModPeptideName', 'PeptideModifiedSequence', 'ModifiedPeptide','IntModifiedPeptide' , 'peptideModSeq' ], 'Modified peptide sequence' ); 
    $idx{seq} = find_index( ['stripped_sequence','PeptideSequence', 'StrippedPeptide', 'peptide_sequence'], 'Stripped peptide sequence' ); 
    $idx{precursor} = find_index( ['Q1','PrecursorMz','q1'], 'Precursor m/z' ); 
    $idx{fragment} = find_index( ['Q3','ProductMz', 'FragmentMz','q3'], 'Fragment m/z' ); 
    $idx{rt} = find_index( ['RT_detected','Tr_recalibrated', 'iRT Value', 'iRT', 'rt_detected', 'irt'], 'Retention time' ); 
    $idx{pre_z} = find_index( [qw(prec_z PrecursorCharge)], 'Precursor ion charge' ); 
    # Missing value OK for these three
    $idx{protstr} = find_index( ['uniprot_id','protein_name','ProteinName','Protein Name','ProteinId', 'UniprotID', 'UniProtIds' ],'Protein name(s)',1 ); 
    $idx{decoy} = find_index( ['decoy'], 'Decoy',1 ); 
    $idx{lossType} = find_index( [qw( FragmentLossType ) ], 'Fragment neutral loss', 1 ); 
    $idx{transition_name} = find_index( ['transition_name'], 'transition_name',1 );
    $idx{IonMobility} = find_index( ['IonMobility'], 'Ion Mobility',1 );
    next;
  }

  # OS decoy skipping column triage
  if ( $opts{skip_decoys} ) {
    if ( $type eq 'os' ) {
      next if !defined $idx{decoy};
      if ( $line[$idx{decoy}] eq 'TRUE' || !$line[$idx{decoy}] ) {
        next;
      }
    }
  }

  ##  Extract assay-related data based on library type

  # Precursor charge
  my $pre_z = $line[$idx{pre_z}];
  $peps{pre_z}->{$pre_z}++;

  my $frg_z; my $ion_n; my $ion_s;
  # Fragment charge
  if($idx{frg_z}){
  	$frg_z = $line[$idx{frg_z}];
  }
  # Ion series
  if($idx{ion_s}){
	$ion_s = $line[$idx{ion_s}];
	if($ion_s =~ /-/){
		if ($ion_s =~ /(-17)/)  {
			$LTOS = '-17';
		} elsif ($ion_s =~ /(-18)/)  {
			$LTOS = '-18';
		} elsif ($ion_s =~ /(-64)/)  {
			$LTOS = '-64';
		}
			
	$ion_s =~ s/-17//g;
	$ion_s =~ s/-18//g;
	$ion_s =~ s/-64//g; #get rid of the -17, -18, -64  you get y10, y16, y10
	}else {
		$LTOS = '';
	}
  }

  # Ion series number
  if($idx{ion_n}){
	$ion_n = $line[$idx{ion_n}];
  }
  if (!defined $frg_z || !defined $ion_s || !defined $ion_n){
	if($idx{transition_name}){
		my $tr_name = $line[$idx{transition_name}];
		$tr_name =~ s/"DECOY_"//g; 
		my ( @tr_array ) = split( /_/, $tr_name );
		if (!defined $frg_z){
			# Fragment charge
			$frg_z = $tr_array[2];
		}
                # Ion series
                my $ion_sno = $tr_array[1];
                ##++++ 
                #OpenSwath library has y-1710, y-1816, y-6410 neutral losses for HN3, Water and oxidized methionine
                # so get rid of the -17, -18, -64 and consider these neutral losses in calculations.
                if($ion_sno =~ /-/){
			if ($ion_sno =~ /("-17")/)  {
				$LTOS = '-17';
			} elsif ($ion_sno =~ /("-18")/)  {
				$LTOS = '-18';
			} elsif ($ion_sno =~ /("-64")/)  {
				$LTOS = '-64';
			}
			$ion_sno =~ s/-17//g;
			$ion_sno =~ s/-18//g;
			$ion_sno =~ s/-64//g; #get rid of the -17, -18, -64  you get y10, y16, y10
		}else {
				$LTOS = '';
			}
		if (!defined $ion_s){
			if ($ion_sno =~ /([b y])/){			
				$ion_s = $1;
			}
		}
		
		if (!defined $ion_n){
			$ion_sno =~ s/[b y]//g;
			$ion_n = $ion_sno;
		}
	}
   }
      
   $peps{frg_z}->{$frg_z}++;
   $peps{frg_s}->{$ion_s}++;
   $peps{frg_n}->{$ion_n}++;		

  # Modified sequence
  my $mseq = $line[$idx{mseq}];
  $mseq =~ s/^_//g;
  $mseq =~ s/_$//g;

  # Stripped sequence
  my $seq = $line[$idx{seq}];

  # Retention time
  my $rt = $line[$idx{rt}];

  # Ion Mobility
  my $im = $line[$idx{IonMobility}];
  push (@imarray, $im);

  # Fragment m/z
  my $fragment = $line[$idx{fragment}];

  #Precursor m/z
  my $precursor = $line[$idx{precursor}];
  $precursorstats{$precursor}++;

  # Protein string (may be multiple).
  my $protstr = $line[$idx{protstr}];
  if ( !defined $idx{protstr} ) {
    $protstr = 'Default';
  }

  my $lossType = ( defined $idx{lossType} ) ? $line[$idx{lossType}] : '';
  if ( $lossType && !defined( $losses{$lossType} ) ) {
    print STDERR "unknown loss type $lossType specified, known values include:\n";
    for my $ls ( sort( keys( %losses ) ) ) {
      print STDERR join( "\t", $ls, $losses{$ls} ) . "\n";
    }
    exit 144;
  }

  if ( $opts{debug} ) {
    print STDERR "Dumping parsed values";
    print STDERR qq~
protstr: $protstr
precursor: $precursor
fragment: $fragment
rt: $rt
seq: $seq
mseq: $mseq
pre_z: $pre_z
frg_s: $ion_s
frg_n: $ion_n
frg_z: $frg_z

~;

    print STDERR "Dumping column mappings\n";
    for my $f ( sort( keys( %cmap ) ) ) {
      print STDERR join( "\t", $cmap{$f}, $f, $line[$cmap{$f}] ) . "\n"; 
    }
    exit;
  }


  # Mod seq + charge
  $pepkey = $mseq . '_' . $pre_z;
  $ion_cnt{$pepkey}++;

  if ( $type ne 'pv' && !$osw_altmod ) {
    $osw_altmod++ if $mseq =~ /\[\d+\]/;
  }

  # Do some clean-up on protein string (defined above)
  $protstr =~ s/\"//g;
  $protstr =~ s/^\d+\/*//;
  $protstr ||= 'DEFAULT';
  if ( $type eq 'os' ) {
    $protstr =~ s/;/,/g;
  }

  my @prot;
  if ( $protstr =~ /,/ ) {
    @prot = split( /,/, $protstr, -1 );
  } else {
    @prot = split( /\//, $protstr, -1 );
  }

# fragment, Q3 plus intensity 
#  my $fragkey = "$line[1] $line[5]";
#  my $fragkey = ( $type eq 'pv' ) ? "$line[1] $line[5]" : "$line[$cmap{ProductMz}] $line[$cmap{}]";
#  my $ion_s = ( $type eq 'pv' ) ? "$line[9]" : "$line[15]";

  $peps{t6_frg_s} ||= {};
  $peps{t6_frg_s}->{$ion_s}++ unless $ion_cnt{$pepkey} > 6;

  if ( $fragment > $precursor ) {
    $stats{fragment_above_precursor}++;
  }
  $stats{fragment}++;

  my $swakey = $pepkey;
  if ( $opts{swath_file} ) {
    if ( !$swaths{$swakey} ) {
      $swaths{$swakey} = [];
      for my $min ( keys( %swadefs ) ) {
        if ( $precursor >= $min && $precursor <= $swadefs{$min} ) {
          push @{$swaths{$swakey}}, { min => $min, max => $swadefs{$min} };
        }
      }
      if ( !scalar( @{$swaths{$swakey}} ) ) {
        $stats{swa_missing}++; 
        $problem_assays{$pepkey}++;
        next;
      } else {
        $stats{swa_defined}++ 
      }
    }
    my $is_bad = 0;
    for my $swa ( @{$swaths{$swakey}} ) {
      if ( $fragment >= $swa->{min} && $fragment <= $swa->{max} ) {
        print STDERR "$fragment ($precursor) in CONFLICT (suppressing further warnings)\n" if $opts{verbose} && !$stats{swa_conflict}; # Complain, but just once
        $stats{swa_conflict}++;
        $problem_assays{$pepkey}++;
        $is_bad++;
      } else {
        $stats{swa_ok}++;
      }
    }
    $stats{swa_conflict_assay} ||= {};
    $stats{swa_conflict_assay}->{$precursor}++ if $is_bad;
#    $stats{swa_conflict_precursor}++ if $is_bad;
    if ( abs( $precursor - $fragment ) <= 5 ) {
      $stats{swa_5}++;
    } elsif ( abs( $precursor - $fragment ) <= 25 ) {
      $stats{swa_25}++;
    }
  }

# Some hand-edited libraries might have blank lines...
  next if $line =~ /^\s*$/;

  if ( !defined $fragment || !defined $precursor ) {
    print STDERR "Unable to find prec a/o fragment ($precursor and $fragment) from $line\n";
    die;
  }

  $extrema{precursor_max} = $precursor if $precursor > $extrema{precursor_max};
  $extrema{fragment_max} = $fragment if $fragment > $extrema{fragment_max};
  $extrema{precursor_min} = $precursor if $precursor < $extrema{precursor_min};
  $extrema{fragment_min} = $fragment if $fragment < $extrema{fragment_min};


  if ( $curr && $curr ne $pepkey ) {

    my $fcnt = 0;
    for my $mz ( sort { $a <=> $b }( keys( %{$dta{frg}} ) ) ) {
      $peps{frg_i}->{cnt}++;
      $peps{frg_i}->{inten} += $dta{frg}->{$mz};
      last if $fcnt++ >= 5;
    }

    if ( $write_data ) {
      my $mw = sprintf( "%0.4f", $dta{pre}->{mz}/$dta{pre}->{z} );
      print MGF qq~BEGIN IONS
PEPMASS=$mw
CHARGE=$dta{pre}->{z}
TITLE=$pepkey
~;
      for my $mz ( sort { $a <=> $b }( keys( %{$dta{frg}} ) ) ) {
        print MGF join( "\t", $mz, $dta{frg}->{$mz} ) . "\n";
      }
      print MGF "END IONS\n";
    }
    if ( scalar(keys(%{$dta{frg}})) > 5 ) {
      $stats{frg_cnt_ok}++;
    } else {
      $stats{frg_cnt_short}++;
    }
    %dta = ( frg => {}, pre => {z => '', mz => ''} );
  }
  $curr = $pepkey;
  $dta{pre}->{mz} = $precursor;
  $dta{pre}->{z} = $pre_z;
  my $mz = sprintf( "%0.4f", $line[1]);
  my $int = $line[5];
  $dta{frg}->{$mz} = $int;

  $tmax{$pepkey} ||= { idx => 1,
                       max => 0,
                       max_idx => 1 };
  if ( $int > $tmax{$pepkey}->{max} ) {
    $tmax{$pepkey}->{max} = $int; 
    $tmax{$pepkey}->{max_idx} = $tmax{$pepkey}->{idx}; 
  }
  $tmax{$pepkey}->{idx}++;

  if ( $opts{assess_massdiff} ) {
    my $masskey = sprintf( "%0.4f", $precursor );
    $massmap{$pepkey} ||= {};
    my $skip = 0;
    my ( $pseq, $pchg ) = split( /_/, $pepkey );
    my $eseq = encode_modifications( $pseq );
#    print STDERR "$pseq becomes $eseq \n" if $pseq ne $eseq;
    my %ion2mz; 
    if ( $massmap{$pepkey}->{$masskey} ) {
      if ( !$massmap{$pepkey}->{$masskey}->{precursor_ok} ) {
        $problem_assays{$pepkey}++;
        
#        Misleading because known fragments labeled as na - revisit if performance becomes an issue
#        # Issue with this precursor already documented, no use checking fragment
#        $skip++ unless $opts{correct_mz};
      }
    } else {

      my $pmz = get_peptide_mass( $eseq, $pchg, 1 );
      my $ions = generate_ions( $pseq );
	  #print "IONHASH: $pseq\t";
      #print Dumper(\$ions);
      for ( my $fchg = 1; $fchg <= 7; $fchg++ ) {    ##7 is maximum fragment charge we consider usually it is upto 3
        for my $yion ( keys( %{$ions->{y}} ) ) {
          my $ion_key = 'y' . $yion . '_' . $fchg;
          my $intermediatey = encode_modifications($ions->{y}->{$yion});
          $ion2mz{$ion_key} = get_peptide_mass($intermediatey, $fchg );  
        }
        for my $bion ( keys( %{$ions->{b}} ) ) {
          my $ion_key = 'b' . $bion . '_' . $fchg;
          my $intermediateb = encode_modifications($ions->{b}->{$bion});
          $ion2mz{$ion_key} = get_peptide_mass( $intermediateb, $fchg ) - (WATER_MASS)/$fchg;
         #  print "\n1.Gen ionkey $ion_key\t $ion2mz{$ion_key}, $ions->{b}->{$bion} seq: $intermediateb,bion  fchg $fchg\n";
      }    
      }

      my $delta = abs( $precursor - $pmz );
      if ( $delta ) {
        my $dawn = $delta * 200000;
        if ( $dawn > $precursor - 0 ) { # More than 5ppm different
          $massmap{$pepkey}->{$masskey}->{precursor_ok} = 0;
          $problem_assays{$pepkey}++;
          if ( $opts{print_mass_error} ) {
             print ME join( "\t", "Precursor", $pepkey, $eseq, $pchg, $precursor, $pmz, $delta, $dawn ) . "\n";
          }
          print STDERR "$pepkey has precursor error with $delta from $pmz and $precursor\n" if $opts{verbose};
          $stats{precursor_bad}++;
        
#        Misleading because known fragments labeled as na - revisit if performance becomes an issue
#        # Issue with this precursor already documented, no use checking fragment
#        $skip++ unless $opts{correct_mz};
        } else {
          $massmap{$pepkey}->{$masskey}->{precursor_ok} = 1;
          $stats{precursor_ok}++;
        }
      } else {
        $massmap{$pepkey}->{$masskey}->{precursor_ok} = 1;
        $stats{precursor_ok}++;
      }
    }
    unless ( $skip ) {
      if ( !defined $massmap{$pepkey}->{$masskey}->{peaks} ) {
        $massmap{$pepkey}->{$masskey}->{peaks} = \%ion2mz;
      }
    }
    my $peak_key = $ion_s . $ion_n . '_' . $frg_z;
#    if ( $massmap{$pepkey}->{$masskey}->{peaks}->{$ion_s . $ion_n . '_' . $frg_z} ) {
    if ( $massmap{$pepkey}->{$masskey}->{peaks}->{$peak_key} ) {
      my $tmz = $massmap{$pepkey}->{$masskey}->{peaks}->{$peak_key};
      if ( $lossType ) {
        $tmz -= $losses{$lossType}/$frg_z;
      }
	   if ( $LTOS ) { 
        $tmz -= $losses{$LTOS}/$frg_z;
      }
	  
      my $tdelta = abs( $tmz - $fragment );
      $stats{delta_cnt}++;
      $stats{delta_sum} += $tdelta;
      if ( $tdelta ) {
        if ( $tdelta * 1000000 > $fragment ) { # More than 1ppm different
          $stats{fragment_bad}++;
          $problem_assays{$pepkey}++;
          if ( $opts{print_mass_error} ) {
             print ME join( "\t", "Fragment", $pepkey, $ion_s . $ion_n . '_' . $frg_z, $tmz, $fragment, $tdelta, $lossType ) . "\n";
          }
        } else {
          $stats{fragment_ok}++;
        }
      } else {
        $stats{fragment_ok}++;
      }
    } else {
      print STDERR "Can't find $peak_key for $pepkey and $masskey\n";
#      die Dumper( $massmap{$pepkey}->{$masskey}->{peaks} );

      $stats{fragment_na}++;
    }
  }


  # Cache trans per pep, pep per prot
  my $mult = ( scalar( @prot ) > 1 ) ? 'mult' : 'sing';

  my $fwd;
  my $decoy;
  for my $prot ( @prot ) {
    $prots{$prot} ||= { sing => {}, mult => {}, stripped => {} };
    $prots{$prot}->{$mult}->{$pepkey}++;
    $prots{$prot}->{stripped}->{$seq}++;
    if ( $prot =~ /DECOY/ ) {
      $decoy++;
    } elsif ( $opts{alt_decoy} && $prot =~ /$opts{alt_decoy}/ ) {
      $decoy++;
    } else {
      $fwd++;
    }
  }
  if ( $decoy ) { 
    if ( $fwd ) { 
      $decoy{mixed}++;
    } else {
      $decoy{decoy}++;
    }
  } else {
    $decoy{fwd}++;
  }

  $peps{$mult}->{$pepkey}++;
  if ( !$ion_s || ($ion_s !~ /y/i && $ion_s !~ /b/i ) ) {
    $peps{ion_s}->{$pepkey}++;
  }

  if ( !$ion_n || $ion_n < 3  ) {
    $peps{ion_n}->{$pepkey}++;
  }

# 0    PrecursorMz [ 778.4129855 ]
# 1    ProductMz [ 498.26707303 ]
# 2    Tr_recalibrated [ 57.3 ]
# 3    transition_name [ 6_AAAAAAAAAAAAAAAASAGGK_2 ]
# 4    CE [ -1 ]
# 5    LibraryIntensity [ 7324.6 ]
# 6    transition_group_id [ 1_AAAAAAAAAAAAAAAASAGGK_2 ]
# 7    decoy [ 0 ]
# 8    PeptideSequence [ AAAAAAAAAAAAAAAASAGGK ]
# 9    ProteinName [ 1/P0CG40 ]
# 10   Annotation [ b7/-0.009,b14^2/-0.009,m10:16/-0.009 ]
# 11   FullUniModPeptideName [ AAAAAAAAAAAAAAAASAGGK ]
# 12   PrecursorCharge [ 2 ]
# 13   GroupLabel [ light ]
# 14   UniprotID [ 1/P0CG40 ]
# 15   FragmentType [ b ]
# 16   FragmentCharge [ 1 ]
# 17   FragmentSeriesNumber [ 7 ]


  if ( $is_irt{$seq} ) {
    $rt{irt}->{$seq}++;
    $irt_cnt++;
  }
  $rt{mseqs}->{$mseq} ||= {};
  $rt{mseqs}->{$mseq}->{$pre_z} ||= $rt;
}
my @imsort= sort{ $a <=> $b } @imarray ; 
my $im_max= sprintf( '%0.4f', (pop@imsort));
my $im_min=sprintf( '%0.4f', (shift@imsort));


for my $ext ( keys( %extrema ) ) {
  $extrema{$ext} = sprintf( "%0.1f", $extrema{$ext} );
}

my $tmax_sum;
my $tmax_cnt;
for my $key ( keys( %tmax ) ) {
  $tmax_cnt++;
  $tmax_sum += $tmax{$key}->{max_idx};
}
my $tmax_avg = sprintf( "%0.1f", $tmax_sum/$tmax_cnt );

# Finish with last accumulated pepion
close ILIB;

if ( $opts{print_mass_error} ) {
  close ME;
}

# RT data calculation
my $n_irt = scalar( keys( %{$rt{irt}} ) );

my @rt_all;
my @rt_two;
my @rt_three;
my %pepions;
for my $mseq ( keys( %{$rt{mseqs}} ) ) {
  for my $ch ( keys( %{$rt{mseqs}->{$mseq}} ) ) {
    $pepions{$mseq . '_' . $ch}++;
    push @rt_all, $rt{mseqs}->{$mseq}->{$ch};
  }
  next unless $rt{mseqs}->{$mseq}->{2} && $rt{mseqs}->{$mseq}->{3};
  push @rt_two, $rt{mseqs}->{$mseq}->{2};
  push @rt_three, $rt{mseqs}->{$mseq}->{3};
}

@rt_all = sort { $a <=> $b } ( @rt_all );
my $rt_min = $rt_all[0];
my $rt_max = $rt_all[$#rt_all];
my $rt_med = $rt_all[int($#rt_all/2)];

my $rsq = 'na';
my $rt_five = 'na';
if ( scalar( @rt_two ) < 8 ) {
  print STDERR "Too few points: ". scalar( @rt_two ) . " to run regression \n";
} else {
  if ( $opts{rt_stats} ) {
    my $rt = $opts{output_file} . ".RT";
    if ( -e $rt ) { 
      print STDERR "RT stats file $rt exists, cannot overwrite\n";
      exit;
    }
    print STDERR "Printing +2/+3 pairs to $rt\n";
    open RT, ">$rt";
  }

  my $reg;
  if ( $has_regression ) {
    $reg = Statistics::Regression->new( 'Y', [ 'Z', 'X' ] );
  }
  my $imperfect = 0;
  my $ok = 0;
  my $ok_delta = 5;
  for ( my $i = 0; $i <= $#rt_two; $i++ ) {
    print RT join( "\t", $rt_two[$i], $rt_three[$i] ) . "\n" if $opts{rt_stats};

    if ( $has_regression ) {
      $reg->include( $rt_two[$i], [ 1.0, $rt_three[$i] ] );
    }
    $imperfect++ unless $rt_two[$i] == $rt_three[$i];
    $ok++ if abs( $rt_two[$i] - $rt_three[$i] ) <= $ok_delta;
  }
  close RT if $opts{rt_stats};
  $rsq = 1.000;
  if ( !$has_regression ) {
    $rsq = 'NoStats';
    print STDERR "Unable to run regression without Statistics::Regression\n" if $opts{verbose};
  } elsif ( $imperfect ) {
    print STDERR "Running regression with " . scalar( @rt_two ) . " Total points\n" if $opts{verbose};
    $rsq = sprintf( "%0.3f", $reg->rsq() );
  } else {
    print STDERR "Skipping regression with " . scalar( @rt_two ) . " identical values\n" if $opts{verbose};
  }
  if ( $rt_max < 600 ) {
    print STDERR "Assuming rt in minutes due to $rt_max\n";
  } else {
    print STDERR "Assuming rt in seconds\n";
    $ok_delta *= 60;
  }
  $rt_five = sprintf( "%0.1f", ($ok/scalar(@rt_two))*100 );
}

my $dtot;
for my $dtype ( qw( fwd decoy mixed ) ) {
  $dtot += $decoy{$dtype};
}
my $decoy_pct = ( $decoy{decoy} ) ? sprintf( "%0.1f", 100*($decoy{decoy}/$dtot)) : 0;
my $mixed_pct = ( $decoy{mixed} ) ? sprintf( "%0.1f", 100*($decoy{mixed}/$dtot)) : 0;
my $fwd_pct = ( $decoy{fwd} ) ? sprintf( "%0.1f", 100*($decoy{fwd}/$dtot)) : 0;


my $fcnt = 0;
for my $mz ( sort { $a <=> $b }( keys( %{$dta{frg}} ) ) ) {
  $peps{frg_i}->{cnt}++;
  $peps{frg_i}->{inten} += $dta{frg}->{$mz};
  last if $fcnt++ >= 5;
}

if ( $write_data ) {
  my $mw = sprintf( "%0.4f", $dta{pre}->{mz} );
  print MGF qq~BEGIN IONS
PEPMASS=$mw
CHARGE=$dta{pre}->{z}
TITLE=$pepkey
~;
  for my $mz ( sort { $a <=> $b }( keys( %{$dta{frg}} ) ) ) {
    print MGF join( "\t", $mz, $dta{frg}->{$mz} ) . "\n";
  }
  print MGF "END IONS\n";
  close MGF;
}


$peps{stripped_len} = {};
my %stripped;
for my $type ( qw( sing mult ) ) {
  for my $pep ( keys( %{$peps{$type}} ) ) {
    $stats{$type}->{pepcnt}++;
    my $tpep = $pep;
    $tpep =~ s/\[[^]]+\]//g;
	  $tpep =~ s/\([^)]+\)//g;
    $tpep =~ s/\d//g;
    $tpep =~ s/_//g;
    $tpep =~ s/-//g;
    if ( $tpep !~ /^\w+$/ ) {
      die "Issue with $pep: stripped value is $tpep\n";
    }
#    $stats{ 'str_' . $type };
    $stripped{$tpep}++;
    my $tlen = length( $tpep );
    $peps{stripped_len}->{$tlen}++;
    $stats{$type}->{len} += $tlen;
    $stats{$type}->{frags} += $peps{$type}->{$pep};
  }
}
my $mcstats = $opts{output_file} . ".mcstats.txt";
open MCSTATS, ">$mcstats";
print MCSTATS "trimmedpeptide\tlength\tMCcount\n";
# Which ppeps were seeen?  
my %seen = ( prots => {}, mprots => {}, peps => {}, mpeps => {} );
for my $sp ( keys( %stripped ) ) {
  if ( $ppeps{peps}->{$sp} ) {
    $seen{peps}->{$sp}++;
  } else {
    $seen{mpeps}->{$sp}++;
  }
  my $splen = length($sp);
  if ( $sp =~ /[KR](?!P)/ ) {
	$sp=~ s/^.//; #remove firstchar
	chop($sp); #remove last char As terminal K and R are not missed cleavages
	$sp =~ s/[RK](?!P)/z/g;
	my $spcount= ($sp =~ tr/z/z/);
	if ($spcount > 0){
 		$stats{mc}++;
		print MCSTATS "$sp\t$splen\t$spcount\n";
	  }
  }
}

#for my $p ( sort( keys( %prots ) ) ) {
#  my $pcnt = scalar( keys( %{$prots{$p}->{stripped}} ) );
#  print join( "\t", $p, $pcnt ) .   "\n";
#}
#exit;

# How many prots were seen
my %pepsperprot;
my @protcnts;
my $prottot;
my %protstats;
for my $prot ( keys( %prots ) ) {
  $pepsperprot{$prot} = scalar( keys( %{$prots{$prot}->{mult}} ) );
  $pepsperprot{$prot} += scalar( keys( %{$prots{$prot}->{sing}} ) );
  if ( $pepsperprot{$prot} < 0 ) {
    print STDERR "$prot\n";
#    print STDERR Dumper( $prots{$prot} );
#    exit;
  }
  $prottot += $pepsperprot{$prot};
  push @protcnts, $pepsperprot{$prot};
  $protstats{$pepsperprot{$prot}}++;
  if ( $ppeps{prots}->{$prot} ) {
    $seen{prots}->{$prot}++;
  } elsif ( $prot =~ /RT-Cal/ ) {
    $seen{prots}->{$prot}++;
  } else {
    $seen{mprots}->{$prot}++;
  }
}


my @sorted = sort( { $a <=> $b } @protcnts );
my $half = int( $#protcnts/2 + 1 );
my $med = $sorted[$half];
my $mean = sprintf( "%0.1f", $prottot/(1+$#protcnts) );
my $devs;
for ( my $i = 0; $i <= $#protcnts; $i++ ) {
  $devs += ( $protcnts[$i] - $mean )**2;
}
my $stddev = sprintf( "%0.1f", ($devs/(1+$#protcnts))**0.5 );

my $overprots = 0;
for my $cnt ( sort { $b <=> $a } keys( %protstats ) ) {
  last if $cnt < ( $mean + ( 3 * $stddev ) );
  $overprots += $protstats{$cnt};
}




$stats{totpeps} = scalar( keys( %stripped ) );
$stats{sprots} = scalar( keys( %{$seen{prots}} ) );
$stats{mprots} = scalar( keys( %{$seen{mprots}} ) );
if ( $stats{libprots} < 1 ) {
  $stats{sprots} = 'N/A';
  $stats{mprots} = 'N/A';
}

if ( $opts{print_ux} ) {
  my $ux = $opts{output_file} . ".UX_prots";
  if ( -e $ux ) { 
    print STDERR "Unexplained proteins file $ux exists, cannot overwrite\n";
    exit;
  }
  print STDERR "Printing 'unexplained' proteins to $ux\n";
  open UX, ">$ux";
  for my $p ( sort( keys( %{$seen{mprots}} ) ) ) {
    print UX "$p\n";
  }
  close UX;
}


#$stats{allprots} = $stats{sprots} + $stats{mprots};
$stats{allprots} = scalar( keys( %prots ) );

#print join( "\n", sort( keys( %{$seen{mprots}} ) ) ) . "\n";
$stats{speps} = scalar( keys( %{$seen{peps}} ) );
$stats{mpeps} = scalar( keys( %{$seen{mpeps}} ) );


#  print join( "\t", qw( Library Pepions Ptyp Mult Ptyp_perc Mult_perc alt_ser low_nr alen afrg afrglen Tot_peps MC_peps Lib_peps Seen_peps Lib_prots Seen_prots UX_prots ) ) . "\n";
# print join( "\t", qw( Library Pepions Ptyp Mult alt_ser low_nr alen afrg afrglen Ptyp_perc Mult_perc alt_ser_perc low_nr_perc ) ) . "\n";

  $stats{sing}->{pepcnt} ||= 0;
  $stats{mult}->{pepcnt} ||= 0;
  $stats{sing}->{len} ||= 0;
  $stats{mult}->{len} ||= 0;
  my $pepcnt = $stats{sing}->{pepcnt} + $stats{mult}->{pepcnt};
  for my $type ( qw( sing mult ) ) {
    $stats{$type}->{perc} = sprintf( "%0.2f", 100 * $stats{$type}->{pepcnt} / $pepcnt );
  }
  $stats{mult}->{frags} ||= 0;
  $stats{sing}->{frags} ||= 0;
  my $tot_len = $stats{sing}->{len} + $stats{mult}->{len};
  my $tot_frags = $stats{sing}->{frags} + $stats{mult}->{frags};

  my $alen = sprintf( "%0.1f", $tot_len/$pepcnt );
  my $afrg = sprintf( "%0.2f", $tot_frags/$pepcnt );
  my $afrglen = sprintf( "%0.2f", $afrg/$alen );


  # depricate s_cnt (alt series), instead do y and b seriest cnt/perc
#  my $s_cnt = scalar(keys(%{$peps{ion_s}} ) ); 
#  my $s_perc = ( $s_cnt ) ? '0.00' : sprintf( "%0.2f", 100*($s_cnt/$pepcnt));
  my $y_cnt = $peps{frg_s}->{y} || 0; 
  my $b_cnt = $peps{frg_s}->{b} || 0; 
  my $y_perc = ( $tot_frags ) ? sprintf( "%0.1f", ($y_cnt/$tot_frags)*100 ) : 0.0;
  my $b_perc = ( $tot_frags ) ? sprintf( "%0.1f", ($b_cnt/$tot_frags)*100 ) : 0.0;

  my $top_y = $peps{t6_frg_s}->{y} || 0;
  my $top_b = $peps{t6_frg_s}->{b} || 0;
  my $T6_y_perc = ( $top_y || $top_b ) ? sprintf( "%0.1f", (100*$top_y/($top_y+$top_b)) ) : 0.0;

  my $ainten = ( $peps{frg_i}->{cnt} ) ? sprintf( "%0.1f", $peps{frg_i}->{inten}/$peps{frg_i}->{cnt} ) : '0.0';

  my $tot_pre = 0;
  for my $z ( keys( %{$peps{pre_z}} ) ) {
    next if $z eq 'peps';
    $tot_pre += $peps{pre_z}->{$z};
  }
  my $pre_2 = ( $tot_pre && $peps{pre_z}->{2} ) ? sprintf( "%0.1f", 100*($peps{pre_z}->{2}/$tot_pre) ) : '0.0';
  my $pre_3 = ( $tot_pre && $peps{pre_z}->{3} ) ? sprintf( "%0.1f", 100*($peps{pre_z}->{3}/$tot_pre) ) : '0.0';

  # The number of missing or < 3 frg numbers
  my $n_cnt = scalar(keys(%{$peps{ion_n}})); 

  my $short_perc = ( $stats{frg_cnt_short} ) ? sprintf( "%0.1f", ($stats{frg_cnt_short}/($stats{frg_cnt_short}+$stats{frg_cnt_ok})*100 )) : 0; 

  my %mods;
  my %mcnt;
  my %mseen;


  for my $mpep ( keys( %{$peps{sing}}), keys(%{$peps{mult}} ) ) {

    # Only count one charge state
    $mpep =~ s/_\d+$//;
    next if $mseen{$mpep}++;
    $mcnt{all}++;
    my $regex = '(\w?\[[^\]]+\])'; 
    if ( $mpep =~ /UniMod/ ) {
      $regex = '(\w?\(UniMod:\d+\))';
    } else {
      next unless $mpep =~ /\[/;
    }
    my $has_mods = 0;
    for my $m ( $mpep =~ /$regex/ig ) {
      $mods{$m}++;
      $has_mods++;
    }
    $mcnt{mod}++ if ( $has_mods );
  }
  my $mcnt = 0;
  $peps{mods} = { none => $mcnt{all} - $mcnt{mod}, any => $mcnt{mod} };
  for my $mod ( sort( keys( %mods ) ) ) {
    $mcnt += $mods{$mod};
    $peps{mods}->{$mod} = $mods{$mod};
  }
  my $mpct = sprintf( "%0.1f", (($mcnt{mod}/$mcnt{all})*100) );
  print STDERR "Saw $mcnt{mod} peptides with modifications, $mcnt total mods in those peptides, and $mcnt{all} distinct modified peptides (includes non-mod)\n";

  my $tot_fragment = $stats{fragment};
  my $fragment_above_pct = ( $stats{fragment_above_precursor} ) ? sprintf( "%0.3f", $stats{fragment_above_precursor}/$stats{fragment}) : 0;

  my @checks;
  if ( $opts{assess_massdiff} ) {
    @checks = ( qw( precursor_ok precursor_bad fragment_ok fragment_bad fragment_na fragment_avg_mdiff ) );
  }


  my $problem_assays = scalar( keys( %problem_assays ) );
  my $n_perc = ( $n_cnt ) ? '0.00' : sprintf( "%0.2f", 100*($n_cnt/$pepcnt));
  my $pepion_cnt = scalar(keys(%pepions)) || 0;
  my @headings = ( qw( library
                        format
                        pepions
                        fragments
                        ptp_percent
                        shared_percent
                        shared_pepions
                        peptides
                        mod_peps
                        mod_percent
                        total_mods
                        chg_2
                        chg_3
                        precursor_min
                        precursor_max
                        fragment_min
                        fragment_max
                        avg_len
                        avg_num_frags
                        avg_frag_len
                        short_perc
                        fragment_above_precursor
                        y_perc
                        b_perc
                        t6_y_perc
                        avg_intensity
                        rt_min
                        rt_max
                        rt_med
                        rt_rsq
                        rt_five
                        n_irt
                        irt_cnt
                        low_nr
                        db_peps
                        seen_peps
                        mc_peps
                        db_prots
                        lib_prots
                        seen_prots
                        ux_prots
                        decoy_pct
                        mixed_pct
                        fwd_pct
                        med_ppp
                        mean_ppp
                        stddev_ppp
                        3_sigma_ppp
                        max_intensity_idx
                        ), @checks, @swaths, 'problem_assays' , 'im_min', 'im_max' );


  my @checkvals = ();
  if ( $opts{assess_massdiff} ) {
    $stats{fragment_avg_mdiff} = ( $stats{delta_cnt} ) ? sprintf( "%0.4f", $stats{delta_sum}/$stats{delta_cnt}) : 0;
    for my $k ( @checks ) {
      push @checkvals, $stats{$k} || 0;
    }
  }

  my @swathvals = ();
  if ( $opts{swath_file} ) {
    for my $skey ( @swaths ) {
      if ( $skey eq 'swa_conflict_assay' ) {
        push @swathvals, scalar(keys(%{$stats{$skey}} ) );
      } else {
        push @swathvals, $stats{$skey};
      }
    }
  }

  $rt_min = sprintf( "%0.1f", $rt_min );
  $rt_max = sprintf( "%0.1f", $rt_max );
  $rt_med = sprintf( "%0.1f", $rt_med );

  my %type2sw = ( os => 'OpenSWATH', 
                  pv => 'Peakview', 
                  sn => 'Spectronaut',
                  pf => 'Prosit/Fusion' );

  my @out_values = ( 
              $libname,
              $type2sw{$type},
              $pepion_cnt,
              $tot_fragment,
              $stats{sing}->{perc},
              $stats{mult}->{perc},
              $stats{mult}->{pepcnt},
              $stats{totpeps},
              $mcnt{all},
              $mpct,
              $mcnt,
              $pre_2,
              $pre_3,
              $extrema{precursor_min},
              $extrema{precursor_max},
              $extrema{fragment_min},
              $extrema{fragment_max},
              $alen,
              $afrg,
              $afrglen,
              $short_perc,
              $fragment_above_pct,
              $y_perc,
              $b_perc,
              $T6_y_perc,
              $ainten,
              $rt_min,
              $rt_max,
              $rt_med,
              $rsq,
              $rt_five,
              $n_irt,
              $irt_cnt,
              $n_cnt,
              $stats{libpeps},
              $stats{speps},
              $stats{mc},
              $stats{libprots},
              $stats{allprots},
              $stats{sprots},
              $stats{mprots},
              $decoy_pct,
              $mixed_pct,
              $fwd_pct,
              $med,
              $mean,
              $stddev,
              $overprots,
# die "over is $overprots med is $med, mean is $mean, and stddev is $stddev\n";
              $tmax_avg,
              @checkvals,
              @swathvals,
              $problem_assays,
	      $im_min,
              $im_max
             );

  my $defs = get_coldefs();
  my @defs;
  for my $head ( @headings ) {
    push @defs, $defs->{$head} || '';
  }

#  if ( !$opts{redirect} && !($opts{filter_assays} || $opts{correct_mz} ) ) {
  if ( !$opts{redirect} ) {
    open OUT, ">$opts{output_file}.QC.tsv";
    print STDERR "Opening $opts{output_file}.QC.tsv for writing!\n";
    if ( $opts{invert_output} ) {
      for ( my $i = 0; $i <= $#headings; $i++ ) {
        print OUT join( "\t", $headings[$i], $out_values[$i] );
        if ( $opts{coldefs} ) {
          print OUT "\t$defs[$i]";
        }
        print OUT "\n";
      }
    } else {
      print OUT join( "\t", @headings ) . "\n";                     
      print OUT join( "\t", @out_values ) . "\n";                     
      if ( $opts{coldefs} ) {
        print OUT join( "\t", @defs ) . "\n";
      }
    }
    close OUT;
#  } elsif ( !($opts{filter_assays} || $opts{correct_mz} ) ) {
  } else {
    if ( $opts{invert_output} ) {
      for ( my $i = 0; $i <= $#headings; $i++ ) {
        print join( "\t", $headings[$i], $out_values[$i] );
        if ( $opts{coldefs} ) {
          print "\t$defs[$i]";
        }
        print "\n";
      }
    } else {
      print join( "\t", @headings ) . "\n";                     
      print join( "\t", @out_values ) . "\n";                     
      if ( $opts{coldefs} ) {
        print join( "\t", @defs ) . "\n";
      }
    }
  }

  if ( !$problem_assays ) {
    print STDERR "No problem assays found!\n";
  } else { 
    print STDERR "$problem_assays problem assays found!\n";
  }

  # User has requested a clean library
  if ( $opts{filter_assays} && $problem_assays ) {

    my $clean = $opts{ion_library} . ".clean";
    my $problem = $opts{ion_library} . ".problem";
    if ( -e $problem ) { 
      print STDERR "problem library $problem exists, cannot overwrite\n";
    } elsif ( -e $clean ) { 
      print STDERR "clean library $clean exists, cannot overwrite\n";
    } else {
      open ILIB, $opts{ion_library}|| die;

      print STDERR "Printing problem assays to $problem\n";
      print STDERR "Printing good assays to $clean\n";
      open CLEANLIB, ">$clean" || die;
      open PROBLEMLIB, ">$problem" || die;
      my $head = 1;
      while ( my $line = <ILIB> ) {
        if ( $head ) {
          print CLEANLIB $line;
          print PROBLEMLIB $line;
          $head = 0;
          next;
        }
        chomp $line;
#        my @line = split( /\t/, $line, -1 );
        my @line = parse_line( $line );

        my $pre_z = $line[$idx{pre_z}];
        my $mseq = $line[$idx{mseq}];

        my $pepkey = $mseq . '_' . $pre_z;

        if ( $problem_assays{$pepkey} ) {
          print PROBLEMLIB "$line\n";
        } else {
          print CLEANLIB "$line\n";
        }
      }
      close ILIB;
      close CLEANLIB;
      close PROBLEMLIB;
    }
  }

  # User has requested an mz-corrected library
  if ( $opts{correct_mz} ) {

    if ( $stats{fragment_na} ) {
      print STDERR "Some ($stats{fragment_na}) unknown fragment types were found, unable to do mz correction\n";
      exit;
    }
    unless ( $stats{fragment_bad} || $stats{precursor_bad} ) {
      print STDERR "No m/z discrepancy with either precursor or fragment ions was found, cntl-d to stop run (10 seconds)"; 
      sleep 10;
    }

    my $mz_corr = $opts{ion_library} . ".mz_corrected";
    if ( -e $mz_corr ) { 
      print STDERR "m/z corrected library $mz_corr exists, cannot overwrite\n";
    } else {
      open ILIB, $opts{ion_library}|| die;

      print STDERR "Printing mz-adjusted assays to $mz_corr\n";
      open MZLIB, ">$mz_corr" || die;
      my $head = 1;
      my %mz_stats;
      while ( my $line = <ILIB> ) {
        if ( $head ) {
          print MZLIB $line;
          $head = 0;
          next;
        }
        chomp $line;
#        my @line = split( /\t/, $line, -1 );
        my @line = parse_line( $line );

        # Pull out needed info
        my $precursor = sprintf( "%0.4f", $line[$idx{precursor}] );
        my $fragment = sprintf( "%0.4f", $line[$idx{fragment}] );
        my $pre_z = $line[$idx{pre_z}];
        my $frg_z = $line[$idx{frg_z}];
        my $ion_s = $line[$idx{ion_s}];
		# Ion series
		##++++ 
		#OpenSwath library has y-17, y-18, y-64 neutral losses for HN3, Water and oxidized methionine
		# so get rid of the -17, -18, -64 and consider these neutral losses in calculations.

		my $LTOS;
		if ($ion_s =~ /(-17)/)  {
			$LTOS = '-17';
		} elsif ($ion_s =~ /(-18)/)  {
			$LTOS = '-18';
		} elsif ($ion_s =~ /(-64)/)  {
			$LTOS = '-64';
		}else {
			$LTOS = '';
		}
		$ion_s =~ s/[-17 -18 -64]//g; #get rid of the -17, -18, -64 

        my $ion_n = $line[$idx{ion_n}];
        my $mseq = $line[$idx{mseq}];
        $mseq =~ s/^_//g;
        $mseq =~ s/_$//g;
        my $ion_key = $ion_s . $ion_n . '_' . $frg_z;
        my $lossType = ( defined $idx{lossType} ) ? $line[$idx{lossType}] : '';

        my $pepkey = $mseq . '_' . $pre_z;
        my $eseq = encode_modifications( $mseq );
        my $pmz = get_peptide_mass( $eseq, $pre_z );
        my $imz = $massmap{$pepkey}->{$precursor}->{peaks}->{$ion_key};
		
        # Adjust for neutral losses (if any)
        $imz -= $losses{$lossType}/$frg_z if ( $lossType );
		$imz -= $losses{$LTOS}/$frg_z if ( $LTOS );
#        print STDERR qq~ PRE: $precursor PRZ: $pre_z FRG: $fragment FRZ: $frg_z MSQ: $mseq PKY: $pepkey INS: $ion_s INN: $ion_n INK: $ion_key pMZ: $pmz fMZ: $imz ~;          

        if ( !$pmz ) {
          die "Unable to calculate precursor m/z from $line\n";
        } elsif ( !$imz ) {
          die "Unable to calculate fragment ion m/z from $line\n";
        }
        if ( $pmz != $precursor ) {
          $mz_stats{prec_changed}++;
        } else {
          $mz_stats{prec_ok}++;
        }
        if ( $imz != $fragment ) {
          $mz_stats{frag_changed}++;
        } else {
          $mz_stats{frag_ok}++;
        }

        # Calc and assign fragment m/z
        $line[$idx{precursor}] = $pmz;
        $line[$idx{fragment}] = $imz;

        if ( $is_CSV ) {
          print MZLIB $quote_char . join( $quote_char . ',' . $quote_char, @line ) . $quote_char . "\n";
        } else {
          print MZLIB join( "\t", @line ) . "\n";
        }
      }
      close ILIB;
      close MZLIB;
      for my $s ( sort( keys( %mz_stats ) ) ) {
        print STDERR "$s\t$mz_stats{$s}\n";
      }
    }
  }

  my $t1 = time();
  my $runtime = $t1 - $t0;
  print STDERR "Finished run in $runtime seconds\n";
  exit unless $opts{full_stats};

  my $fullstats = $opts{output_file} . ".fullstats";
  if ( -e $fullstats ) { 
    print STDERR "full_stats file $fullstats exists, cannot overwrite\n";
    exit;
  }
  print STDERR "Printing full stats to $fullstats\n";

  open FULL, ">$fullstats";
  for my $key ( qw( frg_s mods ) ) {
    my $line = "$key:";    
    my $sep = '';
    for my $val ( sort ( keys( %{$peps{$key}} ) ) ) {
      my $item = "$sep$val=$peps{$key}->{$val}";
      $item =~ s/:/_/g;
      $line .= $item;
      $sep = ',';
    }
    print FULL "$line\n";
  }

  

  for my $key ( qw( pre_z frg_z frg_n stripped_len ) ) {
    my $line = "$key:";    
    my $sep = '';
    my @keys = sort { $a <=> $b } keys( %{$peps{$key}} );
    my $min = $keys[0];
    my $max = $keys[$#keys];
    for ( my $idx = $min; $idx <= $max; $idx++ ) {
      if ( $opts{zero_pad} ) {
        $peps{$key}->{$idx} ||= 0;
      } else {
        next unless defined $peps{$key}->{$idx};
      }
      $line .= "$sep$idx=$peps{$key}->{$idx}";
      $sep = ',';
    }
    print FULL "$line\n";
  } 
  my $pepsing = $peps{sing}; 
  my $pepmult = $peps{mult};
  my %pepchg;
  my @pep = (keys %{$pepsing}, keys %{$pepmult}); #works with all perl versions 
  foreach $_(@pep){ #works with all perl versions
    my ( $mseq, $chg ) = split( /_/, $_ ); #works with all perl versions
    $pepchg{$chg}++;
  }
  my $sep = '';
  my @keys = sort { $a <=> $b } keys( %pepchg );
  my $min = $keys[0];
  my $max = $keys[$#keys];
  my $line = 'Distinct_pre_z:';
  for ( my $idx = $min; $idx <= $max; $idx++ ) {
    if ( $opts{zero_pad} ) {
      $pepchg{$idx} ||= 0;
    } else {
      next unless defined $pepchg{$idx};
    }
    $line .= "$sep$idx=$pepchg{$idx}";
    $sep = ',';
  }
  print FULL "$line\n";

# This section causes the numbers > 5 to be summed into the 5 bin.
#  $protstats{$pepsperprot{$prot}}++;
  my %top5;
  for my $p ( keys( %protstats ) ) {
    if ( $p >= 5 ) {
      $top5{5} += $protstats{$p};
    } else {
      $top5{$p} = $protstats{$p};
    }
  }
#  %protstats = %top5;

 
  my $sep = '';
  my @keys = sort { $a <=> $b } keys( %protstats );
  my $min = $keys[0];
  my $max = $keys[$#keys];
  my $line = 'PepIonsPerProtein:';
  for ( my $idx = $min; $idx <= $max; $idx++ ) {
    if ( $opts{zero_pad} ) {
      $protstats{$idx} ||= 0;
    } else {
      next unless defined $protstats{$idx};
    }
    $line .= "$sep$idx=$protstats{$idx}";
    $sep = ',';
  }
  print FULL "$line\n";

  my %real_ppp;
  for my $prot ( keys( %prots ) ) {
    my $scnt = scalar( keys( %{$prots{$prot}->{stripped}} ) );
    # This statement causes the numbers > 5 to be summed into the 5 bin.
    # $scnt = 5 if $scnt > 5;
    $real_ppp{$scnt}++;
  }
  my $sep = '';
  my @keys = sort { $a <=> $b } keys( %real_ppp );
  my $min = $keys[0];
  my $max = $keys[$#keys];
  my $line = 'PepsPerProtein:';
  for ( my $idx = $min; $idx <= $max; $idx++ ) {
    if ( $opts{zero_pad} ) {
      $real_ppp{$idx} ||= 0;
    } else {
      next unless defined $real_ppp{$idx};
    }
    $line .= "$sep$idx=$real_ppp{$idx}";
    $sep = ',';
  }
  print FULL "$line\n";
      
  
#  $precursorstats{$pepsperprot{$prot}}++;
  my $sep = '';
  my %pre_bin;
  for my $pr ( keys( %precursorstats ) ) {
    for ( my $bin = 50; $bin <= 4000; $bin += 50 ) {
      if ( $bin > $pr ) {
        $pre_bin{$bin}++;
        last;
      }
    }
  }
  my $line = 'PrecusorCounts:';
  for my $pre ( sort { $a <=> $b } keys(%pre_bin)  ) {
    $line .= "$sep$pre=$pre_bin{$pre}";
    $sep = ',';
  }
  print FULL "$line\n";
  
#  $ion_cnt{pepkey}++; 
  my $sep = '';
  my %cnt_bins;
  for my $mpep ( keys( %ion_cnt ) ) {
    $cnt_bins{$ion_cnt{$mpep}}++;
  }
  my @keys = sort { $a <=> $b } keys( %cnt_bins );
  my $min = $keys[0];
  my $max = $keys[$#keys];
  my $line = 'FragmentsPerPrecursor:';
  for ( my $idx = $min; $idx <= $max; $idx++ ) {
    if ( $opts{zero_pad} ) {
      $cnt_bins{$idx} ||= 0;
    } else {
      next unless defined $cnt_bins{$idx};
    }
    $line .= "$sep$idx=$cnt_bins{$idx}";
    $sep = ',';
  }
  print FULL "$line\n";

  close FULL;
# End MAIN program #

#+ 
# Simple CSV parser
#-
sub parse_csv {
  return quotewords( ",", 0, $_[0] );
}

#+ 
# Line parser
#-
sub parse_line {
  my $line = shift;
  if ( $is_CSV ) {
    return parse_csv( $line );
  } else {
    return split( "\t", $line, -1 );
  }
}

#+
# Routine to get index of specific cmap fields 
#-
sub find_index {
  my $names = shift;
  my $target = shift;
  my $fail_ok = shift || 0;
  my $idx;
  for my $name ( @{$names} ) {
    if ( defined( $cmap{$name} ) ) {
      return $cmap{$name};
    }
  }
  return if $fail_ok;
  my $opts = join( ',', @{$names} );
  print STDERR "Unable to find an index for $target in $opts; printing debug info\n\n";
  $opts{debug}++;
}

#+
# routine to read and validate options
#-
sub read_options {

  # Local hash of options
  my %opts;

  GetOptions(\%opts,"ion_library:s",
                    "alt_decoy:s",
                    "peptide_file:s",
                    "write",
                    "output_file:s",
                    "filter_assays",
                    "correct_mz",
                    "full_stats",
                    "rt_stats",
                    "verbose", 
                    'assess_massdiff', 
                    'print_mass_error', 
                    'skip_decoys', 
                    'swath_file:s', 
                    'coldefs',
                    'debug',
                    'print_ux', 
                    'invert_output', 
                    'neutral_loss', 
                    'zero_pad',
                    'help' ) || die "$'";

  $opts{mono} = { G => 57.021464,
                  D => 115.02694,
                  A => 71.037114,
                  Q => 128.05858,
                  S => 87.032029,
                  K => 128.09496,
                  P => 97.052764,
                  E => 129.04259,
                  V => 99.068414,
                  M => 131.04048,
                  T => 101.04768,
                  H => 137.05891,
                  C => 103.00919,
				  U => 150.95363,#U = selenocysteine
                  F => 147.06841,
                  L => 113.08406,
                  R => 156.10111,
                  I => 113.08406,
                  N => 114.04293,
                  Y => 163.06333,
                  W => 186.07931,
				  X => 0
		};
   
  my $err = '';

  # Check for required params
  for my $opt ( qw( ion_library  ) ) {
    $err .= "Missing required option $opt\n" unless $opts{$opt};
  }

  if ( $opts{help} ) {
    print_usage( "" );
  } elsif ( $err ) {
    print_usage( $err );
  }

  # Make sure ion library file exists
  if ( ! -e $opts{ion_library} ) {
    print_usage( "File $opts{ion_library} does not exist\n" );
  }

  # Determine filename base
  if ( !$opts{output_file} ) {
    if ( -t STDOUT ) {
      my $qc = "$opts{ion_library}.QC";
      if ( -e $qc ) {
        print STDERR "QC file $qc exists, cannot overwrite - printing to STDOUT\n";
      } else {
         $opts{output_file} = $qc;
      }
    } else {
      $opts{redirect}++;
    }
    $opts{output_file} = $opts{ion_library};
  }


  if ( $opts{filter_assays} ) {
    for my $o ( qw( correct_mz full_stats rt_stats print_ux ) ) {
      die "filter_assays is incompatible with correct_mz, full_stats, rt_stats, and print_ux" if $opts{$o};
    }
    if ( !$opts{assess_massdiff} ) {
      print STDERR "filter_assays option requires assess_massdiff, enabling\n";
      $opts{assess_massdiff}++;
    }
  } elsif ( $opts{correct_mz} ) {
    for my $o ( qw( full_stats rt_stats print_ux ) ) {
      die "correct_mz is incompatible with full_stats, rt_stats, and print_ux" if $opts{$o};
    }
    if ( !$opts{assess_massdiff} ) {
      print STDERR "correct_mz option requires assess_massdiff, enabling\n";
      $opts{assess_massdiff}++;
    }
  }
   
  return ( %opts );
}


#+
#  Routine to print usage statement if requested or error found
#-
sub print_usage {
  my $msg = shift || "";
  my $prog = basename ( $0 );

  print <<"  EOM";

  $msg
  Usage: $prog --ion_library libfile [ --peptide_file pepfile  --output output_file ]

  # Required
  ion_library    Path to ion library file in Peakview/OpenSWATH/Spectronaut/Prosit format (required)

  # Commonly used
  assess_massdiff   Compare library masses to theoretical monoisotopic values
  c, coldefs        print out column definitions
  f, full_stats     Print file (output_file.fullstats)
  h, help           Print this usage information and exit 
  invert_output     Invert output, 3 columns by X rows
  o, output_file    File (base) to which to write data; given MyOut, generates 
                    MyOut.QC.txt, MyOut.fullstats, etc. If not specified, library name 
                    is used.
  peptide_file      Path to peptide digest file, format is protein(s)\tpeptide
  swath_file        Path to file of SWATH definitions, format min_m/z\tmax_m/z

  # Generate filtered/corrected library
  correct_mz        Will compute precursor and fragment masses for each entry
                    in library and print new library:
                    ion_library.mz_corrected
                    Requires use of --assess, precludes use of 
                    filter_assays, rt_stats, full_stats, print_ux 
  filter_assays     Assesses each assay, and removes any that have a 
                    problem (m/z error, discrepancy from SWATH file).
                    Creates 2 new libraries, one with all passing assays and 
                    another with all failed (problem) assays:
                    ion_library.clean 
                    ion_library.problem
                    Requires use of --assess, precludes use of correct_mz, 
                    rt_stats, full_stats, print_ux

  # Other
  alt_decoy         Alternate decoy prefix (default is DECOY)
  debug             Print parsed/library values in case of format issues
  print_mass_error  Print out ions whose m/z values differ significantly from 
                    theoretical
  print_ux          Print unexplained proteins (seen in ion library but not 
                    in peptide_file if provided
  rt_stats          Print file of +2/+3 RTs for analysis to output_file.RT,
                    only if enough for regression
  skip_decoys       Skip DECOY entries when computing statistics
  v, verbose        Verbose output mode
  write             Write  MGF file with info from each spectrum (rare)
  zero_pad          Print intermediate values=0 with full_stats.
  EOM

  exit;

}



#+
#  Routine to generate output column definitions
#-
sub get_coldefs {
  my %defs = (
library => 'Name of library file being analyzed',
format => 'Library format, one of OpenSWATH, Peakview, or Spectronaut ',
pepions => 'Number of peptide ions (i.e. precursor, sequence + modifications + charge)',
fragments => 'Number of fragment (fragment ions) in library',
ptp_percent => 'Percentage of proteotypic pepions (not shared)',
shared_percent => 'Percentage of shared peptide ions (pepions)',
shared_pepions => 'Number of shared pepions',
peptides => 'Number of distinct peptide sequences',
mod_peps => 'Number of distinct modified peptides (sequences + modifications)',
mod_percent => 'percentage of distinct modified peptides with a mass modification',
total_mods => 'Number of mass modified amino acids',
chg_2 => 'Percentage of charge 2 precursors',
chg_3 => 'Percentage of charge 3 precursors',
avg_len => 'Average peptide Length',
precursor_min => 'Minimum precursor m/z (mass/charge) in library',
precursor_max => 'Maximum precursor m/z in library',
fragment_min => 'Minimum fragment m/z in library',
fragment_max => 'Maximum fragment m/z in library',
avg_num_frags => 'Average number of fragment per assay (precursor)',
avg_frag_len => 'Average fragment sequence length',
short_perc => 'Percentage of assays with 5 or fewer transitions',
y_perc => 'Percentage of y ions',
b_perc => 'Percentage of b ions',
t6_y_perc => 'Percentage of y  ions considering only the top 6 fragments per assay',
avg_intensity => 'Average Intensity',
low_nr => 'Number of fragment annotated as y1,y2,b1, or b2',
fragment_above_precursor => 'Percentage of fragment m/z above precursor m/z',
max_intensity_idx => 'Average index of most intense fragment ',
med_ppp => 'Median number of pepions per protein ',
mean_ppp => 'Mean number of pepions per protein ',
stddev_ppp => 'Standard deviation of the number of pepions per protein ',
'3_sigma_ppp' => 'number of pepions per protein more than 3 standard deviations from the mean ',
rt_min => 'Minimum RT (retention time) in library ',
rt_max => 'Maximum RT in library',
rt_med => 'Median RT in library',
rt_rsq => 'r-squared value of fit between RT of +2 and +3 charge states for the same modified peptide',
rt_five => 'Percentage of +2/+3 charge pairs of the same mod pep within 5 RT units of each other',
n_irt => 'Number of iRT peptides in library',
irt_cnt => 'Number of iRT assays (precursor + fragment)',
db_peps => 'Peptides in reference library, often 7-50 AA',
seen_peps => 'Number of library peptides seen',
mc_peps => 'Number of Missed cleavage peptides',
db_prots => 'Number of Proteins in reference library',
lib_prots => 'Number of proteins with at least one assay in ion library',
seen_prots => 'Number of library proteins seen',
ux_prots => 'Unexplained (not in reference library) proteins',
decoy_pct => 'Percentage of decoy (optionally includes "alternative", non-db decoys) assays',
mixed_pct => 'Percentage of mixed decoy/target (has both decoy and no-decoy annotations) assays',
fwd_pct => 'Percentage of target (non-decoy) assays',
precursor_ok => 'Number of assays where precursor is within 5 PPM (parts per million m/z) of theoretical',
precursor_bad => 'Number of assays where fragment is more than 5 PPM from theoretical',
fragment_ok => 'Number of assays where fragment is within 1 PPM of theoretical',
fragment_bad => 'Number of assays where precursor is more than 1 PPM from theoretical',
fragment_na => 'Number of assays where peak annotation not found in expected b/y series',
fragment_avg_mdiff => 'Average m/z difference between reported and theoretical fragment',
swa_defined => 'Number of peptide ions that fall into a defined SWATH bin ',
swa_missing => 'Number of peptide ions that fail to fall into a defined SWATH bin ( Pepions - swa_defined)',
swa_conflict => 'Number of fragment_ions that fall into same SWATH(s) as precursor',
swa_ok => 'Number of fragment_ions that do not fall into same SWATH(s) as precursor',
swa_conflict_assay => 'Number of precursor that have at least one failing fragment',
swa_5 => 'Number of fragment ions that fall within 5 Th of precursor ion',
swa_25 => 'Number of fragment ions that fall within 25 Th of precursor ion',
problem_assays => 'Assays whose precursor or any fragment m/z values do not match SWATHs file or differ significantly from theoretical values',
im_min => 'Minimum ion mobility reported in library',
im_max => 'Maximum ion mobility reported in library'
      );
  return \%defs;
}


#+
#  Routine to calculate peptide mz
#-
sub get_peptide_mass {
	my $eseq = shift; # single-letter encoded sequence!

  	my $chg = shift || 2;
  	my $is_precursor = shift || 0;
	# print STDERR "Eseq is $eseq, charge is $chg\n" if $is_precursor;

	my $mz_key = $eseq . '_' . $chg;

	#return $mpep2mz{$mz_key} if $mpep2mz{$mz_key}; 
	my $mass = 0;
          
        if ($eseq =~ /[a-z]/) {
            	 for my $keys3 (keys %mods_current)
			{
	 		 $mass += $known_mods{$keys3}{'mz'} * $mods_current{$keys3};
  			}
	} else {
	 	  %mods_current = ();
 	    }

	$eseq =~ s/x//g;
	for my $aa ( split( '', uc($eseq) ) ){
		die "Unknown AA $aa from $eseq" unless $opts{mono}->{$aa};
		$mass += $opts{mono}->{$aa};
		}
	$mass += WATER_MASS;
	my $h_mass = 1.00727646688;
	my $pmz = sprintf( '%0.4f', ( $mass + $chg * H_MASS )/$chg );
	#  print STDERR "eseq is $eseq, ycnt is $ycnt, mass is $mass, pmz is $pmz\n" if $is_precursor;
	return $pmz;
}

#+
#  Routine to translate between modification namespaces
#-
sub encode_modifications {
  my $seq = shift;
  my $regex = '(\w?\[[^\]]+\]-*)'; 

  # Block for debugging mod regex
  if ( 0 ) {
    for my $mod ( keys( %known_mods ) ) {
      my $match = '0';
      if ( $mod =~ /$regex/ ) {
        $match = $1;
      }
      print "$mod\t$match\n";
    }
    exit;
  }


  if ( $seq =~ /UniMod/ ) {
    $regex = '(\w?\(UniMod:\d+\))';
  } 

  %mods_current=();
  for my $m ( $seq =~ /$regex/ig ) {
    $mods_current{$m}++;

    unless( $known_mods{$m}{'aa'} ) {
      print ">$m<";
      die "Unknown modification $m in $seq ($regex)\n";
    }
  }
  my $eseq = $seq;
  for my $m ( keys( %mods_current ) ) {
    my $em = $m;
    $em =~ s/\[/\\[/g;
    $em =~ s/\+/\\+/g;
    $em =~ s/\]/\\]/g;
    $em =~ s/\(/\\(/g;
    $em =~ s/\)/\\)/g;
    $eseq =~ s/$em/$known_mods{$m}{'aa'}/g;
  }
             
  if ( $eseq =~ /(.?\[[^\[\]]+\])/ ) {
    die "1. Unknown mod $1 in $eseq\n";
  } elsif ( $eseq =~ /\d/ ) {
    die "2. Unknown mod in $eseq\n";
  } elsif ( $eseq =~ /\[/ ) {
    die "3. Unknown mod in $eseq\n";
  } elsif ( $eseq =~ /\(/ ) {
    die "4. Unknown mod in $eseq\n";
  }
  return $eseq;
}

#+
#  Routine to generate fragment ions
#-
sub generate_ions {
  my $eseq = shift;
  
  my %b;
  my %y;
  my @aa;
  
  my $regex = '(\w?\[[^\]]+\]-*)'; 
  if ( $eseq =~ /UniMod/ ) {
    $regex = '(\w?\(UniMod:\d+\))';
  } 
  my %modindex;
  
  for my $m ( $eseq =~ /$regex/ig ) {
    my $index = index ($eseq, $m);
	my $em = $m;
	$em =~ s/\[/\\[/g;
	$em =~ s/\+/\\+/g;
	$em =~ s/\]/\\]/g;
	$em =~ s/\(/\\(/g;
	$em =~ s/\)/\\)/g;
	$em =~ s/-\[/\\[/g;
	$eseq =~ s/$em/$known_mods{$m}{'aa'}/;
	$modindex{$index} = $m;
	}
  
   @aa = split( //, $eseq );
	  for my $keymod (keys %modindex)
	{
		$aa[$keymod]= "$modindex{$keymod}";	
	} 
  
  my $bpep = '';
  my $ypep = '';

 #Nterminal modification
 $nterm_mod = $cterm_mod = '';
  if ($eseq =~ /^x/)
  {
	  $nterm_mod = shift @aa; 
  }
  if ($eseq =~ /x$/){
	  $cterm_mod = pop @aa;
  }
  
  my $last2nd = pop @aa;
  if ($last2nd =~ /-/)
   {  
		$cterm_mod = "-$cterm_mod";
   }
   else {
	   push @aa, $last2nd;
   }

 if ($nterm_mod && !$cterm_mod){
	 $bpep = $nterm_mod;
	 for ( my $idx = 0; $idx <= $#aa; $idx++ ) {
		 $bpep .= $aa[$idx];
		 $b{$idx + 1} = $bpep;
		 $ypep = $aa[$#aa - $idx] . $ypep;
		 if ($idx == ($#aa)){
			  $ypep = $nterm_mod.$ypep; 
		   }
		 $y{$idx + 1} = $ypep;
	 }
 } elsif ($cterm_mod && !$nterm_mod) {
	 $ypep = $cterm_mod || '';
	 for ( my $idx = 0; $idx <= $#aa; $idx++ ) {
		 $bpep .= $aa[$idx];
		 $b{$idx + 1} = $bpep;
		 $ypep = $aa[$#aa - $idx] . $ypep;
		 if ($idx == ($#aa)){
			 $ypep = $nterm_mod.$ypep;	
		 }
		 $y{$idx + 1} = $ypep;
		}
  }elsif ($nterm_mod && $cterm_mod){
	  $bpep = $nterm_mod || '';
	  $ypep = $cterm_mod || '';
	  for ( my $idx = 0; $idx <= $#aa; $idx++ ) {
			$bpep .= $aa[$idx];
			if ($idx == ($#aa)){
				$bpep = $bpep.$cterm_mod;
				}
			$b{$idx + 1} = $bpep;
			$ypep = $aa[$#aa - $idx] . $ypep;
			if ($idx == ($#aa)){
				$ypep = $nterm_mod.$ypep;	
			}
			$y{$idx + 1} = $ypep;
		}
	}else
	  {
		  for ( my $idx = 0; $idx <= $#aa; $idx++ ) {
			  $bpep .= $aa[$idx];
			  $b{$idx + 1} = $bpep;
			  $ypep = $aa[$#aa - $idx] . $ypep;
			  $y{$idx + 1} = $ypep;
		  }
	  }
  my %ions = ( b => \%b, y => \%y );
  return \%ions;
}

__DATA__
#
# SWATH File format info
#
# peakview 
# 0    Q1 [ 778.413 ]
# 1    Q3 [ 498.2579 ]
# 2    RT_detected [ 57.3 ]
# 3    protein_name [ 1/P0CG40 ]
# 4    isotype [ light ]
# 5    relative_intensity [ 7324.6 ]
# 6    stripped_sequence [ AAAAAAAAAAAAAAAASAGGK ]
# 7    modification_sequence [ AAAAAAAAAAAAAAAASAGGK ]
# 8    prec_z [ 2 ]
# 9    frg_type [ b ]
# 10   frg_z [ 1 ]
# 11   frg_nr [ 7 ]
# 12   iRT [ 57.3 ]
# 13   uniprot_id [ 1/P0CG40 ]
# 14   decoy [ FALSE ]
# 15   N [ 1 ]
# 16   confidence [ 1 ]
# 17   shared [ FALSE ]

# openswath
# 0    PrecursorMz [ 778.4129855 ]
# 1    ProductMz [ 498.26707303 ]
# 2    Tr_recalibrated [ 57.3 ]
# 3    transition_name [ 6_AAAAAAAAAAAAAAAASAGGK_2 ]
# 4    CE [ -1 ]
# 5    LibraryIntensity [ 7324.6 ]
# 6    transition_group_id [ 1_AAAAAAAAAAAAAAAASAGGK_2 ]
# 7    decoy [ 0 ]
# 8    PeptideSequence [ AAAAAAAAAAAAAAAASAGGK ]
# 9    ProteinName [ 1/P0CG40 ]
# 10   Annotation [ b7/-0.009,b14^2/-0.009,m10:16/-0.009 ]
# 11   FullUniModPeptideName [ AAAAAAAAAAAAAAAASAGGK ]
# 12   PrecursorCharge [ 2 ]
# 13   GroupLabel [ light ]
# 14   UniprotID [ 1/P0CG40 ]
# 15   FragmentType [ b ]
# 16   FragmentCharge [ 1 ]
# 17   FragmentSeriesNumber [ 7 ]
#
# Spectronaut
# 0     PrecursorMz [ 778.4129855 ]
# 1     PrecursorCharge [ 2 ]
# 2     ProteinName [ 1/P0CG40 ]
# 3     Organism [ Human ]
# 4     iRT Value [ 57.3 ]
# 5     PeptideSequence [ AAAAAAAAAAAAAAAASAGGK ]
# 6     PeptideModifiedSequence [ AAAAAAAAAAAAAAAASAGGK ]
# 7     FragmentType [ b ]
# 8     FragmentNumber [ 6 ]
# 9     Fragmentation [ b6 ]
# 10    Losses [  ]
# 11    LossNeutralMass [ 0 ]
# 12    ProductCharge [ 1 ]
# 13    ProductMz [ 427.22995924 ]
# 14    LibraryIntensity [ 8319.7 ]
# 15    ProteinId [ 1/P0CG40 ]
# 

# 0	ReferenceRun [ C_D160304_S251-Hela-2ug-2h_MSG_R01_T0 ]
# 1	PrecursorCharge [ 3 ]
# 2	IntModifiedPeptide [ _SLDTHVTK_ ]
# 3	ModifiedPeptide [ _SLDTHVTK_ ]
# 4	StrippedPeptide [ SLDTHVTK ]
# 5	iRT [ -35.19764 ]
# 6	BGSInferenceId [ Q14166 ]
# 7	IsProteotypic [ True ]
# 8	LabeledPeptide [ _SLDTHVTK_ ]
# 9	PrecursorMz [ 300.831025070542 ]
# 10	ReferenceRunQvalue [ 0.000844789494294673 ]
# 11	ReferenceRunMS1Response [ 9.9044E+07 ]
# 12	FragmentLossType [ noloss ]
# 13	FragmentNumber [ 3 ]
# 14	FragmentType [ y ]
# 15	FragmentCharge [ 1 ]
# 16	FragmentMz [ 347.228896545912 ]
# 17	RelativeIntensity [ 100 ]
# 18	ExcludeFromAssay [ False ]
# 19	Database [ sp ]
# 20	ProteinGroups [ Q14166 ]
# 21	UniProtIds [ Q14166 ]
# 22	Protein Name [ TTL12_HUMAN ]
# 23	ProteinDescription [ Tubulin--tyrosine ligase-like protein 12 ]
# 24	Organisms [ Homo sapiens ]
# 25	Genes [ TTLL12 ]
# 26	Protein Existence [ 1 ]
# 27	Sequence Version [ 2 ]
# 28	FASTAName [ uniprot_sprot_2014-12-11_HUMAN_ISOFORMS ]
#
# PROSIT
# 0     RelativeIntensity [ 0.047176161128631676 ]
# 1     FragmentMz [ 147.112804167 ]
# 2     ModifiedPeptide [ _AAAAAAALQAK_ ]
# 3     LabeledPeptide [ AAAAAAALQAK ]
# 4     StrippedPeptide [ AAAAAAALQAK ]
# 5     PrecursorCharge [ 2 ]
# 6     PrecursorMz [ 478.77981619566503 ]
# 7     iRT [ 9.95382308959961 ]
# 8     FragmentNumber [ 1 ]
# 9     FragmentType [ y ]
# 10    FragmentCharge [ 1 ]
# 11    FragmentLossType [ noloss ]
#
# CHROMAT_LIBRARY (Modified OSW)
# 0     transition_group_id [ HEELMLGDPC[57]LK+2 ]
# 1     transition_name [ HEELMLGDPC[57]LK+2_b3+2 ]
# 2     ProteinId [ sp|P07814|SYEP_HUMAN ]
# 3     PeptideSequence [ HEELMLGDPCLK ]
# 4     ModifiedPeptideSequence [ HEELMLGDPC[57]LK ]
# 5     FullPeptideName [ HEELMLGDPC[57]LK ]
# 6     RetentionTime [ 2220.0 ]
# 7     PrecursorMz [ 721.3443378086647 ]
# 8     PrecursorCharge [ 2 ]
# 9     ProductMz [ 396.151374562 ]
# 10    ProductCharge [ 1 ]
# 11    LibraryIntensity [ 1111.5 ]
# 12    FragmentIonType [ b ]
# 13    FragmentSeriesNumber [ 3 ]
# 14    IsDecoy [ 0 ]
# 15    quantifying_transition [ 1 ]
# 
#EOF
