# DIALib-QC
DIALib-QC, an assessment tool for the systematic evaluation of spectral assay libraries to improved data independent acquisition analysis.

DIALib-QC v1.0 Tool to assess the spectral ion libraries of PEAKVIEW,SPECTRONAUT or OPENSWATH formats.

1)System Requirements:

	- Perl 5 version 16 on Linux systems 
	- Perl Statistics::Regression module dependency (recommended) 
	- Perl Statistics::R module dependency (recommended)
	- R libraries: ggplot2, scales, ggpubr
	
	Tool successfully tested on 
	- Windows 10 Enterprise
	- Oracle Linux Server release 7.7 VERSION="7.7"

2)Installation guide

	DIALib-QC v1.0 is a tool that allows users to analyze DIA libraries in a variety of formats (PeakView, OpenSWATH, Spectonaut and Prosit) based on a host of properties. 
	These include charge, length, retention time, m/z, complexity, etc. This tool can be run online at PeptideAtlas(https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/AssessDIALibrary) 
	for moderately-sized libraries (up to 200 MB). The tool provides an option to download the result as tabular format and a single pdf file containing plots of library characterstics. It is written in the perl programming language, so you will need to have perl 
	installed on your Linux system in order to run it. For ploting feature, R installation along with few libraries is a pre-requisite.
	 
	You can downloaded DIALib-QC tool at http://www.swathatlas.org/DIALibQC.php link and run the tool locally on most Linux systems.
	
(A) Installing Instructions LINUX

	There is one non-core module dependency, Statistics::Regression, you can read more about it at (https://metacpan.org/pod/Statistics::Regression).
	The tool will work without it, but you will lose the test comparing RT values of +2 and +3 ions of the same modified peptide, a metric of the retention time 
	consistency of the ion library being analyzed. 
	Instructions for perl module installation are as follows:
	
	You can install the Statistics::Regression module either 
	
	i) manually or ii) via CPAN, as shown below Optional components in [ square brackets ], 
	sudo allows user to install as if root (probably required for system-wide installations), 
	while --PREFIX allows user to designate an alternative library location on which they have write permissions. 

		i) Manual installation of Statistics::Regression
		wget "https://cpan.metacpan.org/authors/id/I/IA/IAWELCH/Statistics-Regression-0.53.tar.gz"
		tar xvfz Statistics-Regression-0.53.tar.gz
		cd Statistics-Regression-0.53/
		perl Makefile.PL [ --PREFIX="/path/to/user_lib_dir/" ]
		make
		[sudo] make install

		If doing a custom PREFIX install, you'll need to insure that the library can be found, e.g.
		[ export PERL5LIB=$PERL5LIB:/path_to_user_lib_dir/share/perl5/ ]


		ii) CPAN installation of Statistics::Regression
		[sudo] perl -MCPAN -e shell
		cpan[2]> install Statistics::Regression
		exit

	You can install the Statistics::R module as
	
		i) perl -MCPAN -e 'install Statistics::R'
	
	You also need to install ggplot2 and ggpubr in R consule
		
		install.packages("ggplot2")
		install.packages("ggpubr"）
		

(B) Usage after installation: Assessing Library

	perl assess_swathlib.pl --help

	Starting run

  	Missing required option ion_library

	Usage: assess_swathlib.pl --ion_library libfile [ --peptide_file pepfile  --output output_file ]

	# Required
	ion_library    Path to ion library file in Peakview/OpenSWATH/Spectronaut/Prosit format (required)

	# Commonly used
	assess_massdiff   Compare library masses to theoretical monoisotopic values
	c, coldefs        print out column definitions
  	f, full_stats     Print file (output_file.fullstats)
  	h, help           Print this usage information and exit
	invert_output     Invert output, 3 columns by X rows
	o, output_file    File (base) to which to write data; given MyOut, generates
        	          MyOut.QC, MyOut.fullstats, etc. If not specified, library name
                	  is used.
	peptide_file      Path to peptide digest file, format is protein(s)   peptide
	swath_file        Path to file of SWATH definitions, format min_m/z   max_m/z

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
	print_mismatch    Print out ions whose m/z values differ significantly from
        	          theoretical
	print_ux          Print unexplained proteins (seen in ion library but not
        	          in peptide_file if provided
	rt_stats          Print file of +2/+3 RTs for analysis to output_file.RT,
        	          only if enough for regression
	skip_decoys       Skip DECOY entries when computing statistics
	v, verbose        Verbose output mode
	write             Write  MGF file with info from each spectrum (rare)
	zero_pad          Print intermediate values=0 with full_stats.


(C)  Usage after installation: Generating plots
	
	Usage: DIALib-QC_RPlot.pl library.fullstats library.RT
	or 
	Usage: DIALib-QC_RPlot.pl library.fullstats

(D) License
        Copyright (C) 2020 Moritz Lab, Institute for Systems Biology

DIALib-QC v1.0 program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.


You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
