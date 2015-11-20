INSTALLATION
============

Requires the Portable Cray Bioinformatics Library, available at 
http://cbl.sf.net, for alignment of top scorers.

Then simply edit the makefile for the installation directory, and edit the 
#define NUM_SPE 6
statement in sw_db.c for the number of spe's on your system. Then do

$ make
$ sudo make install


USAGE
=====

sw_db [options] <scoring matrix file> <query fasta file> [db fasta file]
      -a output alignments (default = top scores only)
      -e extend gap penalty (eg), default = 1
         1 gap is scored as og + eg
      -n number of top scores to report (default = 10)
      -o open gap penalty (og), default = 8
      -h help 
      <query fasta file> a single query
      [db fasta file] may be specified or read from stdin
	      
example: sw_db -a -e 2 -o 3 pam_1 query nt

See the included pam_1 example scoring matrix file, where comments have a 
"#" in the first column, and the last column of letters is optional.
