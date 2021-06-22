#!/usr/bin/perl -w

#===============================================================================
#        ________             ____________            
#   ________  __ )_________  ___  ___/__(_)__________ 
#   _  _ \_  __  |  __ \_  |/_/____ \__  /___  /_  _ \
#   /  __/  /_/ // /_/ /_>  < ____/ /_  / __  /_/  __/
#   \___//_____/ \____//_/|_| /____/ /_/  _____/\___/ 
#                                                  
#   eBoxSize - calculation of the box size for ligand docking
#
#   This software is distributed WITHOUT ANY WARRANTY (but with best wishes)
#
#   Report bugs and issues to wfeinstein@lsu.edu
#                             michal@brylinski.org
#
#   Computational Systems Biology Group
#   Department of Biological Sciences
#   Center for Computation & Technology
#   Louisiana State University
#   407 Choppin Hall, Baton Rouge, LA 70803, USA
#
#   http://www.brylinski.org
#
#===============================================================================

use strict;
use constant GY_BOX_RATIO => 0.23;

sub GeometryCenter{
	my $atoms=$_[0];
	my @locations=@{$_[1]};
	my @xyzloc=();
	my $x=0;
	my $y=0;
	my $z=0;

	for (my $i=0; $i< $atoms; $i++){
	
		$x += $locations[$i*3+0]*1;
		$y += $locations[$i*3+1]*1;
		$z += $locations[$i*3+2]*1;
	}
	$xyzloc[0]=$x/$atoms;
	$xyzloc[1]=$y/$atoms;
	$xyzloc[2]=$z/$atoms;

	return @xyzloc; 
}
sub Usage{
 print "------------------------------------------------------------\n";
 print "                         eBoxSize\n";
 print "                        version 1.1\n";
 print "       calculation of the box size for ligand docking\n\n";
 print "       report bugs and issues to wfeinstein\@lsu.edu\n";
 print "                                 michal\@brylinski.org\n";
 print "------------------------------------------------------------\n\n";
 
 print "eBoxSize.pl <ligand file in SDF, PDBQT or MOL2>\n";
 die "\n";
}

if ($#ARGV != 0)
{
	Usage();  
}

my $input = $ARGV[0];
my $format;

if ($input =~ /\.sdf$/i){
	$format = "sdf";
}
elsif($input =~ /\.pdbqt$/i){
	$format = "pdbqt";
}
elsif($input =~ /\.mol2$/i){
	$format = "mol2";
}
else {
	Usage();
}

my ($line,@line, @locs);
my $atoms = 0;
my $i=0;

open (INPUT,$input) || die "cant open $input: $!\n";

# .sdf input format
if ($format eq 'sdf'){
	for ($i=0; $i<4; $i++){
		$line = <INPUT>;
        }
        my $ct = substr($line, 0, 3);
        
        for ($i = 0; $i < $ct; $i++){	
		$line = <INPUT>;            
                if (substr($line, 31, 1) ne 'H'){
			@line=split(' ', $line);
			push(@locs, $line[0]); 
			push(@locs, $line[1]);
			push(@locs, $line[2]);
                        $atoms++;
                }
        }
}

# .mol2 input format
elsif ($format eq 'mol2'){
	for ($i=0; $i<3; $i++){
		$line = <INPUT>;
	}
        my $ct = substr($line, 0, 3);
       
	while ($line = <INPUT>){
            if ($line =~ /ATOM/i){
                for ($i = 0; $i < $ct; $i++){
         		$line = <INPUT>;
                        if (substr($line, 9, 1) ne 'H'){
                      	  @line=split(' ', $line);
       		 	        push(@locs, $line[2]); 
				push(@locs, $line[3]);
				push(@locs, $line[4]);
	                        $atoms++;
                        }
                        
                }
                 last;          
             }	
	}
}        

# .pdbqt input format
elsif ($format eq 'pdbqt'){
	my @lines = <INPUT>;
        my @atomlines = grep(/ATOM  |HETATM/, @lines);
            
        foreach $line (@atomlines) {
             if(substr($line, 13, 1) ne 'H'){
       		   my $loc = substr($line, 30, 8);
                   push(@locs, $loc);
                   $loc = substr($line, 38, 8);
                   push(@locs, $loc);
                   $loc = substr($line, 46, 8);
                   push(@locs, $loc);
                   $atoms++;
             }
         }          
}

close(INPUT); 

my @geo_center = GeometryCenter($atoms, \@locs);  
my $dis;

for ($i=0; $i< $atoms; $i++){

	$dis += ($locs[$i*3+0]-$geo_center[0])**2 + ($locs[$i*3+1]-$geo_center[1])**2 + ($locs[$i*3+2]-$geo_center[2])**2;
}

my $gyration_g= sqrt($dis/$atoms)/GY_BOX_RATIO;

printf("%.3f\n", $gyration_g);

exit;