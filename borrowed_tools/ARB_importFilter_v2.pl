#!/usr/bin/env perl
use warnings;
#First 2 lines modified by Chris according to Doug
#Script to reformat spreadsheet and generate custom import filter for import into ARB

system ('cls');

print "\n*********************************
This script reformats spreadsheet data for import into ARB and \ngenerates a corresponding custom import filter 

Written by Brian Oakley.  Please contact with questions or comments:
b.b.c.oakley at warwick.ac.uk
*********************************\n

Press enter to continue...";
<STDIN>;

print "\nTo use,\n
Your data should be formatted in a spreadsheet with any number 
of fields (columns) associated with your sequence data.  
\nThe only critical things are:

1)  Each field (aka column) must have a column heading - this is what 
the field is called in your import filter and ultimately your ARB database.

2)  The first field must be a unique identifier for your sequences and MUST be 
called 'name' or 'full_name' (with no caps) in your spreadsheet.  
The ARB team generally frowns upon using your own name field, but as long as 
names are relatively simple (e.g. 'S1_25' for 'Site 1 clone 25') 
it should be fine and will avoid nonsensical names imposed by ARB.  
Avoid special characters other than underscore.
 
3)  The last field in your database MUST be the sequence data.  Gaps can be 
included and your alignment will be kept as is.  
\nIf you have an existing ARB database and don't want to change the sequence data
but just want to add or update some other fields, just fill in some dummy data 
(e.g. 'ACC') for the last field and then ignore it when you merge with your existing database. 


PLEASE NOTE:  All output files are saved in the same directory as this script.

*******************************************************

\nEnter name of input file to convert when ready...";

chomp (my $input_file = <STDIN>);

open (INPUT, $input_file); 
open (OUT1, ">new_arb_database");
open (OUT2, ">new_import_filter.ift");



@file = <INPUT>; 
$column_headings = shift(@file); 
$column_headings =~ s/ //g;
$column_headings =~ s/\n//g;

@column_headings = split("\t", $column_headings);
$num_of_columns = scalar (@column_headings); 



foreach $line (@file) { 
$line =~ s/ //g;
$line =~ s/\n//g;
	@tabs =split("\t",$line);


	print OUT1 ("BEGIN\n");
$x = -1;
	while ($x < $num_of_columns-2) {
	$x += 1;	
		print OUT1 ("$column_headings[$x]=$tabs[$x]\n");
		
}

	$sequence_data = $tabs[$num_of_columns-1];
	print OUT1 ("warning\n$sequence_data\n");
	print OUT1 ("END\n\n");
}

pop(@column_headings);

print OUT2 ("AUTODETECT\t\"BEGIN\"\n\nBEGIN\t\"BEGIN*\"\n\n");

$y = -1;
	while ($y < $num_of_columns-2) {
	$y += 1;	
		print OUT2 ("MATCH\t\"$column_headings[$y]\\=*\"\n");
		print OUT2 ("\tSRT \"*\\==\"\n");
		print OUT2 ("\tWRITE \"$column_headings[$y]\"\n\n");
}

print OUT2 ("SEQUENCEAFTER\t\"warning\"\nSEQUENCESRT\t\"*\\==\"\nSEQUENCEEND\t\"END\"\n\nEND\t\"END\"\n");


print "\nDone.  Output files created successfully:

1)  'new_arb_database.txt' which contains reformatted sequence data and 
associated information for import into ARB.

2)  'new_import_filter.ift' which is a new import filter custom-made 
for your data.  
\nThis file needs to live in the 'import' directory 
of your ARB installation.  
This is usually at 'usr/arb/lib/import' but may differ depending on how 
ARB is installed.  You will probably need to change permissions to allow 
this file to be copied into the directory.  
\n\n";

close OUT1;
close OUT2;


