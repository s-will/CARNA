#!/usr/bin/perl -w


=head1 NAME

carna.pl

=head1 SYNOPSIS

carna.pl [options] input.fa

=head1 DESCRIPTION

DEPRECATED: rather call carna binary directly (one input files per sequence)

Call carna from a fasta file with two sequences. This script is a
front end for the constraint RNA alignment program carna.

=head1 OPTIONS

=over 8

=item  B<-h, --help>
    
Brief help message

=item  B<-m, ---man>

Full documentation

=item  B<-v, --verbose>

Verbose

=item  B<-vv, --moreverbose>

More Verbose

=item  B<--norun>

Don't run carna, but generate input for carna. Show call with -v option.

=item  B<--args>

Arguments passed to carna. Call carna --help or --args '--help' for
the available options of carna.

=back

=head1 INPUT FORMAT

input.fa is a fasta file containing two sequences. The sequences can
be annotated by dot-bracket strings describing crossing secondary
structure (using different bracket symbols for crossing base pairs).
If no structure annotation is provided, we compute a dot plot using RNAfold -p.

The file can look like this:

 >nameA
 ACGUACGGACGGCAGUGC
 >nameB
 ACGAGGACGCUAGCGAGCGC

or (with structure annotation)

 >nameA
 ACGUACGGACGGCAGUGC
 .((..[[..)).....]]
 >nameB
 ACGAGGACGCUAGCGAGCGC
 .[[.((.]].)......)..

The following bracket pairs are supported: (), [], {}, <>, Aa, Bb, Cc, Dd

Additionally the script supports anchor and structural constraints that are 
specified and have the same semantic as in mlocarna (of the LocARNA package).

Here is example input with anchor and structure constraints:

 >nameA
 GCCAUACGGCAUAC
 (((.[[.))).]]. #S
 ....AA........ #1
 ....12........ #2
 >nameB
 GGUGACCGCCAACAC
 ((([[..)))..]]. #S
 ...AA.......... #1
 ...12.......... #2

First, note that the tag '#S' changes the semantics of the given
structure. A structure with '#S' tag is passed as constraint to
RNAfold and not interpreted as fixed structure.

The lines with #<i> tag give names to sequence positions. One can
specify lines with tags #1,...,#N for arbitrary N. The name of a
sequence position is the concatenation of the corresponding character
in lines with tags #1,...,#N. In the example, the position names are
"A1" and "A2".  The anchor constraint tells that positions with equal
name have to be aligned to each other.

=head1 AUTHORS

Alessandro Dal Palu, Mathias Moehl, Sebastian Will

=head1 BUGS

Once you tell us, we will hunt them.

=cut


use strict;

use FindBin;

##------------------------------------------------------------
## options
use Getopt::Long;
use Pod::Usage;

my $help;
my $man;
#my $quiet;
my $verbose;
my $moreverbose;

my $norun;

my $carnaArgs="";

my $opt_sym;

## Getopt::Long::Configure("no_ignore_case");

GetOptions(
    "v|verbose" => \$verbose,
    "vv|moreverbose" => \$moreverbose,
    #"quiet" => \$quiet,   
    "sym" => \$opt_sym,
    "help"=> \$help,
    "man" => \$man,
    "norun" => \$norun,
    "args=s" => \$carnaArgs
    ) || pod2usage(2);

pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

if ($#ARGV != 0) {
    print STDERR "ERROR: No input file name provided.\n";
    pod2usage(-exitstatus => -1, -verbose => 1);
}
my $inputfilename=$ARGV[0];

if ($moreverbose) {$verbose=1;}

my $tmpsuf="~$$.tmp~";


## trap the control-C signal and clean up
sub INT_handler {
    print "Caught Ctrl-C. Cleanup and Exit.";
    cleanup();
    exit(0);
}

$SIG{'INT'} = 'INT_handler';


## -----------------------------------------------------------
## subs


## copied from the module MLocarna.pm
########################################
## make_normalized_seqname($name, @names list of existing names)
##
## generate a sequence name from $name that
## has at most a length of 16
## and
## does not already exist in @names
## 
########################################
sub make_normalized_seqname {
    my ($name, @names) = @_;
    
    my $maxlen=16;

    chomp $name;
    $name =~ s/[^a-zA-Z\d]/_/g;
    $name = substr $name,0,$maxlen;
    
    my $i=1;
    while (grep /^$name$/, @names) {
	my $arity = int(log($i)/log(10))+1;
	if ($arity >= 4) { # arbitrary limit
	    die "Could not generate unique name";
	}
	
	$name = substr $name,0,$maxlen-$arity-1;
	$name = sprintf("%s_%0$arity"."d",$name,$i);
	$i++;
    }
    return $name;
}

## copied from the module MLocarna.pm
########################################
## check_constraints_in_fasta($fasta)
##
## Checks validity of constraint description in fasta list of hash
## representation.
##
## $fasta list of hashs representation
##
## die with error message when constraint annotation is invalid
##
########################################
sub check_constraints_in_fasta {
    my $fasta=shift;

    my $numseqconstraints = -1;

    for my $seq (@$fasta) {
	if (exists $seq->{"ANNO#S"}) {
	    if (length($seq->{"ANNO#S"}) != length($seq->{seq})) {
		die "Structure constraint length unequals sequence length for "
		    .$seq->{name}.".";
	    }
	}
	
	my $i=1;
	while (exists $seq->{"ANNO#$i"}) {
	    if (length($seq->{"ANNO#$i"}) != length($seq->{seq})) {
		die "Sequence constraint length unequals sequence length for "
		    .$seq->{name}.".";
	    }
	    $i++;
	}
	if ($numseqconstraints == -1) {$numseqconstraints=$i-1;}
	
	if ($numseqconstraints != $i-1) {
	    die "Bad sequence constraints for sequence "
		.$seq->{name}.".";
	}
    }
}

## copied from the module MLocarna.pm
########################################
## sequence_constraint_string($seq hash ref)
##
## Collect sequence constraint string for $seq entry of a fasta list of hashs
##
## returns string "<c1>#<c2>#...<cN>", where ci is the constraints description line i
##
########################################+
sub sequence_constraint_string {
    my $seq=shift;
    my $i=1;
    
    my $str="";
    while (exists $seq->{"ANNO#$i"}) {
	$str .= $seq->{"ANNO#$i"}."#";
	$i++;
    }
    $str =~ s/\#$//;
    
    return $str;
}


## copied from the module MLocarna.pm and MODIFIED to support fix input structure!!! (see below)
########################################
## read_fasta($filename)
## ------------------------------
##
## Read a simplified fasta file. If filename ends with "*.gz" the file
## is automatically uncompressed by gunzip.  Dies if file is not
## readable or gunzippable.  The result is returned as list of the
## sequences.
##
##
## Each sequence is encoded as hash of features
##   name: id of sequence
##   descr: description of the sequence
##   nname: normalized id of the sequence 
##       (generated from id, 
##        can be used as a filename which is unique for the set of
##        sequences in the fasta file)
##   seq:  sequence string
##
## supports special annotation strings as used by locarna
## lines ending in #(\S+) are recognized and 
## returned as features "ANNO#$1" of the sequence.
## 
## return reference to list of hash representation
##
########################################
sub read_fasta {
    my $filename = shift;

    my $fh;
    
    if ($filename =~ /\.gz$/) {
	open($fh,"gunzip -c $filename |") || die "Cannot read or uncompress file $filename.";
    } else {
	open($fh,$filename) || die "Cannot read file $filename.";
    }
 
    my @fasta = ();
    my @names = ();
    
    my $line=<$fh>;
    while(defined($line)) {
	if ($line=~/^>\s*(\S+)\s*(.*)/) {
	    my $name=$1;
	    my $description=$2;
	    
	    my $nname = make_normalized_seqname $name,@names;
	    
	    push @names, $nname;
	    
	    my $seq = { name  => $name,
			nname => $nname, 
			descr => $description };
	    
	    while (defined($line=<$fh>) && ($line !~ /^>/)) {
		chomp $line;
		$line =~ s/\s+//g;
		
		if  ($line =~ /(.+)\#(.+)/) {
		    $seq->{"ANNO\#$2"} .= $1;
		}
		# line is not a sequence and no constraint line --- this case is added to version from MLocarna.pm
		# assume the line encodes a structure
		elsif ($line !~ /^[A-Za-z-]*$/) { 
		    $seq->{"fixstructure"} .= $line;
		}
		    else {
		    $seq->{seq} .= $line;
		}
	    }
	    
	    push @fasta, $seq;
	} else {
	    $line=<$fh>;
	}
    }
    
    close $fh;
    
    return \@fasta;
}


# sub parse_mfasta {
#     my ($file) = @_;
#     my %mfasta;
#     my @names;
    
#     local *PMF_IN;
    
#     open(PMF_IN,$file) || die "Cannot read mfasta file $file\n";

#     my $line=<PMF_IN>;
#     while(defined($line)) {
#         if ($line=~/^>\s*(\S+)/) {
#             my $name=$1;
#             chomp $name;
#             $name =~ s/[^a-zA-Z\d]/_/g;
            
#             if (exists $mfasta{$name}) {
#                 print STDERR "Duplicate name in mfasta: $name\n";
#                 exit -1;
#             }
	    
# 	    push @names,$name;

#             while (defined($line=<PMF_IN>) && ($line !~ /^>/)) {
#                 chomp $line;
#                 $line =~ s/^\s+//; # strip leading and ending blanks
#                 $line =~ s/\s+$//;
# 		if ($line =~ /^[.(){}\[\]<>AaBbCcDd]+$/) {
# 		    $mfasta{$name."#S"} .= $line;
# 		} else {
# 		    if($line !~ /^[ACGTUNacgtun]*$/)  {
# 			print STDERR "WARNING: while parsing mfasta: unrecognized symbols in $line\n";
# 		    }
# 		    $mfasta{$name} .= $line;
# 		} 
# 	    }
#         }
#     }
#     return (\%mfasta,\@names);
# }


sub parse_bracket_structure_single {
    my ($str,$open,$close,$str_array_ref)=@_;
    
    my @str_array = @{ $str_array_ref };
    
    my @stack;
    
    for (my $i=0; $i<length($str); $i++) {
	my $c=substr $str,$i,1;
	
	if ($c eq $open) {
	    push @stack,$i;
	} elsif ($c eq $close) {
	    my $j=pop @stack;
	    $str_array[$i]=$j;
	    $str_array[$j]=$i;
	}
    }

    return @str_array;
}


## @returns array for structure
## entry per sequence position: gives paired position or -1
sub parse_bracket_structure {
    my ($str)=@_;
    my @str_array;

    for (my $i=0; $i<length($str); $i++) {
	$str_array[$i]=-1;
    }

    @str_array=parse_bracket_structure_single $str,"(",")",\@str_array;
    @str_array=parse_bracket_structure_single $str,"{","}",\@str_array;
    @str_array=parse_bracket_structure_single $str,"[","]",\@str_array;
    @str_array=parse_bracket_structure_single $str,"<",">",\@str_array;
    @str_array=parse_bracket_structure_single $str,"A","a",\@str_array;
    @str_array=parse_bracket_structure_single $str,"B","b",\@str_array;
    @str_array=parse_bracket_structure_single $str,"C","c",\@str_array;
    @str_array=parse_bracket_structure_single $str,"D","d",\@str_array;
    
    return @str_array;
}


sub convert_fix_structure_to_pp {
    my ($ppfilename,$name,$seq,$str) = @_;
    
    my @str=parse_bracket_structure($str);
    
    local *OUT;

    open(OUT,">$ppfilename") || die "Cannot write $!";

    print OUT "\n";
    print OUT "$name\t$seq\n";
    print OUT "\n";
    print OUT "\n";
    print OUT "#\n";
    
    for (my $i=0; $i<=$#str; $i++) {
	if ($str[$i]>$i) {
	    print OUT ($i+1)." ".($str[$i]+1)." 1.0\n";
	}
    }
}

sub cleanup {
    unlink "1.pp.$tmpsuf";
    unlink "2.pp.$tmpsuf";
}

## ------------------------------------------------------------
## main part

my $mfasta_ref = read_fasta($inputfilename);

check_constraints_in_fasta($mfasta_ref);

my @mfasta = @{ $mfasta_ref };

if ($verbose) {
    print "Read from $inputfilename:\n";
    foreach my $entry (@mfasta) {
	foreach my $k (keys %{ $entry }) {
	    print "$entry->{name} $k $entry->{$k} \n";
	}
    }
    print "\n";
}

if (2!=@mfasta) {
    print STDERR "ERROR: Need exactly two sequences in input.\n";
    exit(-1);
}

my $i=1;
foreach my $entry (@mfasta) {
    if (exists $entry->{"fixstructure"}) {
	convert_fix_structure_to_pp("$i.pp.$tmpsuf",$entry->{name},$entry->{seq},$entry->{"fixstructure"});
    } else {
	my $cmd="RNAfold -p";
	if (exists $entry->{"ANNO#S"}) {
	    $cmd .= " -C\"".($entry->{"ANNO#S"})."\"";
	}
	system("printf \"$entry->{seq}\" | $cmd")==0 || die "Could not 'RNAfold -p' sequence $i\n";
	rename "dot.ps", "$i.pp.$tmpsuf";
	unlink "rna.ps";
    }
    
    if ($moreverbose) {
	print "Wrote to $i.pp.$tmpsuf:\n";
	system "cat $i.pp.$tmpsuf";
	print "EOF --- $i.pp.$tmpsuf\n";
    }

    
    $i++;
}

my $ppfile1="1.pp.$tmpsuf";
my $ppfile2="2.pp.$tmpsuf";
my $constraintString1 = sequence_constraint_string($mfasta[0]);
my $constraintString2 = sequence_constraint_string($mfasta[1]);

if($opt_sym) {
    #swap($ppfile1,$ppfile2);
    my $tmp=$ppfile1; $ppfile1=$ppfile2; $ppfile2=$tmp;
    #swap($constraintString1,$constraintString2);
    $tmp=$constraintString1; $constraintString1=$constraintString2;
    $constraintString2=$tmp;
 }

my $cmd = "$FindBin::Bin/carna $carnaArgs $ppfile1 $ppfile2";


if (($constraintString1 ne "")
    &&
    ($constraintString2 ne "")) {
    $cmd .= " --anchorA=\"$constraintString1\"";
    $cmd .= " --anchorB=\"$constraintString2\"";
}

$cmd .= " -v" if $moreverbose;

print "Call CARNA: $cmd\n" if $verbose;

if (!$norun) {
    system("time -p $cmd");
    cleanup();
}

## ------------------------------------------------------------
