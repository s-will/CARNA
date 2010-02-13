#!/usr/bin/perl -w

=head1 NAME

cRNAalign.pl

=head1 SYNOPSIS

cRNAalign.pl [options] input.fa

Options:

=over 1

=item  B<-h, --help>                    Brief help message

=item  B<-m, ---man>                    Full documentation

=item  B<-v, --verbose>                 Verbose

=item  B<-q, --quiet>                   Quiet

=item  B<--args>                        Arguments passed to RNAalignment


=back

=head1 DESCRIPTION

cRNAalign.pl is a front end for the constraint RNA alignment program RNAalignment.
Its main use is to provide an easy to use interface to RNAalignment.

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

=cut


use strict;

##------------------------------------------------------------
## options
use Getopt::Long;
use Pod::Usage;

my $help;
my $man;
my $quiet;
my $verbose;
my $RNAalignmentArgs="";

## Getopt::Long::Configure("no_ignore_case");

GetOptions(	   
    "verbose" => \$verbose,
    "quiet" => \$quiet,   
    "help"=> \$help,
    "man" => \$man,
    "args=s" => \$RNAalignmentArgs
    ) || pod2usage(2);

pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

if ($#ARGV != 0) {
    print STDERR "ERROR: No input file name provided.\n";
    pod2usage(-exitstatus => -1, -verbose => 1);
}
my $inputfilename=$ARGV[0];


my $tmpsuf="~$$.tmp~";

## -----------------------------------------------------------
## subs

sub parse_mfasta {
    my ($file) = @_;
    my %mfasta;
    
    local *PMF_IN;
    
    open(PMF_IN,$file) || die "Cannot read mfasta file $file\n";

    my $line=<PMF_IN>;
    while(defined($line)) {
        if ($line=~/^>\s*(\S+)/) {
            my $name=$1;
            chomp $name;
            $name =~ s/[^a-zA-Z\d]/_/g;
            
            if (exists $mfasta{$name}) {
                print STDERR "Duplicate name in mfasta: $name\n";
                exit -1;
            }
            while (defined($line=<PMF_IN>) && ($line !~ /^>/)) {
                chomp $line;
                $line =~ s/^\s+//; # strip leading and ending blanks
                $line =~ s/\s+$//;
		if ($line =~ /^[.(){}\[\]<>AaBbCcDd]+$/) {
		    $mfasta{$name."#S"} .= $line;
		} else {
		    if($line !~ /^[ACGTUNacgtun]*$/)  {
			print STDERR "WARNING: while parsing mfasta: unrecognized symbols in $line\n";
		    }
		    $mfasta{$name} .= $line;
		} 
	    }
        }
    }
    return %mfasta;
}


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
	    print OUT "$i $str[$i] 1.0\n";
	}
    }
}

sub cleanup {
    unlink "1.pp.$tmpsuf";
    unlink "2.pp.$tmpsuf";
}

## ------------------------------------------------------------
## main part

my %mfasta=parse_mfasta($inputfilename);

my @names=grep !/#S$/, keys %mfasta;

if ($#names!=1) {
    print STDERR "ERROR: Need exactly two sequences in input.\n";
    exit(-1);
}

my $i=1;
foreach my $name (@names) {
    print "$name\n";
    if (exists $mfasta{$name."#S"}) {
	convert_fix_structure_to_pp("$i.pp.$tmpsuf",$name,$mfasta{$name},$mfasta{$name."#S"});
    } else {
	print STDERR "Computing structure not implemented yet. Please give structures.\n";
	exit(-1);
    }
    $i++;
}


my $cmd = "src/RNAalignment $RNAalignmentArgs "."1.pp.$tmpsuf 2.pp.$tmpsuf";
print "CALL: $cmd\n";

if ( system($cmd) !=0 ) {
    print STDERR "ERROR: CALL $cmd failed.\n"; exit(-1); 
}

cleanup();

## ------------------------------------------------------------
