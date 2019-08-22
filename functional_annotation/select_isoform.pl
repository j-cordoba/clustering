#!/usr/bin/env perl

use Modern::Perl '2011';
use autodie;

use Tie::IxHash;
use Smart::Comments;
use List::AllUtils qw(sum reduce max);
use File::Basename;
use Path::Class 'file';

unless (@ARGV == 3) {
die <<"EOT";
Usage: $0 <infile.aa.pub.fasta> <expression.table> <infile.pubids>
EOT
}

my $seq_for = read_fasta($protein_fasta);
my $expr_for = read_expr($expression_file);
my $class_for = read_table($evigene_file);

# Obtain gene identifiers from fasta file
my @ids = keys %$seq_for;

# Prepare output file to be written
#my ($basename, $dir, $suffix) = fileparse($protein_fasta, qr{ \.[^.]* }xms );
my $outfile = 'outfile.txt';    #file($dir, $basename . '_best' . $suffix);
open my $out, '>' , $outfile;

my %rank_for_class = (
 main              => 1, # primary transcript with longest high quality cds with known alternatives
 maina2            => 2, # primary transcript with longest high quality cds with known alterantive with an amino acid similiarity (but not cds) of 98%
 althi1            => 3, # higher identity than hi
 althi             => 4, # 98% cutoff
 althim            => 5,
 altmid            => 6, # lower identity matches, genome mapping by DG indicates these may be paralogs/homeologs
 altmida2          => 6, # lower nucleotide identity matches with high protein identity (98%)
 noclass           => 7, # primary transcript without known alternatives
 noclassa2         => 8, # primary transcript without known alternatives, but with a 98% similarity of amino acid to another transcript (but different cds - meaning it could be a different loci/duplication in genome)
 altmidfrag        => 9, # lower identity partial matches, may be paralogs/homeologs
 altmidfraga2      => 9, # lower nucleotide identity partial matches with high protein identity
 dropalthi         => 10,
 dropalthi1        => 10,
 dropaltmid        => 10,
 dropaltmida2      => 10,
 dropaltmidfrag    => 10,
 dropaltmidfraga2  => 10,
 dropnoclass       => 10,
 dropnoclassa2     => 10,
);

# calculate the longest, the most expressed and the best class
GENE:
for my $gene (@ids) {

        # say $gene;
        my $iso_len = $seq_for->{$gene};
        my $iso_expr = $expr_for->{$gene};
        my $iso_class = $class_for->{$gene};

        unless ($iso_expr) { # Issue with IDs equivalence between FINAL_aa_pub and FINAL_tr_pub 607 gene missing!!
                ### Skipping gene: $gene
                next GENE;
        }

        # Length comparison if possible, otherwise best class
        my @len_sort = sort {
                $iso_len->{$b} <=> $iso_len->{$a} # Numerical sort of length values inside gene -> isoforms hash
             || $rank_for_class{$iso_class->{$a}} <=> $rank_for_class{$iso_class->{$b}}
        } keys %$iso_len;

        # Expression comparison if possible, otherwise best class
        my @expr_sort = sort {
                $iso_expr->{$b} <=> $iso_expr->{$a} # Numerical sort of expression values inside gene -> isoforms hash
             || $rank_for_class{$iso_class->{$a}} <=> $rank_for_class{$iso_class->{$b}}
        } keys %$iso_expr;

        # Best class acording to $rank_for_class
        my @class_sort = sort {
                $rank_for_class{$iso_class->{$a}} <=> $rank_for_class{$iso_class->{$b}}
        } keys %$iso_class;


        my $len_id   = shift @len_sort;
        my $expr_id  = shift @expr_sort;
        my $class_id = shift @class_sort;
        ## $len_id
        ## $expr_id
        ## $class_id

        # Computing ratio between len_id and expr_id length
        #my $ratio = $iso_len->{$len_id} / $iso_len->{$expr_id};
        ## $ratio

        # Reading the class asigned to isoform with longest sequence $len_id    
        my $class_len = $iso_class->{$len_id};
        ## $class_len

        # Reading the class asigned to isoform most expressed $expr_id  
        my $class_expr = $iso_class->{$expr_id};
        ## $class_expr

        # Reading the class asigned to isoform most expressed $expr_id  
        my $class_class = $iso_class->{$class_id};
        ## $class_class


        # Writing results to a table
        say {$out} join qq{\t}, $gene,  $len_id, $class_len,  $expr_id, $class_expr, $class_id, $class_class;
}


close $out;


sub read_fasta {
        my $infile = shift;
        ### Reading input file: $infile

        open my $in, '<', $infile;

        my $seq_id;
        my $seq;
        tie my %seq_for, 'Tie::IxHash';         # preserve original seq order
                my $gene;

        LINE:
        while (my $line = <$in>) {
                chomp $line;

                        next LINE if $line =~ m/ ^ \s* \z /xms;                 # skip empty lines
                        next LINE if $line =~ m/ ^ \# /xms;                     # skip comment lines


                # at each '>' char...
                if (substr($line, 0, 1) eq '>') {

                        # add current seq to hash (if any)
                        if ($seq) {
                                # Create a hash of hash [gene -> isoforms -> values]
                                $seq_for{$gene}{$seq_id} = length $seq;
                                $seq = q{};
                        }

                        # extract new seq_id
                        $seq_id = substr($line, 1);
                        $seq_id = (split /\s+/xms, $seq_id)[0];
                        ($gene) = split /t/xms, $seq_id;
                        next LINE;
                }

                # elongate current seq (seqs can be broken on several lines)
                $seq .= $line;
        }

        # add last seq to hash (if any)
        $seq_for{$gene}{$seq_id} = length $seq if $seq;

                return \%seq_for;
}


sub read_table {
        my $infile = shift;
        ### Reading input file: $infile

        my %class_for;

        open my $in, '<', $infile;
        LINE:
        while (my $line = <$in>) {
                chomp $line;

                next LINE if $line =~ m/ ^ \s* \z /xms;                 # skip empty lines
                next LINE if $line =~ m/ ^ \# /xms;                     # skip comment lines

                # Read the line and split in columns (ID + others). Split ID in gene + isoform
                my @elements = split /\t/xms, $line;
                my ($gene) = split /t/xms, $elements[0];

                my $class = $elements[4];

                # Create a hash of hash [gene -> isoforms -> values]
                my $iso = $elements[0];
                $class_for{$gene}{$iso} = $class;

        }

        return \%class_for;
}
sub read_expr {
        my $infile = shift;
        ### Reading input file: $infile

        my %expr_for;

        open my $in, '<', $infile;
        LINE:
        while (my $line = <$in>) {
                chomp $line;

                next LINE if $line =~ m/ ^ \s* \z /xms;                 # skip empty lines
                next LINE if $line =~ m/ ^ \# /xms;                     # skip comment lines
                #       next LINE if $line !~ m/ ^ EUGGR \* /xms;               # skip column headers

                # Read the line and split in ID + expression values. Split ID in isoform and gene
                my @elements = split /\t/xms, $line;
                my $iso = $elements[0];
                my ($gene) = split /t/xms, $iso;

                my $expr = sum @elements[ 1..$#elements ];
                ## $iso
                ## assert: $expr > 0

                # Create a hash of hash [gene -> isoforms -> values]
                $expr_for{$gene}{$iso} = $expr;
        }

        return \%expr_for;
}
