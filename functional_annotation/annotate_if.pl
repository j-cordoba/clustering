#!/usr/bin/env perl

use Modern::Perl '2011';
use autodie;

use Tie::IxHash;
use Smart::Comments;

use List::AllUtils qw(sum reduce max);
use File::Basename;
use Path::Class 'file';


unless (@ARGV == 4) {
die <<"EOT";
Usage: $0 <fasta> <ann1.tab> <ann2.tab> <ann3.tab>
EOT
}

my $fasta = shift;
my $ann1 = shift;
my $ann2 = shift;
my $ann3 = shift;

# .... my ann3 = shift;

my $seq_for  = read_fasta($fasta);
my $ann1_for = read_annotation($ann1);
my $ann2_for = read_annotation($ann2);
my $ann3_for = read_annotation($ann3);

# Obtain gene identifiers from fasta file
my @ids = keys %$seq_for;
### @ids

## Prepare output file to be written

my $outfile = 'annotation.output';
open my $out, '>' , $outfile;


for my $seq_id (@ids) {
	$seq_for->{$seq_id}{GO}   = $ann1_for->{$seq_id}{GO};
	$seq_for->{$seq_id}{KEGG} = $ann1_for->{$seq_id}{KEGG};
	$seq_for->{$seq_id}{GO}   = $ann2_for->{$seq_id}{GO}   unless $seq_for->{$seq_id}{GO};
	$seq_for->{$seq_id}{KEGG} = $ann2_for->{$seq_id}{KEGG} unless $seq_for->{$seq_id}{KEGG};
        $seq_for->{$seq_id}{GO}   = $ann3_for->{$seq_id}{GO}   unless $seq_for->{$seq_id}{GO};
        $seq_for->{$seq_id}{KEGG} = $ann3_for->{$seq_id}{KEGG} unless $seq_for->{$seq_id}{KEGG};

	
	# Writing results to a table
	my $GO   = $seq_for->{$seq_id}{GO};
	my $KEGG = $seq_for->{$seq_id}{KEGG};
	
	say {$out} join qq{\t}, $seq_id, $GO // '' , $KEGG // '';
}

close $out;


sub read_fasta {
        my $infile = shift;
        ### Reading input file: $infile

        open my $in, '<', $infile;

        my $seq_id;
        my $seq;
        tie my %seq_for, 'Tie::IxHash';         # preserve original seq order

        LINE:
        while (my $line = <$in>) {
                chomp $line;
                
       			next LINE if $line =~ m/ ^ \s* \z /xms; 		# skip empty lines
			next LINE if $line =~ m/ ^ \# /xms;			# skip comment lines


                # at each '>' char...
                if (substr($line, 0, 1) eq '>') {

			# add current seq to hash (if any)
			if ($seq) {
				# Create a hash of seq
				$seq_for{$seq_id} = undef;
				$seq = q{};
			}

			# extract new seq_id
			$seq_id = substr($line, 1);
			$seq_id = (split /\s+/xms, $seq_id)[0];
			next LINE;
                }

                # elongate current seq (seqs can be broken on several lines)
                $seq .= $line;
        }

        # add last seq to hash (if any)
        $seq_for{$seq_id} = undef;
		
		return \%seq_for;
}


sub read_annotation{ 
	my $infile = shift;
        ### Reading input file: $infile

	my %ann_for;	

	open my $in, '<', $infile;
        
	LINE:
        while (my $line = <$in>) {
                chomp $line;
                
		next LINE if $line =~ m/ ^ \s* \z /xms; 		# skip empty lines
		next LINE if $line =~ m/ ^ \# /xms;			# skip comment lines
		
	        # Read the line and split in columns (ID + GO + KEGG)
		my @elements = split /\t/xms, $line;
		my $id = $elements[0];
		my $GO = $elements[1];
		my $KEGG = $elements[2];
		# Create a hash of hash [gene -> isoforms -> values]
		$ann_for{$id}{GO} = $GO;
		$ann_for{$id}{KEGG} = $KEGG;
	}

	return \%ann_for;
}

