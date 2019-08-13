#!/usr/bin/env perl

use autodie;
use Smart::Comments;
use Cwd qw(cwd);
my $dir = cwd;


# get new id old id equivalence
my $new_old_ids_tab = shift;

if ($new_old_ids_tab) {
        $new_old_ids_tab = $new_old_ids_tab;

} else {

die <<"EOT";
Usage: $0 <new old ids tab> <blast file in outfmt6>
EOT
}

# get outfmt6 files from blastn
my @files = @ARGV;

if (@files) {
        @files = @files;

} else {

# open directory and get file's names
opendir(DIR, $dir) or die "Could not open '$dir' for reading: $!\n";
@files = grep(/\.50-90.blastn$/,readdir(DIR));
closedir(DIR);
}



# function to parse outfmt6
my %bad_ids_for = read_outfmt6(@files);
my %new_old_ids_for = read_ids($new_old_ids_tab)

### %bad_ids_for
### %new_old_ids_for

my $outfile = 'outfile.txt';
open my $out, '>' , $outfile;

foreach $id (keys %bad_ids_for){
my $conta        = $new_old_ids_for{$id};
my $conta_fun    = $bad_ids_for{$id}{fun}
my $conta_acc    = $bad_ids_for{$id}{acc};
say {$out} join qq{\t}, $conta, $id, $conta_acc, $conta_;
}
close $out;

sub read_ids{
my $file = @_;
my %hash;

        open my $in, '<', $file;
        LINE:
                while (my $line = <$in>) {
                chomp $line;

                next LINE if $line =~ m/ ^ \s* \z /xms;                 # skip empty lines
                next LINE if $line =~ m/ ^ \# /xms;                     # skip comment lines

                # Read the line and split in columns (old ID + new ID)
                my @elements = split "\t", $line;
                my $id = $elements[0];
                my $new_id = $elements[1];

                # Create a hash of hash [old ID -> new ID]
                $hash{$id} = $new_id;
        }
        return %hash;
}

sub read_outfmt6{
my @files = @_;
my %hash;

        foreach $file (@files) {
                ### Reading input file: $file

                # customized outfmt6 from blastn analysis
                # default is => qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
                # cutomized -outfmt "6 qseqid sscinames sskingdoms stitle sseqid staxids sacc"

                open my $in, '<', $file;
                LINE:
                        while (my $line = <$in>) {
                        chomp $line;

                        next LINE if $line =~ m/ ^ \s* \z /xms;                 # skip empty lines
                        next LINE if $line =~ m/ ^ \# /xms;                     # skip comment lines

                        # Read the line and split in columns (ID + sp + acc)
                        my @elements = split "\t", $line;

                        my $id       = $elements[0];
                        my $sp       = $elements[1];
                        next LINE if $sp =~ m/Euglena/xms; # Only lines that are not Euglena!

                        my $function = $elements[4];
                        my $acc      = $elements[6];

                        # Create a hash of hash [ID -> values]
                        $hash{$id}{sp}  = $sp;
                        $hash{$id}{fun} = $function;
                        $hash{$id}{acc} = $acc;
                }
        }
        return %hash;
}



