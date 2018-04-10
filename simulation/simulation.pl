#!/usr/bin/perl

#run: perl simulation.pl windows_5cpgs <depth>

## haplotype frequencies
$f[0] = 0.4;
$f[1] = 0.2;
$f[2] = 0.3;
$f[3] = 0.1;

#$N = 300; # number of reads to sample
$Read_len = 150;
$Frag_len = 500;
$Win_len = 1000;

$N = int(($ARGV[1]*$Win_len/$Read_len) + 0.5);  # calculate the number of reads to be taken by using the depth given in the command line




open(IN,$ARGV[0]);
while(<IN>) {
@fields = split(/\s+/,$_);
$cpg_id = @fields[1];
$pos_cpgs_hap = @fields[2];
@cpg_positions = @fields[3..$#fields];

    if($#cpg_positions >= 2) { ## must have at least two CpGs to get 4 distinct haplotypes; currently hard-coded - must be changed to change the number of haplotypes
    &get_reads(\@cpg_positions);
print STDERR $#cpg_positions+1, " cpgs\n";
        system("python call.py $outfile 0 $pos_cpgs_hap sim_$."); # $. is the number of the line (=sim)
        system("python sim_entropy_methfreq.py $outfile2 sim_$.");
    }
}
close(IN);




## A function that returns a set of reads (note: currently this script is using global variables)

sub get_reads {
    $outfile = "sim_$._Simulated_Data";
    open(OUT,">$outfile");
    $outfile2 = "sim_$._Simulated_Hap_Freqs";
    open(OUT2,">$outfile2");
@cpg_positions = @{$_[0]}; # the argument to the function (which is an array)
#print join(" ",@cpg_positions),"\n";
@hap = ();
@hapstring = ();

    ## determine the haplotype patterns (based on the input number of cpg sites)
    for($h = 0; $h <= $#f; $h++) {
        for($p = 0; $p <= $#cpg_positions; $p++) {
            if(rand(1) < 0.5) {
            $hap[$h][$p] = 0;
            $hapstring[$h] = $hapstring[$h] . "0";
            }
            else {
            $hap[$h][$p] = 1;
            $hapstring[$h] = $hapstring[$h] . "1";
            }
        }
        FOR: for($k = 0; $k < $h; $k++) {
            if($hapstring[$h] eq $hapstring[$k]) {
            delete $hapstring[$h]; # <--------
            $h--; ## means we've produced this haplotype before; try again
            print STDERR "found a repeated haplotype - retrying...\n";
            last FOR;
            }
        }
    }

    ## sample reads
    for($i = 0; $i < $N; $i++) {
        #        print STDERR "$i\n";

    $read = "";
    $p = int(rand($Win_len + $Frag_len)) - $Frag_len; # position of the read - note the read overlaps the window if it starts anywhere on [-Frag_len, Win_len]
    $random_hap = rand(1); # used next to pick a haplotype at random according to the haplotype frequencies
    $sum = 0;
        FOR1: for($k = 0; $k <= $#f; $k++) {
            if($random_hap < $f[$k] + $sum) {
                print OUT2  "$hapstring[$k]	$f[$k]	$ARGV[0]\n";
            $rhap = $k;
#print STDERR "$k $f[$k] $rhap\n";
            last FOR1;
            }
        $sum = $sum + $f[$k];
        }
        for($pos = 0; $pos <= $#cpg_positions; $pos++) {
            if( ($cpg_positions[$pos] >= $p && $cpg_positions[$pos] < $p + $Read_len) || ($cpg_positions[$pos] >= $p + $Frag_len - $Read_len && $cpg_positions[$pos] < $p + $Frag_len)) {
                print OUT "$pos	$hap[$rhap][$pos]	$i	$ARGV[0]\n";
#                print "$pos $hap[$rhap][$pos] $i\n";
            }
        }
    }
}








