
my $geneCoordinatesFile = $ARGV[0];
# Master list set of genes.
# For example bed/knownCanonical.genesWithPOLIIpeaks.bed
my $ucToENSTfile = $ARGV[1];
my $ensGenefile = $ARGV[2];
# For example knownToEnsembl.txt

my %uc_hash;
open(IN, $ucToENSTfile);
while(<IN>){
    chomp $_;
    my($ucID,$tid) = split(/\t/,$_);
#    print STDERR "Storing correspondence $ucID\t$tid\n";
    $uc_hash{$ucID} = $tid;
}

my %enst_hash;
open(IN, $ensGenefile);
while(<IN>){
    chomp $_;
    my ($num,$enst,$chr,$strand,$start1,$start2,$start3,$start4,$index,$exonStarts,$exonStops,$stuff,$ensg) = split(/\t/,$_);
    $enst_hash{$enst} = $ensg;
}

my $minlength = 2000;

open(IN, $geneCoordinatesFile);
while($line = <IN>) {
        chomp $line;
        @line = split("\t", $line);
        $ucID = $line[3];
#	print STDERR "Found UCID $ucID\n";
	$strand = $line[5];
	$start = $line[1];
	$stop = $line[2];
	$chr = $line[0];
	if ($strand eq "+"){
	    $genelength = $stop-$start;
	} elsif ($strand eq "-"){
	    $genelength = $stop-$start;
	}
	my $tid = $uc_hash{$ucID};
	my $gid = $enst_hash{$tid};
	#	print STDERR "Found Conversion: $tid\t$ucID\n";
	#	print STDERR "$tid\t$start\t$stop\t$strand\t$genelength\n";
	if ($gid eq "") {
	    print STDERR "ERR: No match to UCID $ucID\n";
	} else {
	    if ($genelength >= $minlength){
		print "$chr\t$start\t$stop\t$gid\t$tid\t$strand\n";
	    }
	}
       
}
close(IN);
#close(OUT);
