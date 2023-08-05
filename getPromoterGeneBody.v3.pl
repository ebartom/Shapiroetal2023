use File::Basename;
my $geneFile = $ARGV[0];
open(IN, $geneFile);
$prefix = basename($geneFile);
if ($prefix =~ /^([\w\.\_]+).bed$/){ $prefix = $1;}
print STDERR "$prefix\n";
open(P, ">$prefix.promoter.bed");
open(G, ">$prefix.genebody.bed");
while($line = <IN>) {
	chomp $line;
	@line = split("\t", $line);
	$tid = $line[4];
	$gid = $line[3];
	$chr = $line[0];
	$start = $line[1];
	$end = $line[2];
	$strand = $line[5];
	if ($strand eq '+') {
		$pstart = $start - 100;
		$pend = $start + 300;
        $gstart = $start + 300;
        $gend = $start + 2000;
	} else {
	    $pstart = $end - 300;
	    $pend = $end + 100;
	    $gstart = $end - 2000;
	    $gend = $end - 300;
	}
	print P "$chr\t$pstart\t$pend\t$gid\t$tid\t$strand\n";
	print G "$chr\t$gstart\t$gend\t$gid\t$tid\t$strand\n";
}
close(IN);
close(P);
close(G);
