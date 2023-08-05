use File::Basename;
my $geneFile = $ARGV[0];
open(IN, $geneFile);
$prefix = basename($geneFile);
if ($prefix =~ /^([\w\.\_]+).bed$/){ $prefix = $1;}
print STDERR "$prefix\n";
open(G, ">$prefix.fullgenebody.bed");
while($line = <IN>) {
	chomp $line;
	@line = split("\t", $line);
	$tid = $line[4];
	$gid = $line[3];
	$chr = $line[0];
	$start = $line[1];
	$end = $line[2];
	$strand = $line[5];
#	print STDERR "Parsed $chr\t$start\t$end\tTID:$tid\tStrand:$strand\n";
	if ($strand eq '+') {
                $gstart = $start + 300;
                $gend = $end;
	} else {
                $gstart = $start;
                $gend = $end - 300;
	}
	print G "$chr\t$gstart\t$gend\t$gid\t$tid\t$strand\n";
}
close(IN);
close(G);
