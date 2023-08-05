open(IN, $ARGV[1]);
$linecount = 0;
while($line = <IN>) {
	$linecount++;
	if ($linecount == 1) {
		next;
	}
	chomp $line;
	@line = split("\t", $line);
	$tid = $line[4];
	$start{$tid} = $line[1];
	$end{$tid} = $line[2];
}
close(IN);
open(IN, $ARGV[0]);
while($line = <IN>) {
	chomp $line;
	@line = split("\t", $line);
	$tid = $line[4];
	$geneID = $line[3];
	$chr = $line[0];
	$chr_hash{$chr} = 1;
	$start = $start{$tid} - 500;
	$stop = $end{$tid} + 500;
	push @{$gene_starts{$chr}}, $start;
	$gene_stop{$chr}{$start} = $stop;
	$geneID{$chr}{$start} = $geneID;
}
close(IN);
foreach $i (keys %chr_hash) {
	@{$sorted_gene_starts{$i}} = sort {$a <=> $b} @{$gene_starts{$i}};
} 	
foreach $i (keys %chr_hash) {
	foreach $j (0 .. $#{$sorted_gene_starts{$i}}) {
		$switch = 0;
		$k = 0;
		if ($j == $#{$sorted_gene_starts{$i}}) {
			$switch = 1;
		}
		while ($switch == 0) {
			$k++;
			$l = $j + $k;
			if (${$sorted_gene_starts{$i}}[$l] < $gene_stop{$i}{${$sorted_gene_starts{$i}}[$j]}) {
#				print "$i	$j	$k	$l	${$sorted_gene_starts{$i}}[$l]	$gene_stop{$i}{${$sorted_gene_starts{$i}}[$j]}\n";
				$bad_gene{$geneID{$i}{${$sorted_gene_starts{$i}}[$j]}} = 1;
				$bad_gene{$geneID{$i}{${$sorted_gene_starts{$i}}[$l]}} = 1;
				if ($l == $#{$sorted_gene_starts{$i}}) {
					$switch = 1;
				}
			} else {
				$switch = 1;
			}
		}
	}
}
#open(OUT, ">final_transcripts.bed");
open(IN, $ARGV[0]);
while($line = <IN>) {
	chomp $line;
	@line = split("\t", $line);
	$geneID = $line[3];
	if ($bad_gene{$geneID} == 1) {
		next;
	} else {
		print STDOUT "$line\n";
	}
}
