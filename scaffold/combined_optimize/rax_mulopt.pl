#!/usr/bin/perl -w
#script to run combined data optimizations on scaffold constrained trees in raxml
#takes a phylip combined data file (data file below) or noct_combi.phylip
#then takes a tree file output from raxml MK output or any series of trees in Newick format (in this case, scaf.bestTree)
##then runs raxml combined data model on the tree using a partition file (parts.txt)


#rax_mulopt.pl

print "Please enter data file: ";
$data=<STDIN>;
chomp($data);

print "Please enter tree file: ";
$treefile=<STDIN>;
chomp($treefile);
	open (IN,"<$treefile") or die "Cannot open $treefile\n";
		my @text=();
			while($line=<IN>)	{
				chomp($line);
					push(@text,$line);
	 				}
	 					close (IN);
	 				
	 	my $N = scalar @text;
	 	my $count=0;
		for	(my $i = 0; $i < $N; $i++)	{
			open (OUT,">temptree.txt");
				print OUT "$text[$i]";
					close(OUT);
						@shell = `raxmlHPC-PTHREADS-SSE3 -f e -t temptree.txt -m MULTIGAMMA -KMK -s $data -q parts.txt -T2 -n$data$i -oThyroptera_tricolor`;
							print "@shell\n";
	 	 				}