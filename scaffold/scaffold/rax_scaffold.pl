#!/usr/bin/perl -w
#script to run backbone constraint scaffolds in raxml
#takes a phylip morphology data file (data file below) or noct_morp.phylip
#then takes a tree file output from MrBayes posterior or any series of trees in Newixk format (in this case, samptre.tree)
##then runs raxml Markov model with backbone constraints
#rax_scaffold.pl

#require raxml in path

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
			open (OUT,">backboneConstraint.txt");
				print OUT "$text[$i]";
					close(OUT);
						@shell = `raxmlHPC-PTHREADS-SSE3 -p12345 -r backboneConstraint.txt -m MULTIGAMMA -KMK -s $data -T2 -n$i -oThyroptera_tricolor`;
							print "@shell\n";
	 	 				}