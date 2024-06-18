#!/usr/bin/perl
use Parallel::ForkManager;
use strict;
use warnings;
use Getopt::Long;


my %hash;

sub startchar(\@\$){
	##riportami GLI index degli elementi che cominciano come voglio io.
	my ($ref_one, $ref_two) = @_;
	my @opts = @{$ref_one};
	my $target = ${$ref_two};

	my $indx = "NA";
	my @indx;

	for (0 .. $#opts) { 
		my $tst = $opts[$_];		
		if ($tst =~ m/^\Q$target/){
			$indx = $_ ;
			push @indx, $_ ;
		}
	}
	return @indx;
}

#my $revcomp = reverse $origin_seq;
#$revcomp =~ tr/ATGCatgc/TACGtacg/;
#
#print $revcomp

#my $strand;
#if ($flag & 0x10){
#    $strand  =  "-";
#} else {
#    $strand  = "+";
#};

# QUESTO è OK!
#if ($flag & 0x10){
#   $seq  =  reverse $seq;
#	$seq =~ tr/ATGCatgc/TACGtacg/;
#} 

##### anche questo buono già traslato
#####my $seq = "ACAGCACCACCAAAGGCAGCACAATGAGAACAGCATTCTCCTCAACAGACAAGCTGGGAGTATCTAGACACCCGACCTCAATAGCTCCAGAACAGCCCTAAAACATTTCCTCCCTAACCACCACTCAAGTCACCAGCTTGGAAAGTATTAAGAAAACCCAAATC";
#####my $mm = "MM:Z:C+h,22;C+m,22;";
#####my $ml = "ML:B:C,0,255";
#####
#####
#####sub startchar(\@\$){
#####	##riportami l'index dell'elemento dell'array fornito, che comincia come voglio io.
#####	my ($ref_one, $ref_two) = @_;
#####	my @opts = @{$ref_one};
#####	my $target = ${$ref_two};
#####
#####	my $indx = "NA";
#####	my @indx;
#####
#####	for (0 .. $#opts) { 
#####		my $tst = $opts[$_];
#####		#print $tst, $_ ,  "\n";
#####		
#####		if ($tst =~ m/^\Q$target/){
#####			$indx = $_ ;
#####			push @indx, $_ ;
#####			#print "dioscannat \n";
#####		}
#####
#####		#last unless $indx eq "NA";
#####	}
#####	#return $indx;
#####	return @indx;
#####	#print $indx, " hurray! \n";
#####}
#####
#####my @sandro = ("C+h,22", "C+m,22", "asdf", "asdflksdfk");
#####my $tar = "C+h";
#####
######print startchar(@sandro, my $g = "gibbo"), startchar(@sandro, my $g = "gigio");
#####my $g ; 
#####print join(",",startchar(@sandro, $g = "C+h")  )  , "\n";
#####print join(",",startchar(@sandro, $g = "asdf")  ) , "\n";
#####print join(",",startchar(@sandro, $g = "C+m")    ), "\n";
#####print join(",",startchar(@sandro, $g = "C+")    ) , "\n";
#####
#####my $drico = (startchar(@sandro, $g = "asdf"))[0];
#####print $drico, "\n";

my %inmod =  ( 
	               'C+m'   => 'C+m', 
    	           '5mc'   => 'C+m', 
    	           '5mC'   => 'C+m', 
  
				   'C+h '  => 'C+h',
	               '5hmc'  => 'C+h',
	               '5hmC'  => 'C+h');

my @modlist = ("C+m");

for my $mod (@modlist){
	$mod = $inmod{$mod}; 
}  #= ("C+m");

my $idMM = 12;
my $idML = 13;
#my $printout = "$splits[3]\tNA\tNA\tNA\n";
my $threshold = 170;
my $gianc = 0;
while (my $bline = <>){
$gianc++;
print "line: $gianc \n";
				my @splitb = split("\t", $bline);
				my $flag = $splitb[1];
				my $seq = $splitb[9];

				if ($flag & 0x10){
				    $seq  =  reverse $seq;
					$seq =~ tr/ATGCatgc/TACGtacg/;
				} 

				##### find context positions
				my $context_offset = -1;
				my @CGpos;
				my %CGhash;
				my $counter_CG = -1;
				while($seq =~ m/CG/g){
					$counter_CG++;
                    push @CGpos, pos($seq)+$context_offset; # 1 based position in read
                    #push @CGpos, pos($seq)+$context_offset-1; # così è la zero based position, ma a noi non interessa?
					$CGhash{(pos($seq)+$context_offset)} = $counter_CG; # 1 based position
                }
				my $CGn = scalar @CGpos;

				##### find base positions
				my @Cpos;
				while($seq =~ m/C/g){
                    push @Cpos, pos($seq); # 1 based position in read
                    #push @Cpos, pos($seq)-1; # così è la zero based position, ma a noi non interessa?
                }
				my $Cn = scalar @Cpos;

				print "CG STUFF: ", join(",", @CGpos), " giao:", $CGn,  "\n";
				print "C STUFF: ", join(",", @Cpos), " giao:", $Cn,  "\n";

				my $mm = $splitb[$idMM]; $mm =~ s/MM:Z://;              # MM:Z:C+h,18,1;C+m,18,1; --  # C+h,18,1;C+m,18,1;	
				my $ml = $splitb[$idML]; $ml =~ s/ML:B:C,//;            # ML:B:C,0,1,89,7 -- # 0,1,89,7

				for my $mod (@modlist){ 
					my $Nmod = 0;
                	my $printmod;                                       # = "OUT: " . join("\t",@splitb[(0,2,3)]) . "\t$CGn\t$Nmod\n";
					
					if ($CGn > 0 && $mm) {
						#print "bellafregna " , join(" , ", $flag, $seq, $mm, $ml);

						my @premm = split(";", $mm);                    # (C+h,18,1) , (C+m,18,1) , ()	
						my @posmm = startchar(@premm, $mod);            # array with index of the correct element in @premm (in this case 1!)

						die "Fatal Error: Too many fields in MM which match with $mod at: $bline" if $#posmm > 0;
						
						$premm[$posmm[0]] =~ s/\Q$mod,//;                # (18,1)
						my @splitmm = split(",",  $premm[$posmm[0]]);    # (18), (1)
						my @splitml = split(",", $ml);                   # (0),(1),(89),(7)

						#print join(" mm ", @splitmm), " sep ",join(" ml " , @splitml), "\n";
						#print "rangeeee:  ", join("::",(0..$#splitmm)), " \n";
						for my $sm (0..$#splitmm){
							#print "rangooooo $sm \n";
							##### l'elemento di splitmm fondamentalmente è l'indice del vettore con le C. solo che è incrementale, e quindi, a tutti quelli successivi al primo va aggiunto il precedente + 1 (perchè 0 vuol dire la successiva)
							if ($sm > 0){$splitmm[$sm] = $splitmm[$sm]+$splitmm[$sm-1]+1 } #increment
							my $Cposval = $Cpos[$splitmm[$sm]];
							
							die "Fatal error: less C in seq compared to MM field at: $bline" unless $Cposval;

							if (exists($CGhash{$Cposval})){
								if( $splitml[(($sm * (scalar @premm))+$posmm[0])] >= $threshold ){$Nmod++} # prendo l'indice degli slementi in MM (che sarebbe $sm), lo moltiplico per lo scalar di @premm (che indica quante possibili mod ci sono ) perchè ML ha serie di N numeri dove N è il numero di possibili mod. poi ci aggiungo l'indice della mod in MM che rappresenta di "quanto shiftare a destra" quando si sceglie l'indice di ML.
								delete($CGhash{$Cposval});
							} 

						}
					
						#print "scalar CGHASH: ",scalar(keys(%CGhash)), " \n";
						if (scalar(keys(%CGhash)) > 0){print "accidenti troppe CG:", scalar(keys(%CGhash)), "positioned as follows ", join("::", %CGhash), " at $bline"}

					} elsif ($CGn > 0 && !($mm)){
						$printmod = join("\t",@splitb[(0,2,3)]) . "\t$CGn\tMM_Undetected\n";
					}

					$printmod = join("\t",@splitb[(0,2,3)]) . "\t$CGn\t$Nmod\n";

					print $printmod;

				}
}

				###### siamo rimasti quì, è da andare avanti con questa questione volevo fare if ($mm) {} else {aggiungi al printing l'output vuoto (tipo 0?)}
				###### cerchiamo l'indice della mod, e nell'ml il valore di prob sarà uguale al numero di elementi/2+indice della mod? no ma che cazzo dico, sarà i numeri pari compresi tra 0 e il numero di elementi-1+indicemod.
#					MM:Z:   ML:B:C
#   MM:Z:C+h,14,0,2;C+m,14,0,2;     ML:B:C,6,44,55,125,193,155

exit

