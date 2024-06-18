#!/usr/bin/perl
use Parallel::ForkManager;
use strict;
use warnings;
use Cwd qw(abs_path);
use Cwd qw(cwd);
use File::Spec;
use File::Basename;
use Getopt::Long;

my %inmod =  ( 
	               'C+m'   => 'C+m', 
    	           '5mc'   => 'C+m', 
    	           '5mC'   => 'C+m', 
  
				   'C+h'  => 'C+h',
	               '5hmc'  => 'C+h',
	               '5hmC'  => 'C+h');
my %outmod =  ( 
	               "C+h"   => '5hmC', 
    	           "C+m"   => '5mC'); 

my $version= "alpha";

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


my $commandl = "\@PG\tID:filotools\tPN:mod_dem\tVN:$version\tCL:" . $0 . " ". (join " ", @ARGV) . " ::FULL PATHS:: ";

#"Fragmentomics Interactive Library ONT tools"

##### PREPROCESSING 


my @modbams;
my @modlist = ("C+m");
my $exc_head;
my $threshold = 170;

GetOptions(
	"modbam|b=s" => \@modbams,
	"mods|m=s" => \@modlist
	#"threshold|t=s" => \$threshold,
	#"exclude_header|h=s" => $exc_head
	); 

die "Missing modbam path, use -b\n" unless @modbams; 

my @header;
my %bam;

#print(@modbams)
print(@modlist)

my $cb=0;
for my $bampath (@modbams){
	$cb++;
	$bampath = abs_path($bampath);
	
	next "duplicated modbam in input, skipping $bampath, iteration $cb \n" if exists($bam{"$bampath"});

	open(BAM, "samtools view -h $bampath |");

	my $clb=0;
	while (my $line = <BAM>){
		if ($line =~ m/^@/){ 
			if ($cb == 1){
				##### Header is inherited only from first modbam
				push @header, $line;
			}
		} else {
			$clb++;
			chomp $line;
			my @split = split("\t", $line);
			$bam{"$bampath"}{"$split[0]"}=$line;
			#unless (exists($bam{"$split[0]"})){ $bam{"$split[0]"}=$line } # questo lo possiamo anche mettere di default, non è significativamente più lento
			
			if ($clb==1){
				my $temp ; 
				my $splitpos = join(",",((startchar(@split, $temp = "MM:Z:"))[0], (startchar(@split, $temp = "ML:B:C"))[0]));
				#print $splitpos, "daje gaso $bampath \n";

				if (exists($bam{"$bampath"}{"splitpos"}) && $bam{"$bampath"}{"splitpos"} ne $splitpos){
					die "Fatal error: inconsistent bam formatting \n";
				} else {
					$bam{"$bampath"}{"splitpos"} = $splitpos ;
				}
			} elsif ($clb==200000){
				$clb=0;
			}
		}
	}
}

my @stats = @ARGV;
print "INPUTS: ", join("\n",@stats), "\n";

my $forks=12;
my $pm = Parallel::ForkManager->new($forks);
 
$pm->run_on_finish( sub {
    my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $out) = @_;
	print "done, ", join(',', ($pid, $exit_code, $exit_signal, $core_dump, $$out)), "\n"; 
});
 

for my $mod (@modlist){
	$mod = $inmod{$mod}; 
}  #= ("C+m");


foreach my $q (@stats) {
    my $pid = $pm->start and next;
	$q = abs_path($q);
	open(STATS, $q);
	

	##### DEFINE OUTPUT DIRS  (SAME PATH AS STATS INPUT FILES + /MOD) AND NAMES
	my $statsdir = dirname($q) . "/MOD/";
	if (! -e $statsdir) {
 		mkdir($statsdir) or die "Can't create $statsdir:$!\n";
	}

	my $cropname = basename($q);
	if ($cropname =~ m/\.stats$/){
		$cropname =~ s/\.stats$//
	} else {
		$cropname =~ s/\.[^.]+$//;
	}


	while (my $sline = <STATS>){
		chomp $sline;
		my @splits = split("\t", $sline);
		
		my $check_found = 0;
		for my $bampath (@modbams){

			my ($idMM, $idML) = split(",", $bam{$bampath}{"splitpos"});
			

			if (exists($bam{$bampath}{"$splits[3]"})){
				$check_found=1;

				my $bline = $bam{$bampath}{"$splits[3]"};
				#####copia da quì
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
				#print "CG STUFF: ", join(",", @CGpos), " giao:", $CGn,  "\n";
				#print "C STUFF: ", join(",", @Cpos), " giao:", $Cn,  "\n";
				my $mm = $splitb[$idMM]; $mm =~ s/MM:Z://;              # MM:Z:C+h,18,1;C+m,18,1; --  # C+h,18,1;C+m,18,1;	
				my $ml = $splitb[$idML]; $ml =~ s/ML:B:C,//;            # ML:B:C,0,1,89,7 -- # 0,1,89,7

				for my $mod (@modlist){ 
					my $Nmod = 0;
                	my $printmod;                                       # = "OUT: " . join("\t",@splitb[(0,2,3)]) . "\t$CGn\t$Nmod\n";
					
					if ($CGn > 0 && $mm) {

						my @premm = split(";", $mm);                    # (C+h,18,1) , (C+m,18,1) , ()	
						my @posmm = startchar(@premm, $mod);            # array with index of the correct element in @premm (in this case 1!)

						die "Fatal Error: Too many fields in MM which match with $mod at: $bline" if $#posmm > 0;
						
						$premm[$posmm[0]] =~ s/\Q$mod,//;                # (18,1)
						my @splitmm = split(",",  $premm[$posmm[0]]);    # (18), (1)
						my @splitml = split(",", $ml);                   # (0),(1),(89),(7)

						for my $sm (0..$#splitmm){
							##### l'elemento di splitmm fondamentalmente è l'indice del vettore con le C. solo che è incrementale, e quindi, a tutti quelli successivi al primo va aggiunto il precedente + 1 (perchè 0 vuol dire la successiva)
							if ($sm > 0){$splitmm[$sm] = $splitmm[$sm]+$splitmm[$sm-1]+1 } #increment
							my $Cposval = $Cpos[$splitmm[$sm]];
							
							die "Fatal error: less C in seq compared to MM field at: $bline" unless $Cposval;

							if (exists($CGhash{$Cposval})){
								if( $splitml[(($sm * (scalar @premm))+$posmm[0])] >= $threshold ){$Nmod++} # prendo l'indice degli slementi in MM (che sarebbe $sm), lo moltiplico per lo scalar di @premm (che indica quante possibili mod ci sono ) perchè ML ha serie di N numeri dove N è il numero di possibili mod. poi ci aggiungo l'indice della mod in MM che rappresenta di "quanto shiftare a destra" quando si sceglie l'indice di ML.
								delete($CGhash{$Cposval});
							} 

						}
					
						if (scalar(keys(%CGhash)) > 0){print "Too many CG:", scalar(keys(%CGhash)), "positioned as follows ", join("::", %CGhash), " at $bline"}

					} elsif ($CGn > 0 && !($mm)){
						$printmod = join("\t",@splitb[(0,2,3)]) . "\t$CGn\tMM_Undetected\n";
					}

					$printmod = join("\t",@splitb[(0,2,3)]) . "\t$CGn\t$Nmod\n";

					open(OUT, ">>$statsdir$cropname.$outmod{$mod}");
					print OUT $printmod;
					close(OUT);

				}
				###### copia fino a quì
			} 
		}


		if ($check_found==0){
			for my $mod (@modlist){
				# PRINTA LE BARACCHE ETC
				open(OUT, ">>$statsdir$cropname.$outmod{$mod}");
				print OUT "$splits[3]\tNA\tNA\tNA\tNA\n";
				close(OUT);
			}
		}


	}
    $pm->finish(0, \$q);
}
$pm->wait_all_children;




exit