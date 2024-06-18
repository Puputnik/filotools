#!/usr/bin/perl
use Parallel::ForkManager;
use strict;
use warnings;
use Cwd qw(abs_path);
use Cwd qw(cwd);
use File::Spec;
use File::Basename;
use Getopt::Long;
use Data::Dumper;

my %inmod =  ( 
	               'C+m'   => 'C+m?', 
	               'C+m?'   => 'C+m?', 
    	           '5mc'   => 'C+m?', 
    	           '5mC'   => 'C+m?', 
  
				   'C+h'  => 'C+h?',
				   'C+h?'  => 'C+h?',
	               '5hmc'  => 'C+h?',
	               '5hmC'  => 'C+h?');
my %outmod =  ( 
	               "C+h?"   => '5hmC', 
    	           "C+m?"   => '5mC'); 

my $version= "dorado";

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

sub get_config(\$){
	my ($ref_one) = @_;
	
	my $cpath = ${$ref_one};
	my %conf;

	open(CONF, $cpath);
	for my $cline (<CONF>){
		$cline =~ s/ +/ /;
		my @csplit = split(" ", $cline);
		$conf{"$csplit[0]"} = $csplit[1];
	}
	return %conf;
}

sub set_defaults(\%){
	my ($ref_one) = @_;
	my %conf = %{$ref_one};

    my @fields = qw(stats_path       
			        motif_path       
			        motif_count_path 
			        readl_path       
			        raw_readl_path   
			        merged_fastq_path);

    my %hashd =     qw(stats_path            bam_path/STATS 
                 	   motif_path          stats_path/MOTIF
                 	   motif_count_path    stats_path/MOTIF_COUNTS
                 	   readl_path          stats_path/READLENGTH_COUNTS
                 	   raw_readl_path      fastq_path/READL_RAW
                 	   merged_fastq_path   fastq_path/MERGED);

	for my $f (@fields){
		my $check = 0;

		if(!(exists $conf{"$f"} )){
			$check = 1;
		}  elsif ($conf{"$f"} eq "default"){
			$check = 1;
		}
		if ($check == 1){
			my @splitf = split("/",$hashd{$f});
			$conf{$f} = $conf{$splitf[0]} . "/" . $splitf[1];
		}

	}
	return %conf	
}

my %config = get_config(my $config_path=$ENV{config_path});

print Dumper(\%config);

%config = set_defaults(%config);

print Dumper(\%config);

my $commandl = "\@PG\tID:filotools\tPN:mod_dem\tVN:$version\tCL:" . $0 . " ". (join " ", @ARGV) . " ::CONFIG PATH:: " . abs_path($config_path);

#"Fragmentomics Interactive Library ONT tools"

##### PREPROCESSING 


my @modbams;
#my @modlist = ("C+m");
my @modlist ;
my $exc_head;
my $threshold  = 170;
my $force_config = 0;

my $mod_out="default";
my $modbam_out="default";
my @outputs_set;

my $auto=0;
my $help_var=0;

GetOptions(
	#### general 
	"mods|m=s"        => \@modlist,
	"threshold|t=s"   => \$threshold,
	"ouputs|w=s"      => \@outputs_set,
	"auto|a"          => \$auto,

	#### if not config
		#### inputs
	"modbam|b=s"      => \@modbams,
		#### outputs
	"modbam_out|o=s"  => \$modbam_out,
	"mod_out|u=s"     => \$mod_out,
	
	#### if config
	"force_config|C"  => \$force_config, #### force inputs modbam, and outputs modbam_out mod_out
	"help|h" => \$help_var
	#"exclude_header|h=s" => $exc_head
	); 


my $usage = "
			USAGE:

			standard mode (all stats files must belong to the same run and modbam file)
			filotools mod_dem -m 5mC -m 5hmC -b <path/to/modbam.bam> -o <path/to/store/demultiplexed_output_modbam/> -u <path/to/store/per_read_output/>  *.stats

			auto mode (can be launched on stats files from multiple runs, as soon as everything is named in Filo's format)
			filotools mod_dem -m 5mC -m 5hmC -a -o <path/to/store/demultiplexed_output_modbam/> -u <path/to/store/per_read_output/>  *.stats
			auto mode with config
			filotools mod_dem -m 5mC -m 5hmC -a -C  *.stats



			can work with 3 kind of inputs (inferred from input extension)
		    1) .sam/.bam
			2) .stats
			3) else (treated as space/tab separated text file, first column should contain read IDs)
			
			generates 2 kind of ouputs (by default )
			1) demultiplexed modbams                    (default output folder : in_folder/MODBAM)
			2) demultiplexed per-read-modfication files (default output folder : in_folder/MOD)
			
			general flags:
			--ouputs     | -w : can be modbam or per_read, can be specified multiple times re-entering -w (default: modbam AND per_read)
			--mods       | -m : modifications to analyze (mod available: C+m/5mc/5mC and C+h/5hmc/5hmC/) (default: C+m)
			--threshold  | -t : likelihood threshold for mod filtering, meaningful only for per_read output (default: 170)
			--auto       | -a : automatically detect library name, use only if filenames are in Filo's format and basecall_path is defined in config file. 
			
			input paths:
			--modbam     | -b : path to modbam  
			
			output paths:
			--modbam_out | -o : path to store demultiplexed modbam files
			--mod_out    | -u : path to store per_read modification files
			
			if --auto is set, CONFIG is used only to infer modbam input, basecall_path, fastq_path and merged_fastq_path should be defined in config file
			if --force_config is set, CONFIG is used for output paths (--modbam_out, --mod_out) (note: if --auto is not set, you still need to provide an input path with --modbam flag)


			";

if ($help_var == 1){
		print $usage;
		exit ;
}

if (scalar @modlist == 0){
	@modlist = ("C+m");
}

my %outset;
if (scalar @outputs_set == 0){
	@outputs_set=qw(modbam per_read);
}
for my $o (@outputs_set){
	if ($o ne "modbam" && $o ne "per_read"){
		die "unrecognized output $o; quitting \n";
	} else {
		$outset{$o}=1;
	}
}

#print($force_config);
#print($threshold);


my @header;
my %bam;

#print(@modbams)
print join(" ", @modlist) . "\n";
my @stats = @ARGV;


if ($auto == 1){
	die "modbams provided while --auto is set, quitting\n" if @modbams; 
	unless(exists $config{"basecall_path"} && exists $config{"fastq_path"} && exists $config{"merged_fastq_path"}){
		die "--auto is set, missing config arguments, check basecall_path, fastq_path, merged_fastq_path\n";
	}
	OUTER: foreach my $q (@stats){
		my $cropped = basename($q);
		$cropped =~ s/\.[^.]+$//;
		my @splitcrop;
		#### if exists path fastq 
		if ( -e "$config{fastq_path}/$cropped.fastq" || -e "$config{fastq_path}/$cropped.fastq.gz" ){
			@splitcrop = split("-",$cropped);
			my $libraryname = $splitcrop[4] . "_" . $splitcrop[5] . "_count" . $splitcrop[6];
			push(@modbams, "$config{basecall_path}/$libraryname/*/*/DORAMOD/modcalls.bam")
		} else {
			#### if exists path merged fastq 
			if ( -e "$config{merged_fastq_path}/logs_merge/$cropped.log"){
				
				open(LOG,"$config{merged_fastq_path}/logs_merge/$cropped.log");

				for my $log (<LOG>){
					chomp $log;
					my @splitlog = split("\t", $log);					
					
					unless ($splitlog[4] eq "Complete"){
						next OUTER;
					}
					my @splitfields = split(",", $splitlog[3]);

					for my $field (@splitfields){
						@splitcrop = split("-", $field);
						my $libraryname = $splitcrop[4] . "_" . $splitcrop[5] . "_count" . $splitcrop[6];
						push(@modbams, "$config{basecall_path}/$libraryname/*/*/DORAMOD/modcalls.bam")
					}
				}
			}
		}
	}
} 

print join("\n",@modbams) . "\n";

die "Missing modbam path, use -b\n" unless @modbams; 

my $cb=0;

$commandl = $commandl . " ::MODBAM PATHS:: ";


##### riempi l'hash del bam (%bam) con le linee dei bam.
FORBAM: for my $bampath (@modbams){
	$cb++;
	$bampath = abs_path($bampath);
	

	if (exists $bam{"$bampath"}){
		print "duplicated modbam in input, skipping $bampath, iteration $cb \n" ;
		next FORBAM;
	} else {
		$commandl = $commandl . " " . $bampath;
	}
	
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
				#### this is just a control that happens every 200k reads, to check if the bam is consistent in the position of the MM and ML columns
				#### $splitpos contiene l'index della colonna i MM e ML

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

print "keys: ", join("\n",keys %bam), "\n";

$commandl = $commandl . " ::INPUT PATHS:: " . join(" ",@stats);
#for my $q (@stats){
#	$commandl = $commandl . " " . abs_path($q);
#}

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

	my $format;
	
	##### evaluate input format
	if ($q =~ /.stats$/){
		$format = "stats";
		open(STATS, $q);

	} elsif ($q =~ /.bam$/ || $q =~ /.sam$/ ){
		$format = "sambam";
		open(STATS, "samtools view $q |");

	} else {
		$format = "ids";
		open(STATS, $q);
	}


	##### DEFINE OUTPUT DIRS  (SAME PATH AS STATS INPUT FILES + /MOD) AND NAMES
	##### outputs for modbam
	my $modoutdir;
	if ($mod_out ne "default"){
		$modoutdir = $mod_out;
	} else {
		$modoutdir = dirname($q) . "/MOD/";
		if ($force_config == 1 && exists $config{"mod_path"}){
			if ($config{"mod_path"} ne "default" && $config{"mod_path"} ne ""){
				$modoutdir=$config{"mod_path"}
			}
		}
	}

	my $modbamoutdir;
	if ($modbam_out ne "default"){
		$modbamoutdir = $modbam_out;
	} else {
		$modbamoutdir = dirname($q) . "/MODBAM/";
		if ($force_config == 1 && exists $config{"modbam_path"}){
			if ($config{"modbam_path"} ne "default" && $config{"modbam_path"} ne ""){
				$modbamoutdir=$config{"modbam_path"}
			}
		}
	}


	#if ($mod_out eq "default" && $force_config == 0){
	#	$modoutdir = dirname($q) . "/MOD/";
	#} elsif ($force_config == 1){
	#	$modoutdir=$config{"mod_path"};
	#} else {
	#}

	#my $modbamoutdir
	#if ($modbam_out eq "default"){
	#	$modbamoutdir = dirname($q) . "/MODBAM/";
	#} elsif ($force_config == 1){
	#	$modbamoutdir=$config{"modbam_path"};
	#} else {
	#	$modbamoutdir = $modbam_out;
	#}


	my $cropname = basename($q);
	if ($cropname =~ m/\.stats$/){
		$cropname =~ s/\.stats$//
	} else {
		$cropname =~ s/\.[^.]+$//;
	}


##### CREATE DIRECTORIES AND EMPTY OUTPUT FILES
	if (exists $outset{"per_read"}){
		if (! -e $modoutdir) {
 			mkdir($modoutdir) or die "Can't create $modoutdir:$!\n";
		}
		for my $mod (@modlist){ 
			#open(OUT, ">$modoutdir/$cropname.bam");
			open(OUT, ">$modoutdir/$cropname.$outmod{$mod}");
			print OUT "#" . "$commandl\n";
			close(OUT);
		}
	}

	if (exists $outset{"modbam"}){
		if (! -e $modbamoutdir) {
 			mkdir($modbamoutdir) or die "Can't create $modbamoutdir:$!\n";
		}
		#### print header
		open(OUT, ">$modbamoutdir/$cropname\_mod_mappings.sam");
		print OUT @header;
		print OUT "$commandl\n";
		close(OUT);
	}

##### PARSE INPUTS
	while (my $sline = <STATS>){
		chomp $sline;
		my @splits;
		my $readid;
		##### CHECK FORMAT
		if ($format eq "stats"){
			@splits = split("\t", $sline);
			$readid = $splits[3];
		} elsif ($format eq "sambam"){
			@splits = split("\t", $sline);
			$readid = $splits[0];
		} else {
			@splits = split(' ', $sline);
			$readid = $splits[0];
		}

		my $check_found = 0;
#		for my $bampath (@modbams){
		for my $bampath (keys %bam){

			my ($idMM, $idML) = split(",", $bam{$bampath}{"splitpos"});
			

			if (exists($bam{$bampath}{"$readid"})){
				$check_found=1;
				my $bline = $bam{$bampath}{"$readid"};
				
				if (exists $outset{"modbam"}){
					open(OUT, ">>$modbamoutdir/$cropname\_mod_mappings.sam");
					print OUT "$bline\n";
					close(OUT);
				}

				#####copia da quì
				my @splitb = split("\t", $bline);
				my $flag = $splitb[1];
				my $seq = $splitb[9];

				my $rightclip = 0 ;
				my $leftclip = 0 ;
				my $len = 0 ;

				unless ( $splitb[2] eq '*'){
					my @cig = ( $splitb[5] =~ /(\d+\w)/g);
        			for my $c (@cig){
        			    if ($c =~ /(\d+)[MX=DN]/g){
        			        $len = $len + $1;
        			    }
        			}

        			my $cigo = $splitb[5];

        			my $H1 = 0;
        			my $H2 = 0;

        			my $S1 = 0;
        			my $S2 = 0;


        			#### calculating length of leftmost hardclip
        			if ($cigo =~ m/^(\d+)H/){
        			    $H1 = $1;
        			    $cigo =~ s/^$H1\H//;
        			} 

        			#### calculating length of rightmost hardclip
        			if ($cigo =~ m/(\d+)H$/){
        			    $H2 = $1;
        			    $cigo =~ s/$H2\H$//;
        			}

        			#### calculating length of leftmost softclip
        			if ($cigo =~ m/^(\d+)S/){
        			    $S1 = $1;
        			} 

        			#### calculating length of rightmost softclip
        			if ($cigo =~ m/(\d+)S$/){
        			    $S2 = $1;

        			}
					
					if ($flag & 0x10){
				    $seq  =  reverse $seq;
					$seq =~ tr/ATGCatgc/TACGtacg/;

					$leftclip =  $H2 + $S2 ; 
					$rightclip =  $H1 + $S1 ;

					} else {
					
					$leftclip =  $H1 + $S1 ; 
					$rightclip =  $H2 + $S2 ;
					
					}
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
					print pos($seq), "\n";
					$CGhash{(pos($seq)+$context_offset)} = $counter_CG; # 1 based position
                }
				my $CGn = scalar @CGpos;
				print "CCCCCCCCCCCCCC\n";
				##### find base positions
				my @Cpos;
				while($seq =~ m/C/g){
                    push @Cpos, pos($seq); # 1 based position in read
					print pos($seq), "\n";
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
					my @summary_modpos;
					my @summary_modvalues;
					my @summary_modpos_positives;

					my @summary_modpos_within;
					my @summary_modvalues_within;
					my @summary_modpos_positives_within;

					if ($CGn > 0 && $mm) {

						my @premm = split(";", $mm);                    # (C+h,18,1) , (C+m,18,1) , ()	
						my @posmm = startchar(@premm, $mod);            # array with index of the correct element in @premm (in this case 1!)

						die "Fatal Error: Too many fields in MM which match with $mod at: $bline" if $#posmm > 0;
						
						$premm[$posmm[0]] =~ s/\Q$mod,//;                # (18,1)
						my @splitmm = split(",",  $premm[$posmm[0]]);    # (18), (1)
						my @splitml = split(",", $ml);                   # (0),(1),(89),(7)
						print "$readid ", join("::", %CGhash), " at test $bline \n";

						for my $sm (0..$#splitmm){
							##### l'elemento di splitmm fondamentalmente è l'indice del vettore con le C. solo che è incrementale, e quindi, a tutti quelli successivi al primo va aggiunto il precedente + 1 (perchè 0 vuol dire la successiva)
							if ($sm > 0){$splitmm[$sm] = $splitmm[$sm]+$splitmm[$sm-1]+1 } #increment
							my $Cposval = $Cpos[$splitmm[$sm]];
							
							die "Fatal error: less C in seq compared to MM field at: $bline" unless $Cposval;

							if (exists($CGhash{$Cposval})){
								
								my $within = "false";
								if ($Cposval > $leftclip & $Cposval < length($seq)-$rightclip){
									$within = "true"; 
								}

								my $modvalue = $splitml[(($sm * (scalar @premm))+$posmm[0])] ;
								if( $modvalue >= $threshold ){
									$Nmod++;
									@summary_modpos_positives = (@summary_modpos_positives, $Cposval);
									#### HERE
									if ($within eq "true"){ @summary_modpos_positives_within = (@summary_modpos_positives_within, $Cposval); }	

								} # prendo l'indice degli slementi in MM (che sarebbe $sm), lo moltiplico per lo scalar di @premm (che indica quante possibili mod ci sono ) perchè ML ha serie di N numeri dove N è il numero di possibili mod. poi ci aggiungo l'indice della mod in MM che rappresenta di "quanto shiftare a destra" quando si sceglie l'indice di ML.
								delete($CGhash{$Cposval});
								
								@summary_modvalues = (@summary_modvalues, $modvalue);
								@summary_modpos = (@summary_modpos, $Cposval);

								if ($within eq "true"){ 
									@summary_modvalues_within = (@summary_modvalues_within, $modvalue);
									@summary_modpos_within = (@summary_modpos_within, $Cposval);
								}

								
							} 

						}
						if (scalar(keys(%CGhash)) > 0){print "Too many CG:", scalar(keys(%CGhash)), "positioned as follows ", join("::", %CGhash), " at $bline \n"}

					} elsif ($CGn > 0 && !($mm)){
						$printmod = join("\t",@splitb[(0,2,3)]) . "\t$CGn\tMM_Undetected\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
					}
					#### columns: c("mod_readID", "mod_chr", "mod_pos", "CGn", "Nmod", "Cn", "mod_pos_Cmod", "mod_pos_allCG", "modvalues_CG" )
					
					my $mod_pos_Cmod = "NONE";
					if ((scalar @summary_modpos_positives) > 0){
						$mod_pos_Cmod = join(",",@summary_modpos_positives);
					}

					my $mod_pos_allCG = "NONE";
					if ((scalar @summary_modpos) > 0){
						$mod_pos_allCG = join(",",@summary_modpos);
					}

					my $modvalues_CG = "NONE";
					if ((scalar @summary_modvalues) > 0){
						$modvalues_CG = join(",",@summary_modvalues);
					}

					#### within
					my $mod_pos_Cmod_within = "NONE";
					if ((scalar @summary_modpos_positives_within) > 0){
						$mod_pos_Cmod_within = join(",",@summary_modpos_positives_within);
					}

					my $mod_pos_allCG_within = "NONE";
					if ((scalar @summary_modpos_within) > 0){
						$mod_pos_allCG_within = join(",",@summary_modpos_within);
					}

					my $modvalues_CG_within = "NONE";
					if ((scalar @summary_modvalues_within) > 0){
						$modvalues_CG_within = join(",",@summary_modvalues_within);
					}

					$printmod = join("\t",@splitb[(0,2,3)]) . "\t$CGn\t$Nmod\t$Cn\t$mod_pos_Cmod\t$mod_pos_allCG\t$modvalues_CG\t$mod_pos_Cmod_within\t$mod_pos_allCG_within\t$modvalues_CG_within\t$flag\t$leftclip\t$len\t$rightclip\n";

					open(OUT, ">>$modoutdir/$cropname.$outmod{$mod}");
					print OUT $printmod;
					close(OUT);

				}
				###### copia fino a quì
			} 
		}


		if ($check_found==0){
			for my $mod (@modlist){
				# PRINTA LE BARACCHE ETC
				open(OUT, ">>$modoutdir/$cropname.$outmod{$mod}");
				print OUT "$readid\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
				close(OUT);
			}
		}


	}
	if (exists $outset{"modbam"}){
		#open(OUT, ">>$modbamoutdir/$cropname\_mod_mappings.sam");
		my $cmd  = "samtools view -h  -O BAM " . " $modbamoutdir/$cropname\_mod_mappings.sam " . " -o " . " $modbamoutdir/$cropname\_mod_mappings.bam " ;
		system($cmd);
		$cmd = "rm " .  " $modbamoutdir/$cropname\_mod_mappings.sam " ;
		system($cmd);
	}

    $pm->finish(0, \$q);
}
$pm->wait_all_children;




exit
