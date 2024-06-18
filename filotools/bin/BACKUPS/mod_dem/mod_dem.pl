#!/usr/bin/perl
use Parallel::ForkManager;
use strict;
use warnings;
use Cwd qw(abs_path);
use Cwd qw(cwd);
use File::Spec;
use Getopt::Long;

my %inmod =  ( 
	               C+m   => 'C+m', 
    	           5mc   => 'C+m', 
    	           5mC   => 'C+m', 
  
				   C+h   => 'C+h',
	               5hmc  => 'C+h',
	               5hmC  => 'C+h');
my %outmod =  ( 
	               C+h   => '5hmC', 
    	           C+m   => '5mC'); 

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
GetOptions(
	"modbam|b=s" => \@modbams
	#"mods|m=s" => \@modlist,
	#"exclude_header|h=s" => $exc_head
	); 

die "Missing modbam path, use -b\n" unless @modbams; 

my @header;
my %bam;

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


#my $bampath = "/tank/LIQUID_BIOPSY_GRANT/MEGALODON/L02_FAT85244_count0/mod_mappings.bam";
#my $bampath = shift @ARGV;

#"(File::Spec->splitpath($bampath))[1]/"

#my @stats = qw( /tank/USB3/LB_ANALYSIS_2022/BAM/STATS/PLB03-T00-1C00-B01-L02-FAT85244-1.stats);
my @stats = qw(
/tank/USB3/LB_ANALYSIS_2022/BAM/STATS/PLB03-T00-1C00-B01-L02-FAT85244-0.stats
/tank/USB3/LB_ANALYSIS_2022/BAM/STATS/PLB03-T00-1C00-B02-L02-FAT85244-0.stats
/tank/USB3/LB_ANALYSIS_2022/BAM/STATS/PLB03-T00-2B00-B03-L02-FAT85244-0.stats
/tank/USB3/LB_ANALYSIS_2022/BAM/STATS/PLB03-T00-2C00-B04-L02-FAT85244-0.stats
/tank/USB3/LB_ANALYSIS_2022/BAM/STATS/PLB03-T00-2D00-B05-L02-FAT85244-0.stats
/tank/USB3/LB_ANALYSIS_2022/BAM/STATS/PLB03-T00-EX00-B12-L02-FAT85244-0.stats
/tank/USB3/LB_ANALYSIS_2022/BAM/STATS/PLB16-T00-1100-B06-L02-FAT85244-0.stats
/tank/USB3/LB_ANALYSIS_2022/BAM/STATS/PLB17-T00-1100-B07-L02-FAT85244-0.stats
/tank/USB3/LB_ANALYSIS_2022/BAM/STATS/PLB18-T00-1100-B08-L02-FAT85244-0.stats
/tank/USB3/LB_ANALYSIS_2022/BAM/STATS/PLB19-T00-1100-B09-L02-FAT85244-0.stats
/tank/USB3/LB_ANALYSIS_2022/BAM/STATS/PLB20-T00-1100-B10-L02-FAT85244-0.stats
/tank/USB3/LB_ANALYSIS_2022/BAM/STATS/PLB21-T00-1100-B11-L02-FAT85244-0.stats);
print @stats;

my $forks=12;
my $pm = Parallel::ForkManager->new($forks);
 
$pm->run_on_finish( sub {
    my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $out) = @_;
    #my $q = $data_structure_reference->{input};
    #$results{$q} = $data_structure_reference->{result};
	#print @_;
	print "done, ", join(',', ($pid, $exit_code, $exit_signal, $core_dump, $$out)), "\n"; 
});
 

for my $mod (@modlist){
	$mod = $inmod{$mod}; 
}  #= ("C+m");


foreach my $q (@stats) {
    my $pid = $pm->start and next;
    #print "ciao mamma guarda come parso $q\n";
	$q = abs_path($q);
	open(STATS, $q);
#	my $c = 0;
	while (my $sline = <STATS>){
		chomp $sline;
		my @splits = split("\t", $sline);
		my $printout = "$splits[3]\tNA\tNA\tNA\n";
		for my $bampath (@modbams){

			my ($idMM, $idML) = split(","; $bam{$bampath}{"splitpos"});
			

			if (exists($bam{$bampath}{"$splits[3]"})){

				my $bline = $bam{$bampath}{"$splits[3]"};
				#####copia da quì
				my @splitb = split("\t", $bline);
				
				###### faccio quì la questione della stringa? direi di sì in ogni caso mi serve per sapere 1) quante CG, 2) quante C
				#######if ($flag & 0x10){
				#######   $seq  =  reverse $seq;
				#######	$seq =~ tr/ATGCatgc/TACGtacg/;
				#######} 


				my $mm $splitb[$idMM];
				my $ml $splitb[$idML];

				$mm =~ s/MM:Z://;
				$ml =~ s/ML:B:C,//;

				for my $mod (@modlist){
				###### siamo rimasti quì, è da andare avanti con questa questione volevo fare if ($mm) {} else {aggiungi al printing l'output vuoto (tipo 0?)}
				###### cerchiamo l'indice della mod, e nell'ml il valore di prob sarà uguale al numero di elementi/2+indice della mod? no ma che cazzo dico, sarà i numeri pari compresi tra 0 e il numero di elementi-1+indicemod.
#					MM:Z:   ML:B:C
#   MM:Z:C+h,14,0,2;C+m,14,0,2;     ML:B:C,6,44,55,125,193,155
				}

				###### copia fino a quì
			} 
#		$c++;
#		print "ciao $c $q $sline\n" if $c == 1;
		}
	}
    $pm->finish(0, \$q);
}
$pm->wait_all_children;




exit