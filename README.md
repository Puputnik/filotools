# filotools
Bioinformatic suite for Nanopore cfDNA analysis

The pipeline is optimized to work with fastq obtained from Nanopore multiplex DNA runs (NBD114.24 or NBD114.96), basecalled with dorado; but it's adaptable to older chemistries or basecallers.

Dorado can produce two types of outputs:

- Single bam files: 
  All the reads are stored in a single bam file. these files are usually pretty large (especially for promethion experiments) and the demultiplexing handles them badly.
  Split the bam file in multiple .fastq files, and store them in a folder named "fastq"
```
  mkdir fastq
  cd fastq
  samtools fastq /path/to/modcalls.bam | split -l 4000 --additional-suffix=.fastq
```
  where "/path/to/modcalls.bam" is the path in which is stored the .bam file produced by dorado.


- Splitted bam files (typically if the basecalling is performed during the run via MinKNOW):
  can be either:
    - Demultiplexed folders (bam_pass/barcode*/*bam) if "Barcoding" was set ON during the run.
   ```   
      mkdir fastq
      for i in /path/to/bam_pass/*/*bam ; do samtools fastq $i > fastq/$(basename $i).fastq
   ```   
    - Mixed bam files (bam_pass/*bam) if "Barcoding" was set OFF during the run
   ```   
      mkdir fastq
      for i in /path/to/bam_pass/*bam ; do samtools fastq $i > fastq/$(basename $i).fastq
```
  where "/path/to/bam_pass/" is the folder in which MinKNOW is saving alignments.

If you plan to perform methylation analysis you'll need aligned bams (generated using "Alignment ON" during a MinKNOW run, or specify the --reference flag if running dorado offline) 
please use GCF_000001405.39_GRCh38.p13 as reference genome (you can download it here https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.39/ ).


the pipeline requires basecalled fastqs
