# Introduction
This workflow is intended to wrap the functionality documented in the [ALLHiC wiki](https://github.com/tangerzhang/ALLHiC/wiki).

This functionality is broken into 3 main sections: mapping, filtering, and allhic. This is covered by 3 scripts: run_map.sh, run_filter.sh, and run_allhic.sh. Each script does need some customization by setting a few variables, most of which can be found towards the top of the script. The key inputs that need to be supplied are a reference assembly and the HiC reads (FASTQ format). The HiC reads are assumed to be paired-end. The FASTQ files are assumed to be GZIP compressed. Additional inputs for each script are generally provided by the script run beforehand. The order to run the scripts in is run_map.sh, run_filter.sh, and finally run_allhic.sh. Detail for each script follows. For more general information consult the aforementioned [ALLHiC wiki](https://github.com/tangerzhang/ALLHiC/wiki).

# Script Details

## run_map.sh
This script needs to know the reference assembly file (FASTA) and HiC paired-end reads files (FASTQ). It will begin by attempting to trim the reads using a tool such as [TrimGalore](https://github.com/FelixKrueger/TrimGalore). Next it attempts to index the reads using [BWA](https://github.com/lh3/bwa) and [samtools](https://github.com/samtools/samtools). Next the HiC reads are mapped, independently, to the reference assembly using [BWA](https://github.com/lh3/bwa) and the aln method. Finally, the independent mapping SAM files are merged using [BWA](https://github.com/lh3/bwa) and the sampe method.

You will need to ensure that bwa and samptools can be found within PATH, or else update the command invocations. You will also need to define where to find the trimming tool, as well as tailor the command-line arguments for it.

## run_filter.sh
This script will take the merged mapping SAM file from the prior script, run a pair of pre-processing script to filter down the set of reads data before converting the mapping SAM file to BAM. Needed are the [ALLHiC scripts](https://github.com/tangerzhang/ALLHiC), as well as [BWA](https://github.com/lh3/bwa) and [samtools](https://github.com/samtools/samtools).

You will need to adjust the file name of the merged mapped reads. You will also need to define where to find the ALLHiC scripts. As in the prior script bwa and samtools will need to be found in PATH.

## run_allhic.sh
Now the actual processing with ALLHiC begins. The first step in this script is to run ALLHiC_prune against a [constructed allele table](https://github.com/tangerzhang/ALLHiC/wiki/ALLHiC:-identify-allelic-contigs) to remove alleles and weak signals from the cleaned BAM file from the previous script. This produces a pruned BAM file. Next the reference asembly and the pruned BAM file are used with several restriction enzyme cut sites to partition up the contigs into K groups (K should be the number of known or anticipated chromosomes). Next the groups are ordered and oriented. Finally, the tour outputs are converted to fasta sequences and an agp mapping file. An optional step to generate a contact map in PDF format is available, but requires a file listing the group names and lengths.

