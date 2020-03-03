#!/bin/bash
#SBATCH --job-name=zmays-allhic-map
#SBATCH --account=cpu-s1-pgl-0
#SBATCH --partition=cpu-s1-pgl-0
#SBATCH --chdir=/data/gpfs/assoc/pgl/dcurdie/zmays-allhic
#SBATCH --ntasks=1
#SBATCH --mem=65536M 
#SBATCH --cpus-per-task=16
#SBATCH --output=sbatch.out
#SBATCH --error=sbatch.err

OUTFILE="run_filter.log"
REFSEQ="canu.contigs.fasta"

SCRATCH=/data/gpfs/home/dcurdie/scratch
PERL=$SCRATCH/miniconda-3.7/bin/perl
ALLHIC=$SCRATCH/ALLHiC

export PATH=$PATH:$ALLHIC/bin:$ALLHIC/scripts

function log_error
{
    echo "$(date '+%F %T') ERROR: $*" >&2 |& tee -a $OUTFILE
}

function log_info
{
    echo "$(date '+%F %T') INFO: $*" | tee -a $OUTFILE
}

echo > $OUTFILE < /dev/null

ALN=Ki13_S13-bwa_aln.sam
BASE=$(basename -s.sam $ALN)
log_info "Preprocess merged alignments (PreprocessSAMs.pl)"
$PERL $ALLHIC/scripts/PreprocessSAMs.pl $ALN $REFSEQ MBOI |& tee -a $OUTFILE
if [ $? != 0 ]; then
    log_error "PreprocessSAMs.pl exited non-zero[$?]"
    exit 1
fi
for F in ${BASE}.bam ${REFSEQ}.{near_GATC.500.bed,pos_of_GATC.txt}; do
    if [ ! -s $F ]; then
        log_error "PreprocessSAMs.pl output file $F does not exist or is zero-length"
        exit 1
    fi
done

log_info "Filter merged alignments (filterBAM_forHiC.pl)"
$PERL $ALLHIC/scripts/filterBAM_forHiC.pl ${BASE}.REduced.paired_only.bam ${BASE}.clean.sam |& tee -a $OUTFILE
if [ $? != 0 ]; then
    log_error "filterBAM_forHiC.pl exited non-zero[$?]"
    exit 1
fi
for F in ${BASE}.REduced.{bam,paired_only.bam,paired_only.flagstat} ${BASE}.clean.sam; do
    if [ ! -s $F ]; then
        log_error "filterBAM_forHiC.pl output file $F does not exist or is zero-length"
        exit 1
    fi
done

log_info "Convert to BAM (samtools view)"
samtools view -o ${BASE}.clean.bam -bt ${REFSEQ}.fai ${BASE}.clean.sam |& tee -a $OUTFILE
if [ $? != 0 ]; then
    log_error "samtools view exited non-zero[$?]"
    exit 1
fi
if [ ! -s ${BASE}.clean.bam ]; then
    log_error "samtools output file ${BASE}.clean.bam does not exist or is zero-length"
    exit 1
fi

touch run.done
exit 0

