#!/bin/bash
#SBATCH --job-name=zmays-allhic
#SBATCH --account=cpu-s1-pgl-0
#SBATCH --partition=cpu-s1-pgl-0
#SBATCH --chdir=/data/gpfs/assoc/pgl/dcurdie/zmays-allhic
#SBATCH --ntasks=1
#SBATCH --mem=65536M 
#SBATCH --cpus-per-task=16
#SBATCH --output=sbatch.out
#SBATCH --error=sbatch.err

OUTFILE="run_allhic.log"
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

##
# Prune step
BAM=Ki13_S13-bwa_aln.clean.bam
ALLTBL=Allele.ctg.table
log_info "Pruning $BAM..."
$ALLHIC/bin/ALLHiC_prune -i $ALLTBL -b $BAM -r $REFSEQ | tee -a $OUTFILE
if [ $? != 0 ]; then
    log_error "ALLHiC_prune exited non-zero[$?]"
    exit 1
fi
if [ ! -s prunning.bam ]; then
    log_error "ALLHiC_prune output file prunning.bam does not exist or is zero-length"
    exit 1
fi

##
# Partition step
RE_SITES="GATCGATC,GAATGATC,GATTGATC,GACTGATC,GAGTGATC,GAATAATC,GAATATTC,GAATACTC,GAATAGTC,GTATAATC,GTATATTC,GTATACTC,GTATAGTC,GCATAATC,GCATATTC,GCATACTC,GCATAGTC,GGATAATC,GGATATTC,GGATACTC,GGATAGTC,GATCAATC,GATCATTC,GATCACTC,GATCAGTC"
NOCTGS=10
log_info "Partitioning pruned BAM "
$ALLHIC/bin/ALLHiC_partition -b prunning.bam -r $REFSEQ -e $RE_SITES -k $NOCTGS |& tee -a $OUTFILE
if [ $? != 0 ]; then
    log_error "ALLHiC_partition exited non-zero[$?]"
    exit 1
fi
for FILE_NAME in prunning.{clm,distribution.txt,pairs.txt,counts_${RE_SITES//,/_}.txt} ; do
    if [ ! -s $FILE_NAME ]; then
        log_error "ALLHiC_partition output file $FILE_NAME does not exist or is zero-length"
        exit 1
    fi
done
## HACK: need to fix ALLHiC_partition or allhic so that output files match
if [ -e prunning.counts_${RE_SITES//,/_}.txt -a ! -e prunning.counts_${RE_SITES}.txt ]; then
    ln -f -s prunning.counts_${RE_SITES//,/_}.txt prunning.counts_${RE_SITES}.txt
fi

##
# Rescue step (Skipping)
#log_info "Rescuing $BAM"
#CLUSTERS=
#COUNTS=
#$ALLHIC/bin/ALLHiC_rescue -i $BAM -r $REFSEQ -c $CLUSTERS -i $COUNTS |& tee -a $OUTFILE

##
# Optimize step
log_info "Extracting BAM using restriction enzyme sites $RE_SITES"
$ALLHIC/bin/allhic extract $BAM $REFSEQ --RE $RE_SITES |& tee -a $OUTFILE
if [ $? != 0 ]; then
    log_error "allhic extract exited non-zero[$?]"
    exit 1
fi

log_info "End of known outputs. Terminating script."
exit 0

if [ ! -s SOME_FILE_NAME_HERE ]; then
    log_error "allhic extract ouput file SOME_FILE_NAME_HERE does not exist or is zero-length"
    exit 1
fi
IDX=0
CLM=${BAM/bam/clm}
while [ $IDX -le $NOCTGS ]; do
    IDX=$(($IDX+1))
    log_info "Running optimize pass $IDX..."
    $ALLHIC/bin/allhic optimize group${IDX}.txt $CLM |& tee -a $OUTFILE
    if [ $? != 0 ]; then
        log_error "allhic optimize exited non-zero[$?]"
        exit 1
    fi
    if [ ! -s SOME_FILE_NAME_HERE ]; then
        log_error "allhic optimize output file SOME_FILE_NAME_HERE does not exist or is zero-length"
        exit 1
    fi
done

##
# Build step
log_info "Building group locations"
$ALLHIC/bin/ALLHiC_build $REFSEQ |& tee -a $OUTFILE
if [ $? != 0 ]; then
    log_error "ALLHiC_build exited non-zero[$?]"
    exit 1
fi
if [ ! -s SOME_FILE_NAME_HERE ]; then
    log_error "ALLHiC_build output file SOME_FILE_NAME_HERE does not exist or is zero-length"
    exit 1
fi

##
# Plot step
GROUPS=
CHRN=
SIZE=500k
FORMAT=pdf
log_info "Generating the chromatin contact matrix to evaluate genome scaffolding"
$ALLHIC/bin/ALLHiC_plot $BAM $GROUPS $CHRN $SIZE $FORMAT |& tee -a $OUTFILE
if [ $? != 0 ]; then
    log_error "ALLHiC_plot exited non-zero[$?]" 
    exit 1
fi
if [ ! -s SOME_FILE_NAME_HERE ]; then
    log_error "ALLHiC_plot output file SOME_FILE_NAME_HERE does not exist or is zero-length"
    exit 1
fi

log_info "Run complete"

touch run.done
exit 0

