#!/bin/bash
#SBATCH --job-name=_your_job_name_here_
#SBATCH --account=_your_account_name_here_
#SBATCH --partition=_your_partition_name_here_
#SBATCH --chdir=_your_working_directory_name_here_
#SBATCH --ntasks=1
#SBATCH --mem=65536M 
#SBATCH --cpus-per-task=16
#SBATCH --output=sbatch.out
#SBATCH --error=sbatch.err

OUTFILE="run_map.log"
REFSEQ="canu.contigs.fasta"
PGLBIN=/data/gpfs/assoc/pgl/bin
TRIMGALORE=$PGLBIN/TrimGalore/trim_galore

function log_error
{
    echo "$(date '+%F %T') ERROR: $*" >&2 |& tee -a $OUTFILE
}

function log_info
{
    echo "$(date '+%F %T') INFO: $*" | tee -a $OUTFILE
}

echo > $OUTFILE < /dev/null

R1=Ki11_S13_L002_R1_001.fastq.gz
R2=${R1/R1/R2}
log_info "Trimming reads file $R1 with TrimGalore"
# paired-end reads, so R1 and R2
$TRIMGALORE --paired --retain_unpaired --cores 16 --max_n 40 --gzip $R1 $R2 |& tee -a $OUTFILE
if [ $? != 0 ]; then
    log_error "TrimGalore exited non-zero[$?]"
    exit 1
fi
R1=${R1/001.fastq/001_val_1.fq}
R2=${R1/_R1_001_val_1/_R2_001_val_2}
for F in $R1 $R2; do
    if [ ! -e $F ]; then
        log_error "TrimGalore output file $F not found"
        exit 1
    fi
done

log_info "Index reference (bwa index)" 
bwa index -a bwtsw $REFSEQ |& tee -a $OUTFILE
if [ $? != 0 ]; then
    log_error "bwa index exit non-zero[$?]"
    exit 1
fi
log_info echo "Generate reference FAIDX (samtools faidx)"
samtools faidx $REFSEQ |& tee -a $OUTFILE
if [ $? != 0 ]; then
    log_error "samtools faidx exited non-zero[$?]"
    exit 1
fi


for RFILE in Ki11_S13_L002_R[12]_001_val_[12].fq.gz; do
    log_info "Aligning reads from $RFILE to reference $REFSEQ" 
    ALIGN_OUT="$(basename -s.fq.gz $RFILE).sai"
    bwa aln -t 24 -f $ALIGN_OUT $REFSEQ $RFILE |& tee -a $OUTFILE
    if [ $? != 0 ]; then
        log_error "bwa aln exited non-zero[$?]"
        exit 1
    fi
    if [ ! -s $ALIGN_OUT ]; then
        log_error "bwa aln output file $ALIGN_OUT does not exist or is zero-length"
        exit 1
    fi
done

log_info "Merging results (bwa sampe)"
A1="$(basename -s.fq.gz $R1).sai"
A2="$(basename -s.fq.gz $R2).sai"
MERGE_FILE=Ki13_13-bwa_aln.sam
bwa sampe -f $MERGE_FILE $REFSEQ $A1 $A2 $R1 $R2 |& tee -a $OUTFILE
if [ $? != 0 ]; then
    log_error "bwa sampe exited non-zero[$?]"
    exit 1
fi
if [ ! -s $MERGE_FILE ]; then
    log_error "bwa sampe output file $MERGE_FILE does not exist or is zero-length"
    exit 1
fi

touch run.done
exit 0

