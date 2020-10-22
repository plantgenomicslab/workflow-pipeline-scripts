#!/bin/sh
#SBATCH --job-name=opuntia-repeatannotation
#SBATCH --account=cpu-s1-pgl-0
#SBATCH --partition=cpu-s1-pgl-0
#########SBATCH --chdir=/data/gpfs/assoc/pgl/dcurdie/opuntia-repeatannotation
#SBATCH --ntasks=1
####SBATCH --mem=65536M
#SBATCH --mem=120g
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=ALL
#SBATCH --mail-user=wyim@unr.edu
#SBATCH --output=sbatch.out
#SBATCH --error=sbatch.err
 #conda activate base
##
# Based upon the process documented in:
# http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/Repeat_Library_Construction-Advanced
##

##
# How many node CPUs to use
NCPU=32
# Import file system locations
# Root of all data
PGL_ROOT=/data/gpfs/assoc/pgl
# Home
GPFS_HOME=/data/gpfs/home
# Scratch
SCRATCH=/data/gpfs/home/wyim/scratch/
# Working directory
WORKDIR=/data/gpfs/home/wyim/scratch/data/iceplant/repeat
#WORKDIR=$SCRATCH/data/opuntia_coche/repeat_annotation
# Do we need some temp space..?
TMPDIR=$SCRATCH/tmp
# Shared reference data
#  All transposase protein database: http://www.hrt.msu.edu/uploads/535/78637/Tpases020812.gz
#   now repbase20.05_aaSeq_cleaned_TE.fa from wyim, based on Repbase
#  DNA transposase protein database: http://www.hrt.msu.edu/uploads/535/78637/Tpases020812DNA.gz
#   now repbase20.05_aaSeq_classII_TE.fa from wyim, based on RepBase
#  Plant protein database: http://www.hrt.msu.edu/uploads/535/78637/alluniRefprexp070416.gz
#   now uniprot_sprot.fasta from UniProtKB https://www.uniprot.org/downloads
DATADIR=$SCRATCH/bin/REPET/dataref
##
# Location of important tools / binaries
 Location of important tools / binaries
 PERL=$SCRATCH/bin/miniconda3/bin/perl
 DIR_MITE=$SCRATCH/bin/REPET/MITE-Tracker
 DIR_GT=$SCRATCH/bin/gt-1.5.10-Linux_x86_64-64bit-complete/bin  # http://genometools.org/ version 1.5.5
 DIR_CRL=$SCRATCH/bin/REPET/Custom-Repeat-Library  # https://github.com/plantgenomicslab/Custom-Repeat-Library
 DIR_RM=$SCRATCH/bin/REPET/RepeatMasker-4.1.0/RepeatMasker # http://www.repeatmasker.org/
 DIR_BLAST=$SCRATCH/bin/ncbi-blast-2.2.28+/bin/  # NCBI BLAST+
 DIR_RD=$SCRATCH/bin/REPET/RepeatModeler-2.0.1/ # http://www.repeatmasker.org/RepeatModeler.html
 DIR_PE=$SCRATCH/bin/REPET/ProtExcluder1.2  # http://www.hrt.msu.edu/uploads/535/78637/ProtExcluder1.2.tar.gz
 export PATH=$PATH:$DIR_PE
 export VSEARCH_BINDIR=$SCRATCH/bin/
 export NSEG_COMMAND=$SCRATCH/bin/REPET/binary/bin/nseg  # this is needed for RepeatScout filter-stage-1.prl script
##

#
# Let's get to work...
pushd $WORKDIR > /dev/null
genome=iceplant
seqfile1=$WORKDIR/genome.fa
if [ ! -e $seqfile1 ]; then
    echo "ERROR: Sequence file $seqfile1 not found" >&2
    exit 1
fi
fasta_rename.py  -i genome.fa --pre $genome -o genome_rename.fa
seqfile=$WORKDIR/genome_rename.fa
if [ ! -e $seqfile ]; then
    echo "ERROR: Rename $seqfile not found" >&2
    exit 1
fi
seqfileindex=SequenceIndex
##
# Run MITE-hunter
# - where is it? Er, use MITE-Tracker instead: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2376-y
#  MITE-Tracker needs to know where to find vsearch
seqfile=$(realpath $seqfile)
if [ ! -d $DIR_MITE ]; then
    echo "ERROR: MITE Tracker directory $MITE_DIR not found" >&2
    exit 1
fi
if [ ! -e $DIR_MITE/MITETracker.py ]; then
    echo "ERROR: MITE Tracker script not found in $MITE_DIR" >&2
    exit 1
fi
mkdir -p repeats/results/$genome
ln -f -s $(pwd)/repeats/results/$genome $DIR_MITE/results/$genome
pushd repeats > /dev/null
python $DIR_MITE/MITETracker.py -g $seqfile -j $genome -w $NCPU
if [ "$?" != "0" ]; then
    echo "ERROR: MITE Tracker exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -s results/$genome/all.fasta ]; then
    echo "ERROR: MITE Tracker output file results/$genome/all.fasta not found" >&2
    exit 1
fi
cp results/$genome/all.fasta > ../MITE.lib
popd >/dev/null 2>&1

##
# Classify LTR retrotransposons

# First, collection of candidate elements with LTRs that are 99% or more in similarity using LTRharvest
$DIR_GT/gt suffixerator -db $seqfile -indexname $seqfileindex -tis -suf -lcp -des -ssp -dna
$DIR_GT/gt suffixerator -ssp -des -dna -indexname $seqfileindex -db $seqfile -suf -lcp
if [ "$?" != "0" ]; then
    echo "ERROR: Genometools suffixerator exited non-zero[$?]" >&2
    exit 1
fi
for F in SequenceIndex.{des,esq,lcp,llv,md5,prj,sds,suf}; do
    if [ ! -e $F ]; then
        echo "ERROR: Genometools suffixerator output file $F not found" >&2
        exit 1
    fi
done
$DIR_GT/gt ltrharvest -index $seqfileindex -out ${seqfile}.out99 -outinner ${seqfile}.outinner99 \
     -gff3 ${seqfile}.gff99 -minlenltr 100 -maxlenltr 6000 -mindistltr 1500 -maxdistltr 25000 \
     -mintsd 5 -maxtsd 5 -motif tgca -similar 99 -vic 10 > ${seqfile}.result99
if [ "$?" != "0" ]; then
    echo "ERROR: Genometools ltrharvest exited non-zero[$?]" >&2
    exit 1
fi
for F in ${seqfile}.{out,outinner,gff}99; do
    if [ ! -e $F ]; then
        echo "ERROR: Genometools ltrharvest output file $F not found" >&2
        exit 1
    fi
done

# Second, using LTRdigest to find elements with PPT (poly purine tract) or PBS (primer binding site)
$DIR_GT/gt gff3 -sort ${seqfile}.gff99 > ${seqfile}.gff99.sort
if [ "$?" != "0" ]; then
    echo "ERROR: Genometools gff3 exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -s ${seqfile}.gff99.sort ]; then
    echo "ERROR: Zero-length output from genometools gff3" >&2
    exit 1
fi
$DIR_GT/gt ltrdigest -trnas $DATADIR/eukaryotic-tRNAs.fa ${seqfile}.gff99.sort ${seqfileindex} > ${seqfile}.gff99.dgt
if [ "$?" != "0" ]; then
    echo "ERROR: Genometools ltrdigest exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -s ${seqfile}.gff99.dgt ]; then
    echo "ERROR: Zero-length output from genometools ltrdigest" >&2
    exit 1
fi

# make sure BioPerl is active/available...
$PERL $DIR_CRL/CRL_Step1.pl --gff ${seqfile}.gff99.dgt
if [ "$?" != "0" ]; then
    echo "ERROR: CRL_Step1.pl exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -e CRL_Step1_Passed_Elements.txt ]; then
    echo "ERROR: CRL_Step1.pl output file CRL_Step1_Passed_Elements.txt not found" >&2
    exit 1
fi

# Third, further filtering of the candidate elements
$PERL $DIR_CRL/CRL_Step2.pl --step1 CRL_Step1_Passed_Elements.txt --repeatfile ${seqfile}.out99 \
    --resultfile ${seqfile}.result99 --sequencefile $seqfile --removed_repeats CRL_Step2_Passed_Elements.fasta
if [ "$?" != "0" ]; then
    echo "ERROR: CRL_Step2.pl exited non-zero[$?]" >&2
    exit 1
fi
for F in {CRL_Step2_Passed_Elements,Repeat_down1,Repeat_up1}.fasta; do
    if [ ! -e $F ]; then
        echo "ERROR: CRL_Step2.pl output file $F not found" >&2
        exit 1
    fi
done

mkdir fasta_files
mv Repeat_*.fasta fasta_files
mv CRL_Step2_Passed_Elements.fasta fasta_files
pushd fasta_files > /dev/null

# this requires MUSCLE to be available
$PERL $DIR_CRL/CRL_Step3.pl --directory $(pwd) --step2 CRL_Step2_Passed_Elements.fasta --pidentity 60 --seq_c 25
if [ "$?" != "0" ]; then
    echo "ERROR: CRL_Step3.pl exit non-zero[$?]" >&2
    exit 1
fi
[ -e CRL_Step3_Passed_Elements.fasta ] && mv CRL_Step3_Passed_Elements.fasta ..
rm -f Repeat_*.fasta
popd > /dev/null
if [ ! -e CRL_Step3_Passed_Elements.fasta ]; then
    echo "ERROR: CRL_Step3.pl output file CRL_Step3_Passed_Elements.fasta not found" >&2
    exit 1
fi
if [ ! -s CRL_Step3_Passed_Elements.fasta ]; then
    echo "ERROR: Zero-length output from CRL_Step3.pl" >&2
    exit 1
fi

# Fourth, identify elements with nested insertions
$PERL $DIR_CRL/ltr_library.pl --resultfile ${seqfile}.result99 --step3 CRL_Step3_Passed_Elements.fasta --sequencefile $seqfile
if [ "$?" != "0" ]; then
    echo "ERROR: ltr_library.pl exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -e lLTR_Only.lib ]; then
    echo "ERROR: ltr_library.pl output file lLTR_Only.lib not found" >&2
    exit 1
fi
if [ ! -s lLTR_Only.lib ]; then
    echo "ERROR: Zero-length output from ltrlibrary.pl" >&2
    exit 1
fi
cat lLTR_Only.lib MITE.lib > repeats_to_mask_LTR99.fasta
$DIR_RM/RepeatMasker -pa $(($NCPU/4)) -lib repeats_to_mask_LTR99.fasta -nolow -dir . ${seqfile}.outinner99
if [ "$?" != "0" ]; then
    echo "ERROR: RepeatMasker exited non-zero[$?]" >&2
    exit 1
fi
for F in ${seqfile}.outinner99.{masked,out}; do
    if [ ! -e $F ]; then
        echo "ERROR: RepeatMasker output file $F not found." >&2
        exit 1
    fi
done

$PERL $DIR_CRL/cleanRM.pl ${seqfile}.outinner99.out ${seqfile}.outinner99.masked > ${seqfile}.outinner99.unmasked
if [ "$?" != "0" ]; then
    echo "ERROR: cleanRM.pl exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -s ${seqfile}.outinner99.unmasked ]; then
    echo "ERROR: Zero-length output file from cleanRM.pl" >&2
    exit 1
fi
$PERL $DIR_CRL/rmshortinner.pl ${seqfile}.outinner99.unmasked 50 > ${seqfile}.outinner99.clean
if [ "$?" != "0" ]; then
    echo "ERROR: rmshortinner.pl exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -s ${seqfile}.outinner99.clean ]; then
    echo "ERROR: Zero-length output from rmshortinner.pl" >&2
    exit 1
fi

if [ ! -e $DATADIR/repbase20.05_aaSeq_classII_TE.fa.pin ]; then
    $DIR_BLAST/makeblastdb -in $DATADIR/repbase20.05_aaSeq_classII_TE.fa -dbtype prot
    if [ "$?" != "0" ]; then
        echo "ERROR: makeblastdb exited non-zero[$?]" >&2
        exit 1
    fi
fi
$DIR_BLAST/blastx -query ${seqfile}.outinner99.clean -db $DATADIR/repbase20.05_aaSeq_classII_TE.fa -evalue 1e-10 -num_descriptions 10 -num_threads $NCPU \
    -out ${seqfile}.outinner99.clean_blastx.out.txt
if [ "$?" != "0" ]; then
    echo "ERROR: blastx exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -s ${seqfile}.outinner99.clean_blastx.out.txt ]; then
    echo "ERROR: Zero-length output from blastx" >&2
    exit 1
fi
$PERL $DIR_CRL/outinner_blastx_parse.pl --blastx ${seqfile}.outinner99.clean_blastx.out.txt --outinner  ${seqfile}.outinner99
if [ "$?" != "0" ]; then
    echo "ERROR: outinner_blastx_parse.pl exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -s passed_outinner_sequence.fasta ]; then
    echo "ERROR: Zero-length output from outinner_blastx_parse.pl" >&2
    exit 1
fi

# Fifth, building examplars
$PERL $DIR_CRL/CRL_Step4.pl --step3 CRL_Step3_Passed_Elements.fasta --resultfile ${seqfile}.result99 \
    --innerfile passed_outinner_sequence.fasta --sequencefile ${seqfile}
if [ "$?" != "0" ]; then
    echo "ERROR: CRL_Step4.pl exited non-zero[$?]" >&2
    exit 1
fi
for F in {Inner,lLTRs}_Seq_For_BLAST.fasta; do
    if [ ! -e $F ]; then
        echo "ERROR CRL_Step4.pl output file $F not found" >&2
        exit 1
    fi
done

$DIR_BLAST/makeblastdb -in lLTRs_Seq_For_BLAST.fasta -dbtype nucl
if [ "$?" != "0" ]; then
    echo "ERROR: makeblastdb exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -e lLTRs_Seq_For_BLAST.fasta.nin ]; then
    echo "ERROR: makeblastdb output file lLTRs_Seq_For_BLAST.fasta.pin not found" >&2
    exit 1
fi
$DIR_BLAST/blastn -query lLTRs_Seq_For_BLAST.fasta -db lLTRs_Seq_For_BLAST.fasta -evalue 1e-10 -num_descriptions 1000 -num_threads $NCPU \
    -out lLTRs_Seq_For_BLAST.fasta.out
if [ "$?" != "0" ]; then
    echo "ERROR: blastn exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -e lLTRs_Seq_For_BLAST.fasta.out ]; then
    echo "ERROR: blastn output file lLTRs_Seq_For_BLAST.fasta.out not found" >&2
    exit 1
fi
$DIR_BLAST/makeblastdb -in Inner_Seq_For_BLAST.fasta -dbtype nucl
if [ "$?" != "0" ]; then
    echo "ERROR: makeblastdb exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -e Inner_Seq_For_BLAST.fasta.nin ]; then
    echo "ERROR: makeblastdb output file Inner_Seq_For_BLAST.fasta.pin not found" >&2
    exit 1
fi
$DIR_BLAST/blastn -query Inner_Seq_For_BLAST.fasta -db Inner_Seq_For_BLAST.fasta -evalue 1e-10 -num_descriptions 1000 -num_threads $NCPU \
    -out Inner_Seq_For_BLAST.fasta.out
if [ "$?" != "0" ]; then
    echo "ERROR: blastn exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -e Inner_Seq_For_BLAST.fasta.out ]; then
    echo "ERROR: blastn output file Inner_Seq_For_BLAST.fasta.out not found" >&2
    exit 1
fi

$PERL $DIR_CRL/CRL_Step5.pl --LTR_blast lLTRs_Seq_For_BLAST.fasta.out --inner_blast Inner_Seq_For_BLAST.fasta.out \
    --step3 CRL_Step3_Passed_Elements.fasta --final LTR99.lib --pcoverage 90 --pidentity 80

# Cleanup before proceeding
rm -f CRL_Step1_Passed_Elements.txt CRL_Step2_Passed_Elements.fasta CRL_Step3_Passed_Elements.fasta lLTR_Only.lib \
    passed_outinner_sequence.fasta lLTRs_Seq_For_BLAST.fasta* Inner_Seq_For_BLAST.fasta*

##
# Collection of relatively old LTR retrotransposons

# Repeating Classify LTR above, but without motif and 99% similarity
$DIR_GT/gt ltrharvest -index $seqfileindex -out ${seqfile}.out85 -outinner ${seqfile}.outinner85 \
    -gff3 ${seqfile}.gff85 -minlenltr 100 -maxlenltr 6000 -mindistltr 1500 -maxdistltr 25000 \
    -mintsd 5 -maxtsd 5 -vic 10 > ${seqfile}.result85
if [ "$?" != "0" ]; then
    echo "ERROR: Genometools ltrharvest exited non-zero[$?]" >&2
    exit 1
fi
for F in ${seqfile}.{out,outinner,gff}85; do
    if [ ! -e $F ]; then
        echo "ERROR: Genometools ltrharvest output file $F not found" >&2
        exit 1
    fi
done
$DIR_GT/gt gff3 -sort ${seqfile}.gff85 > ${seqfile}.gff85.sort
if [ "$?" != "0" ]; then
    echo "ERROR: Genometools gff3 exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -s ${seqfile}.gff85.sort ]; then
    echo "ERROR: Zero-length output from genometools gff3" >&2
    exit 1
fi
$DIR_GT/gt ltrdigest -trnas $DATADIR/eukaryotic-tRNAs.fa ${seqfile}.gff85.sort ${seqfileindex} > ${seqfile}.gff85.dgt
if [ "$?" != "0" ]; then
    echo "ERROR: Genometools ltrdigest exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -s ${seqfile}.gff85.dgt ]; then
    echo "ERROR: Zero-length output from genometools ltrdigest" >&2
    exit 1
fi

$PERL $DIR_CRL/CRL_Step1.pl --gff ${seqfile}.gff85.dgt
if [ "$?" != "0" ]; then
    echo "ERROR: CRL_Step1.pl exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -e CRL_Step1_Passed_Elements.txt ]; then
    echo "ERROR: CRL_Step1.pl output file CRL_Step1_Passed_Elements.txt not found" >&2
    exit 1
fi

# further filtering of the candidate elements
$PERL $DIR_CRL/CRL_Step2.pl --step1 CRL_Step1_Passed_Elements.txt --repeatfile ${seqfile}.out85 \
    --resultfile ${seqfile}.result85 --sequencefile $seqfile --removed_repeats CRL_Step2_Passed_Elements.fasta
if [ "$?" != "0" ]; then
    echo "ERROR: CRL_Step2.pl exited non-zero[$?]" >&2
    exit 1
fi
for F in {CRL_Step2_Passed_Elements,Repeat_down1,Repeat_up1}.fasta; do
    if [ ! -e $F ]; then
        echo "ERROR: CRL_Step2.pl output file $F not found" >&2
        exit 1
    fi
done

mkdir -p fasta_files  # probably got removed by prior run of CRL_Step3.pl
mv Repeat_*.fasta fasta_files
mv CRL_Step2_Passed_Elements.fasta fasta_files
pushd fasta_files > /dev/null

# this requires MUSCLE to be available
$PERL $DIR_CRL/CRL_Step3.pl --directory $(pwd) --step2 CRL_Step2_Passed_Elements.fasta --pidentity 60 --seq_c 25
if [ "$?" != "0" ]; then
    echo "ERROR: CRL_Step3.pl exited non-zero[$?]" >&2
    exit 1
fi
[ -f CRL_Step3_Passed_Elements.fasta ] && mv CRL_Step3_Passed_Elements.fasta ..
rm -f Repeat_*.fasta
popd > /dev/null
if [ ! -e CRL_Step3_Passed_Elements.fasta ]; then
    echo "ERROR: CRL_Step3.pl output file CRL_Step3_Passed_Elements.fasta not found" >&2
    exit 1
fi
if [ ! -s CRL_Step3_Passed_Elements.fasta ]; then
    echo "ERROR: Zero-length output from CRL_Step3.pl" >&2
    exit 1
fi

# identify elements with nested insertions
$PERL $DIR_CRL/ltr_library.pl --resultfile ${seqfile}.result85 --step3 CRL_Step3_Passed_Elements.fasta --sequencefile $seqfile
if [ "$?" != "0" ]; then
    echo "ERROR: ltr_library.pl exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -e lLTR_Only.lib ]; then
    echo "ERROR: ltr_library.pl output file lLTR_Only.lib not found" >&2
    exit 1
fi
if [ ! -s lLTR_Only.lib ]; then
    echo "ERROR: Zero-length output from ltrlibrary.pl" >&2
    exit 1
fi
cat lLTR_Only.lib MITE.lib > repeats_to_mask_LTR85.fasta
$DIR_RM/RepeatMasker -pa $(($NCPU/4)) -lib repeats_to_mask_LTR85.fasta -nolow -dir . ${seqfile}.outinner85
if [ "$?" != "0" ]; then
    echo "ERROR: RepeatMasker exited non-zero[$?]" >&2
    exit 1
fi
for F in ${seqfile}.outinner85.{masked,out}; do
    if [ ! -e $F ]; then
        echo "ERROR: RepeatMasker output file $F not found." >&2
        exit 1
    fi
done

$PERL $DIR_CRL/cleanRM.pl ${seqfile}.outinner85.out ${seqfile}.outinner85.masked > ${seqfile}.outinner85.unmasked
if [ "$?" != "0" ]; then
    echo "ERROR: cleanRM.pl exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -s ${seqfile}.outinner85.unmasked ]; then
    echo "ERROR: Zero-length output from cleanPM.pl" >&2
    exit 1
fi

$PERL $DIR_CRL/rmshortinner.pl ${seqfile}.outinner85.unmasked 50 > ${seqfile}.outinner85.clean
if [ "$?" != "0" ]; then
    echo "ERROR: rmshortinner.pl exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -s ${seqfile}.outinner85.clean ]; then
    echo "ERROR: Zero-length output from rmshortinner.pl" >&2
    exit 1
fi

$DIR_BLAST/blastx -query ${seqfile}.outinner85.clean -db $DATADIR/repbase20.05_aaSeq_classII_TE.fa -evalue 1e-10 -num_descriptions 10 -num_threads $NCPU \
    -out ${seqfile}.outinner85.clean_blastx.out.txt
if [ "$?" != "0" ]; then
    echo "ERROR: blastx exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -s ${seqfile}.outinner85.clean_blastx.out.txt ]; then
    echo "ERROR: Zero-length output from blastx" >&2
    exit 1
fi

$PERL $DIR_CRL/outinner_blastx_parse.pl --blastx ${seqfile}.outinner85.clean_blastx.out.txt --outinner  ${seqfile}.outinner85
if [ "$?" != "0" ]; then
    echo "ERROR: outinner_blastx_parse.pl exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -s passed_outinner_sequence.fasta ]; then
    echo "ERROR: Zero-length output from outinner_blastx_parse.pl" >&2
    exit 1
fi

# building examplars
$PERL $DIR_CRL/CRL_Step4.pl --step3 CRL_Step3_Passed_Elements.fasta --resultfile ${seqfile}.result85 \
    --innerfile passed_outinner_sequence.fasta --sequencefile ${seqfile}
if [ "$?" != "0" ]; then
    echo "ERROR: CRL_Step4.pl exited non-zero[$?]" >&2
    exit 1
fi
for F in {Inner,lLTRs}_Seq_For_BLAST.fasta; do
    if [ ! -e $F ]; then
        echo "ERROR CRL_Step4.pl output file $F not found" >&2
        exit 1
    fi
done

$DIR_BLAST/makeblastdb -in lLTRs_Seq_For_BLAST.fasta -dbtype nucl
if [ "$?" != "0" ]; then
    echo "ERROR: makeblastdb exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -e lLTRs_Seq_For_BLAST.fasta.nin ]; then
    echo "ERROR: makeblastdb output file lLTRs_Seq_For_BLAST.fasta.nin not found" >&2
    exit 1
fi
$DIR_BLAST/blastn -query lLTRs_Seq_For_BLAST.fasta -db lLTRs_Seq_For_BLAST.fasta -evalue 1e-10 -num_descriptions 1000 -num_threads $NCPU \
    -out lLTRs_Seq_For_BLAST.fasta.out
if [ "$?" != "0" ]; then
    echo "ERROR: blastn exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -e lLTRs_Seq_For_BLAST.fasta.out ]; then
    echo "ERROR: blastn output file lLTRs_Seq_For_BLAST.fasta.out not found" >&2
    exit 1
fi
$DIR_BLAST/makeblastdb -in Inner_Seq_For_BLAST.fasta -dbtype nucl
if [ "$?" != "0" ]; then
    echo "ERROR: makeblastdb exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -e Inner_Seq_For_BLAST.fasta.nin ]; then
    echo "ERROR: makeblastdb output file Inner_Seq_For_BLAST.fasta.nin not found" >&2
    exit 1
fi
$DIR_BLAST/blastn -query Inner_Seq_For_BLAST.fasta -db Inner_Seq_For_BLAST.fasta -evalue 1e-10 -num_descriptions 1000 -num_threads $NCPU \
    -out Inner_Seq_For_BLAST.fasta.out
if [ "$?" != "0" ]; then
    echo "ERROR: blastn exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -e Inner_Seq_For_BLAST.fasta.out ]; then
    echo "ERROR: blastn output file Inner_Seq_For_BLAST.fasta.out not found" >&2
    exit 1
fi

$PERL $DIR_CRL/CRL_Step5.pl --LTR_blast lLTRs_Seq_For_BLAST.fasta.out --inner_blast Inner_Seq_For_BLAST.fasta.out \
    --step3 CRL_Step3_Passed_Elements.fasta --final LTR85.lib --pcoverage 90 --pidentity 80
#if [ "$?" != "0" ]; then
#    echo "ERROR: CRL_Step5.pl exited non-zero" >&2
#    exit 1
#fi
if [ ! -e LTR85.lib ]; then
    echo "ERROR: CRL_Step5.pl output file LTR85.lib not found" >&2
    exit 1
fi

# Cleanup before proceeding
rm -f CRL_Step1_Passed_Elements.txt CRL_Step2_Passed_Elements.fasta CRL_Step3_Passed_Elements.fasta lLTR_Only.lib \
    passed_outinner_sequence.fasta lLTRs_Seq_For_BLAST.fasta* Inner_Seq_For_BLAST.fasta*

# filter first set of LTRs
$DIR_RM/RepeatMasker -pa $(($NCPU/4)) -lib LTR99.lib -dir . LTR85.lib
if [ "$?" != "0" ]; then
    echo "ERROR: RepeatMasker exited non-zero[$?]" >&2
    exit 1
fi
for F in LTR85.lib.{masked,out}; do
    if [ ! -e $F ]; then
        echo "ERROR: RepeatMasker output file $F not found" >&2
        exit 1
    fi
done

$PERL $DIR_CRL/remove_masked_sequence.pl --masked_elements LTR85.lib.masked --outfile FinalLTR85.lib
if [ "$?" != "0" ]; then
    echo "ERROR: remove_masked_sequence.pl exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -e FinalLTR85.lib ]; then
    echo "ERROR: remove_masked_sequence.pl output file FinalLTR85.lib not found" >&2
    exit 1
fi

cat LTR99.lib FinalLTR85.lib > allLTR.lib

##
# Collecting repetitive sequences by RepeatModeler

cat allLTR.lib MITE.lib > allMITE_LTR.lib
$DIR_RM/RepeatMasker -pa $(($NCPU/4)) -lib allMITE_LTR.lib -dir . ${seqfile}
if [ "$?" != "0" ]; then
    echo "ERROR: RepeatMasker exited non-zero[$?]" >&2
    exit 1
fi
for F in ${seqfile}.{masked,out}; do
    if [ ! -e $F ]; then
        echo "ERROR: RepeatMasker output file $F not found" >&2
        exit 1
    fi
done

$PERL $DIR_CRL/rmaskedpart.pl ${seqfile}.masked 50 > umseqfile
if [ "$?" != "0" ]; then
    echo "ERROR: rmaskedpart.pl exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -s umseqfile ]; then
    echo "ERROR: Zero-length output file from rmaskedpart.pl" >&2
    exit 1
fi

$PERL $DIR_RD/BuildDatabase -name umseqfiledb -engine ncbi umseqfile
if [ "$?" != "0" ]; then
    echo "ERROR: BuildDatabase exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -e umseqfiledb.nin ]; then
    echo "ERROR: BuildDatabase output file umseqfile.nin not found" >&2
    exit 1
fi

$PERL $DIR_RD/RepeatModeler -pa $NCPU -database umseqfiledb >umseqfile.out
if [ "$?" != "0" ]; then
    echo "ERROR: RepeatModeler exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -s umseqfile.out ]; then
    echo "ERROR: RepeatModeler output file umseqfile.out does not exist or is zero-length" >&2
    exit 1
fi

echo "Done with RepeatModeler. Please update me with the appropriate directory name to continue."
RM_OUT_DIR=( RM_[0-9]*.[A-Z][a-z][a-z][A-Z][a-z][a-z][0-9]* )
if [ -d $RM_OUT_DIR ]; then
    echo "I think the RepeatModeler output dir is $RM_OUT_DIR"
else
    echo "ERROR: Unable to guess the RepeatModeler output dir" >&2
    exit 1
fi
exit 0
#RD_OUT_DIR=RM_80459.TueMar241837122020
CLASSIFIED=$(find . -name consensi.fa.classified)


##
# the file consensi.fa.classified will be in a RM_<PID>.<TIMESTAMP> subdirectory. need to find it...
$PERL $DIR_CRL/repeatmodeler_parse.pl --fastafile $CLASSIFIED --unknowns repeatmodeler_unknowns.fasta \
    --identities repeatmodeler_identities.fasta
#if [ "$?" != "0" ]; then
 #   echo "ERROR: repeatmodeler_parse.pl exited not-zero[$?]" >&2
 #   exit 1
#fi
#for F in repeatmodeler_{unknowns,identities}.fasta; do
 #   if [ ! -e $F ]; then
 #       echo "ERROR: repeatmodeler_parse.pl output file $F not found" >&2
 #       exit 1
 #   fi
#done

if [ ! -e $DATADIR/repbase20.05_aaSeq_cleaned_TE.fa.pin ]; then
    $DIR_BLAST/makeblastdb -in $DATADIR/repbase20.05_aaSeq_cleaned_TE.fa -dbtype prot
    if [ "$?" != "0" ]; then
        echo "ERROR: makeblastdb exited non-zero[$?]" >&2
        exit 1
    fi
fi
$DIR_BLAST/blastx -query repeatmodeler_unknowns.fasta -db $DATADIR/repbase20.05_aaSeq_cleaned_TE.fa -evalue 1e-10 -num_descriptions 10 -num_threads $NCPU \
    -out modelerunknown_blast_results.txt
if [ "$?" != "0" ]; then
    echo "ERROR: blastx exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -s modelerunknown_blast_results.txt ]; then
    echo "ERROR: zero-length output file from blastx: modelerunknown_blast_results.txt" >&2
    exit 1
fi
$PERL $DIR_CRL/transposon_blast_parse.pl --blastx modelerunknown_blast_results.txt --modelerunknown repeatmodeler_unknowns.fasta
if [ "$?" != "0" ]; then
    echo "ERROR: transposon_blast_parse.pl exited non-zero[$?]" >&2
    exit 1
fi
if [ ! -s unknown_elements.txt ]; then
    echo "ERROR: zero-length output file from transposon_blast_parse.pl" >&2
    exit 1
fi

mv unknown_elements.txt ModelerUnknown.lib
cat identified_elements.txt repeatmodeler_identities.fasta > ModelerID.lib

##
# Exclusion of gene fragments
if [ ! -e $DATADIR/uniprot_sprot.fasta.pin ]; then
    $DIR_BLAST/makeblastdb -in $DATADIR/uniprot_sprot.fasta -type prot
    if [ "$?" != "0" ]; then
        echo "ERROR: makeblastdb exited non-zero[$?]" >&2
        exit 1
    fi
fi

for LIB in ModelerID.lib ModelerUnknown.lib allMITE_LTR.lib; do
    $DIR_BLAST/blastx -query $LIB -db $DATADIR/uniprot_sprot.fasta -evalue 1e-10 -num_descriptions 10 -num_threads $NCPU \
        -out ${LIB}_blast_results.txt
    if [ "$?" != "0" ]; then
        echo "ERROR: blastx of $LIB exited non-zero[$?]" >&2
        exit 1
    fi
    if [ ! -s ${LIB}_blast_results.txt ]; then
        echo "ERROR: zero-length output file from blastx: ${LIB}_blast_results.txt" >&2
        exit 1
    fi

    $PERL $DIR_PE/ProtExcluder.pl ${LIB}_blast_results.txt ${LIB}
    if [ "$?" != "0" ]; then
        echo "ERROR: ProtExcluder.pl exited non-zero[$?]" >&2
        exit 1
    fi
    if [ ! -e ${LIB}noProtFinal ]; then
        echo "ERROR: ProtExcluder.pl output file ${LIB}noProtFinal not found" >&2
        exit 1
    fi
done
## allMITE_LTR classify

$PERL $DIR_RD/RepeatClassifier -pa $NCPU  -consensi allMITE_LTR.libnoProtFinal
    if [ ! -e  allMITE_LTR.libnoProtFinal.classified ]; then
        echo "ERROR: RepeatClassifier output file not found" >&2
        exit 1
    fi
cat allMITE_LTR.libnoProtFinal.classified  ModelerUnknown.libnoProtFinal ModelerID.libnoProtFinal >> allRepeats.lib

echo "Run complete"
