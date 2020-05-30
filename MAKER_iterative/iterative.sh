#!/bin/bash


module load mpich/ge/gcc/64/3.3

MPINCPU=68
NCPU=60

##

species=Bcarinata
seqfile=/data/gpfs/home/wyim/scratch/data/Brassica/BCEM_maker_after_Gap_pilon/BCEM_pilon.fasta


#SNAP=
#AUGUSTUS=
maker=~/scratch/bin/maker_mpich3/maker/bin/maker
gff3_merge=~/scratch/bin/maker_mpich3/maker/bin/gff3_merge
fasta_merge=~/scratch/bin/maker_mpich3/maker/bin/fasta_merge
maker2zff=~/scratch/bin/maker_mpich3/maker/bin/maker2zff
fathom=~/scratch/bin/maker/exe/snap/fathom
forge=~/scratch/bin/maker/exe/snap/forge
hmmassembler=/data/gpfs/home/wyim/scratch/bin/maker/exe/snap/hmm-assembler.pl
bedtools=~/scratch/bin/bedtools2/bin/bedtools
BUSCO=~/scratch/bin/busco-master/scripts/run_BUSCO.py
BUSCO_lineage=~/scratch/bin/busco-master/db/viridiplantae_odb10/
AUGUSTUS_species=$(echo $species)
AUGUSTUS_CONFIG_PATH=/data/gpfs/assoc/pgl/bin/augustus-3.3.2/config/


if [ ! -e $seqfile ]; then
    echo "ERROR: Sequence file $seqfile not found" >&2
    exit 1
fi

#if [ ! -e $AUGUSTUS_CONFIG_PATH/species/$AUGUSTUS_species ]; then
#    echo "ERROR: Sequence file $AUGUSTUS_CONFIG_PATH/species/$AUGUSTUS_species not found" >&2
#    exit 1
#fi


#sbatch -n 128 -J BCH709 -A cpu-s1-pgl-0 -p cpu-s1-pgl-0  --cpus-per-task 1 --mem-per-cpu=3800M -e maker_mpi_1a.e -o maker_mpi_1a.o --wrap="mpiexec -n 128 ~/scratch/bin/maker_mpich3/maker/bin/maker -base carinata_1 -fix_nucleotides -genome BcaBC.fasta"


for ((i=4;i<=6;i+=1))
do

echo $i


MAKER_JID=$(sbatch -n $MPINCPU --parsable -J MAKER -A cpu-s1-pgl-0 -p cpu-s1-pgl-0 -e maker_mpi_${i}.e -o maker_mpi_${i}a.o --cpus-per-task 1 --mem-per-cpu=4800M --wrap="/cm/shared/apps/mpich/ge/gcc/64/3.3/bin/mpiexec -n $MPINCPU  $maker -base  ${species}_${i} -fix_nucleotides -genome $seqfile round${i}_maker_opts.ctl maker_bopts.ctl maker_exe.ctl" )
#MAKER_JID=$(sbatch -n $MPINCPU --parsable -J MAKER -A cpu-s2-bch709-0 -p cpu-s2-core-0 -e maker_mpi_${i}.e -o maker_mpi_${i}a.o --cpus-per-task 1 --mem-per-cpu=3800M --wrap="/cm/shared/apps/mpich/ge/gcc/64/3.3/bin/mpiexec -n $MPINCPU  $maker -base  ${species}_${i} -fix_nucleotides -genome $seqfile round${i}_maker_opts.ctl maker_bopts.ctl maker_exe.ctl" )
echo $MAKER_JID

srun -A cpu-s1-pgl-0 -p cpu-s1-pgl-0 -c 1 --mem-per-cpu=400m  --dependency=afterok:${MAKER_JID} echo done




${gff3_merge} -d ${species}_${i}.maker.output/${species}_${i}_master_datastore_index.log
${gff3_merge} -n -s -d ${species}_${i}.maker.output/${species}_${i}_master_datastore_index.log > ${species}_${i}.all.maker_nofasta.gff
${fasta_merge} -d ${species}_${i}.maker.output/${species}_${i}_master_datastore_index.log

## Run SNAP training
mkdir -p snap/round${i} && cd snap/round${i}
${maker2zff} -x 0.25 -l 50 ../../${species}_${i}.all.gff
${fathom} genome.ann genome.dna -gene-stats > gene-stats.log 2>&1
${fathom} genome.ann genome.dna -validate > validate.log 2>&1
${fathom} genome.ann genome.dna -categorize 1000 > categorize.log 2>&1
${fathom} uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1

${forge} export.ann export.dna > forge.log 2>&1
#### assembly the HMM
${hmmassembler} ${species}_round${i} . > ${species}_round${i}.hmm
cd ../../




### AUGUSTUS
mkdir -p augustus/round${i} && cd augustus/round${i}

awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' ../../${species}_${i}.all.maker_nofasta.gff | awk -v OFS="\t" '{ if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000 }' | bedtools getfasta -fi $seqfile -bed - -fo  ${species}_${i}.transcripts1000.fasta

### CONDA
source ~/scratch/bin/miniconda3/etc/profile.d/conda.sh
conda activate busco

### BUSCO

AUGUSTUS_CONFIG_PATH=/data/gpfs/assoc/pgl/bin/augustus-3.3.2/config/

BUSCOJID=$(sbatch -n $NCPU --parsable -J BUSCO -A cpu-s1-pgl-0 -p cpu-s1-pgl-0 -e BUSCO_${i}.e -o BUSCO_${i}a.o  --wrap="$BUSCO -i ${species}_${i}.transcripts1000.fasta  -o  ${species}_${i}_maker -l $BUSCO_lineage -m genome -c $NCPU --long -sp carinata1  -z --augustus_parameters='--progress=true'")

echo $BUSCOJID

                      srun -A cpu-s1-pgl-0 -p cpu-s1-pgl-0 -c 1 --mem-per-cpu=400m  --dependency=afterok:${BUSCOJID} echo done


cd run_${species}_${i}_maker/augustus_output/retraining_parameters

PID=$( ls -1  -d BUSCO*.txt | sed 's/.*maker//g; s/[^0-9]*//g')

rename "s/BUSCO_${species}_${i}_maker_${PID}/${species}${i}/g" *

sed -i "s/BUSCO_${species}_${i}_maker_${PID}/${species}${i}/g" ${species}${i}_parameters.cfg
sed -i "s/BUSCO_${species}_${i}_maker_${PID}/${species}${i}/g" ${species}${i}_parameters.cfg.orig1

mkdir -p $AUGUSTUS_CONFIG_PATH/species/${species}${i}
cp  *  $AUGUSTUS_CONFIG_PATH/species/${species}${i}/
conda deactivate
cd ../../../../../




# transcript alignments
#awk '{ if ($2 == "est2genome") print $0 }'  ${species}_1.all.maker_nofasta.gff > ${species}${i}.all.maker.est2genome.gff
# protein alignments
#awk '{ if ($2 == "protein2genome") print $0 }' ${species}_1.all.maker_nofasta.gff > ${species}${i}.all.maker.protein2genome.gff

cp ${species}1.all.maker.est2genome.gff ${species}${i}.all.maker.est2genome.gff
cp  ${species}1.all.maker.protein2genome.gff  ${species}${i}.all.maker.protein2genome.gff
cp ${species}1.all.maker.repeats.gff ${species}${i}.all.maker.repeats.gff

#awk '{ if ($2 ~ "repeat") print $0 }' ${species}_1.all.maker_nofasta.gff  | grep -v -e "Satellite" -e "Low_complexity" -e "Simple_repeat" | perl -ane '$id; if(!/^\#/){@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\t", @F)."\n"} print $_' \ >   ${species}${i}.all.maker.repeats.gff


### CTL file
k=$((i+1))
cp round${i}_maker_opts.ctl round${k}_maker_opts.ctl

sed -i "/^snaphmm=/ s/snaphmm=.*/snaphmm=snap\/round${i}\/${species}_round${i}.hmm/" round${k}_maker_opts.ctl
sed -i "/^augustus_species=/ s/augustus_species=.*/augustus_species=${species}${i}/" round${k}_maker_opts.ctl
sed -i "/^protein_gff=/ s/protein_gff=.*/protein_gff=${species}${i}\.all\.maker\.protein2genome\.gff/" round${k}_maker_opts.ctl
sed -i "/^est_gff=/ s/est_gff=.*/est_gff=${species}${i}\.all\.maker\.est2genome\.gff/" round${k}_maker_opts.ctl
#sed -i "/^rm_gff=/ s/rm_gff=.*/rm_gff=${species}${i}\.all\.maker\.repeats\.gff/" round${k}_maker_opts.ctl

#sed -i "/^maker_gff=/ s/maker_gff=.*/maker_gff=${species}_1.all.maker_nofasta.gff/" round${k}_maker_opts.ctl

#sed -i "/^rmlib=/ s/rmlib=.*/rmlib=/" round${k}_maker_opts.ctl

#sed -i "/^model_org=/ s/model_org=all/model_org=simple/" round${k}_maker_opts.ctl

sed -i "/^est=/ s/est=.*/est=/" round${k}_maker_opts.ctl
sed -i "/^protein=/ s/protein=.*/protein=/" round${k}_maker_opts.ctl


sed -i "/^est_pass=0/ s/est_pass=0/est_pass=1/" round${k}_maker_opts.ctl
sed -i "/^protein_pass=0/ s/protein_pass=0/protein_pass=1/" round${k}_maker_opts.ctl
sed -i "/^protein_pass=0/ s/protein_pass=0/protein_pass=1/" round${k}_maker_opts.ctl

sed -i "/^gmhmm=/ s/gmhmm=.*/gmhmm=/" round${k}_maker_opts.ctl

sed -i "/^est2genome=/ s/est2genome=1/est2genome=0/" round${k}_maker_opts.ctl
sed -i "/^protein2genome=/ s/protein2genome=1/protein2genome=0/" round${k}_maker_opts.ctl



done
