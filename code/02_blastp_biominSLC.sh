#!/bin/bash
#SBATCH -t 18:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=500GB
#SBATCH --export=NONE
#SBATCH --error=outs_errs/"%x_error.%j" #write out slurm error reports
#SBATCH --output=outs_errs/"%x_output.%j" #write out any program outpus
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails

module load blast-plus/2.14.1

cd ../references/blast_dbs/

for database in *_prot.pin; do
    species=$(basename "$database" _prot.pin)

    blastp \
    -query "../biomineralization/SLC4.fasta" \
    -db ${species}_prot \
    -out "../biomineralization/${species}_SLC4_BLAST_out.tab" \
    -evalue 1E-05 \
    -num_threads 48 \
    -max_target_seqs 5 \
    -max_hsps 1 \
    -outfmt 6
done

for database in *_prot.pin; do
    species=$(basename "$database" _prot.pin)

    blastp \
    -query "../biomineralization/SLC4.fasta" \
    -db ${species}_prot \
    -out "../biomineralization/${species}_SLC4_BLAST.tab" \
    -evalue 1E-05 \
    -num_threads 48 \
    -max_target_seqs 5 
done

echo "Blast complete!" $(date)
