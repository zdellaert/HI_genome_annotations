#!/usr/bin/env bash
#SBATCH --export=NONE
#SBATCH --nodes=1 --ntasks-per-node=24
#SBATCH --signal=2
#SBATCH --no-requeue
#SBATCH --mem=32GB
#SBATCH --error="%x_error.%j"
#SBATCH --output="%x_output.%j"
#SBATCH -t 01:00:00
#SBATCH --mail-type=BEGIN,END,FAIL #email you when job starts, stops and/or fails

SINGULARITY_IMAGE="docker://memesuite/memesuite:latest"

cd ../annotation/promoters
mkdir -p fimo_output

module load apptainer/latest

species=(Mcap Pacuta Pcomp)
fastas=(
  "Montipora_capitata_HIv3_promoters_500_upstream.fasta"
  "Pocillopora_acuta_HIv2_promoters_500_upstream.fasta"
  "Porites_compressa_HIv1_promoters_500_upstream.fasta"
)

# run FIMO 

for i in "${!species[@]}"; do
    sp="${species[$i]}"
    fasta="${fastas[$i]}"

    echo "Starting analysis of $sp using $fasta"
    
    # run FIMO with default settings (--thresh 0.0001, --max-stored-scores 100000)
    singularity exec --cleanenv $SINGULARITY_IMAGE fimo \
        -oc "fimo_output/${sp}_stress_TFs" \
        --thresh 0.0001 \
        --max-stored-scores 100000 \
        ../../references/motif_databases/stress_TFs.meme \
        "$fasta"
done