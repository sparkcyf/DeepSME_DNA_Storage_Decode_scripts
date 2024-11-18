# clean old files
bash clean.sh

# set path of fastq files (you need to copy fastq to "fastq" folder)
export FASTQ_FOLDER="./fastq/" 
# set path for muscle
export MUSCLE_PATH="./bin/muscle5.1.linux_intel64"
# set path for minimap2
export MINIMAP2_PATH="./bin/minimap2-2.28_x64-linux.tar.bz2"
python3 mapping.py
python3 paf_to_tsv.py
python3 readsdic_gen.py
python3 seq_grouping.py

SEQ_BATCH_SIZES=(16 32 64 128)
RANDOM_SEEDS=(0 1 2 3 4)
OUTPUT_FASTA_FN="deepsme_verify"


# assembly
mkdir -p assembly_log


for SEQ_BATCH_SIZE in "${SEQ_BATCH_SIZES[@]}"; do
  for RANDOM_SEED in "${RANDOM_SEEDS[@]}"; do
    # log fn
    LOG_FILENAME="${OUTPUT_FASTA_FN}-bs${SEQ_BATCH_SIZE}-rs${RANDOM_SEED}.log"

    export SEQ_BATCH_SIZE
    export RANDOM_SEED
    export OUTPUT_FASTA_FN
    
    echo "Running with SEQ_BATCH_SIZE=$SEQ_BATCH_SIZE, RANDOM_SEED=$RANDOM_SEED"

    echo "Running with SEQ_BATCH_SIZE=$SEQ_BATCH_SIZE, RANDOM_SEED=$RANDOM_SEED" > assembly_log/$LOG_FILENAME
    
    python3 -u alignmentNassembly.py >> assembly_log/$LOG_FILENAME 2>&1
  done
done

echo "All jobs completed."

# decode
DECODE_RES_DIR="decode_res/"
INPUT_FASTA_DIR="assembly_results/"

mkdir -p $DECODE_RES_DIR

for SEQ_BATCH_SIZE in "${SEQ_BATCH_SIZES[@]}"; do
  for RANDOM_SEED in "${RANDOM_SEEDS[@]}"; do

    LOG_FILENAME="${OUTPUT_FASTA_FN}-bs${SEQ_BATCH_SIZE}-rs${RANDOM_SEED}.log"
    FAILED_INDEX_TXT="${DECODE_RES_DIR}/${OUTPUT_FASTA_FN}-bs${SEQ_BATCH_SIZE}-rs${RANDOM_SEED}.txt"
    INPUT_FASTA_PATH="${INPUT_FASTA_DIR}/${OUTPUT_FASTA_FN}-bs${SEQ_BATCH_SIZE}-rs${RANDOM_SEED}.fasta"
    DECODE_LOG="${DECODE_RES_DIR}/${OUTPUT_FASTA_FN}-bs${SEQ_BATCH_SIZE}-rs${RANDOM_SEED}.log"
    
    echo "Running decoding with SEQ_BATCH_SIZE=$SEQ_BATCH_SIZE, RANDOM_SEED=$RANDOM_SEED"

    FAILED_INDEX_TXT=$FAILED_INDEX_TXT INPUT_FASTA_PATH=$INPUT_FASTA_PATH python3 -u sustech_decode.py >> $DECODE_LOG 2>&1
  done
done

echo "All decoding jobs completed."
