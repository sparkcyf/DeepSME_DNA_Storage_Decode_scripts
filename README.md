# Decode Scripts for 5hmc-modified DNA Storage Datasets Basecalled by DeepSME

## Prerequisites

1. Install conda environment for decoding (derived from [Composite-Hedges-Nanopores](https://github.com/ysfhtxn/Composite-Hedges-Nanopores) )

```bash
conda env create -f environment.yml && conda activate CHN

# optional, if you want to use mpi4py for parallel decoding
# you need to have MPI installed on your system first
pip install mpi4py
```

2. Install minimap2 and [muscle5](https://github.com/rcedgar/muscle) for multiple sequence alignment

3. Set the environment variable `MINIMAP2_PATH`, `MUSCLE_PATH` in `run.sh` to the path of minimap2 and muscle5.

4. Create a folder named `fastq` and put the basecalled fastq files in it. (you can download the sample fastq from NCBI SRA: https://www.ncbi.nlm.nih.gov/sra/SRX28067610 )

```bash
mkdir fastq
cp /path/to/fastq/xxx.fastq fastq/basecalling.fastq
```

## Usage

> [!NOTE]
> Expected run time for demo decoding script on a server with AMD Ryzen 9 5950X, 128G of RAM and SSD storage is about **160** minutes.

```bash
# run the decoding script
time bash run.sh

# if you want to use mpi4py for parallel decoding
time bash run_mpi.sh
```

## Output

After the decode, you can check the txt files output in the `decode_res` folder. These files contain the decoded accuracy at each sequence batch size (like the example below). We have put some logs in the `decode_res` folder for reference.

```
Decoding Summary:
Total segments processed: 487
Failed segments: 9 (1.85%)
Total bytes processed during decoding: 8753
Total mismatched bytes during decoding: 162/8753 (1.85%)

File Comparison Results:
Picture accuracy: 7629/7775 (98.12%)
Text accuracy: 978/978 (100.00%)
Weighted average accuracy: 98.33%
Failed indices: [106, 220, 227, 235, 304, 311, 333, 391, 414]
decode finished!
```

> [!NOTE]
> We treated the text and image as binary data and calculated the accuracy by comparing the decoded binary with the ground truth.


