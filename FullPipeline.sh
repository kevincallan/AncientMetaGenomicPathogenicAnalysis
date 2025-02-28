#!/bin/bash

# Automated M. tuberculosis Analysis Pipeline
# This script automates the download and analysis of both ancient and modern M. tuberculosis samples

# Default values
TAXID=1773  # M. tuberculosis taxonomy ID
NUM_GENOMES=50
THREADS=20
OUTPUT_DIR="mtb_analysis"
ANCIENT_SAMPLES=(
    "ERR650569" "ERR650974" "ERR650975" "ERR650978" "ERR650979"
    "ERR650981" "ERR650982" "ERR650983" "ERR651000" "ERR651001"
    "ERR651002" "ERR651003" "ERR651004"
)

# Help function
usage() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -t, --taxid INT        Taxonomy ID for modern samples (default: 1773)"
    echo "  -n, --num-genomes INT  Number of modern genomes to download (default: 50)"
    echo "  -p, --threads INT      Number of threads to use (default: 20)"
    echo "  -o, --output DIR       Output directory (default: mtb_analysis)"
    echo "  -h, --help            Show this help message"
    exit 1
}

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -t|--taxid) TAXID="$2"; shift ;;
        -n|--num-genomes) NUM_GENOMES="$2"; shift ;;
        -p|--threads) THREADS="$2"; shift ;;
        -o|--output) OUTPUT_DIR="$2"; shift ;;
        -h|--help) usage ;;
        *) echo "Unknown parameter: $1"; usage ;;
    esac
    shift
done

# Create directory structure
mkdir -p "${OUTPUT_DIR}"/{ancient,modern,results}/{raw,processed,aligned,stats}
mkdir -p "${OUTPUT_DIR}/logs"

# Log function
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "${OUTPUT_DIR}/logs/pipeline.log"
}

# Error handling
set -e
trap 'log "Error occurred on line $LINENO. Exit code: $?"' ERR

# Function to download ancient samples
download_ancient_samples() {
    log "Starting ancient sample downloads..."
    
    for sample in "${ANCIENT_SAMPLES[@]}"; do
        log "Downloading ${sample}..."
        for read in 1 2; do
            wget -P "${OUTPUT_DIR}/ancient/raw" \
                "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${sample:0:6}/${sample}/${sample}_${read}.fastq.gz"
        done
    done
}

# Function to download modern samples
download_modern_samples() {
    log "Starting modern sample downloads..."
    
    # Create temporary directory for esearch results
    mkdir -p "${OUTPUT_DIR}/tmp"
    
    # Search for and download modern samples
    esearch -db assembly -query "txid${TAXID}[Organism] AND (latest[filter] AND (complete genome[filter] OR chromosome level[filter]))" \
        | esummary \
        | xtract -pattern DocumentSummary -element FtpPath_RefSeq \
        > "${OUTPUT_DIR}/tmp/ftppaths.txt"
    
    # Process each FTP path
    count=0
    while read -r url && [ $count -lt $NUM_GENOMES ]; do
        if [ -z "$url" ]; then continue; fi
        
        fname=$(basename "$url")
        log "Downloading ${fname}..."
        
        wget -P "${OUTPUT_DIR}/modern/raw" "${url}/${fname}_genomic.fna.gz"
        wget -P "${OUTPUT_DIR}/modern/raw" "${url}/${fname}_genomic.gbff.gz"
        
        count=$((count + 1))
    done < "${OUTPUT_DIR}/tmp/ftppaths.txt"
}

# Function to process ancient samples
process_ancient_samples() {
    log "Processing ancient samples..."
    
    # Run adapter removal in parallel
    for sample in "${ANCIENT_SAMPLES[@]}"; do
        log "Running adapter removal for ${sample}..."
        AdapterRemoval \
            --threads $THREADS \
            --file1 "${OUTPUT_DIR}/ancient/raw/${sample}_1.fastq.gz" \
            --file2 "${OUTPUT_DIR}/ancient/raw/${sample}_2.fastq.gz" \
            --collapse \
            --minadapteroverlap 1 \
            --minlength 25 \
            --minquality 25 \
            --gzip \
            --trimns \
            --trimqualities \
            --basename "${OUTPUT_DIR}/ancient/processed/${sample}.trimmed25bp25q" &
    done
    wait
    
    # Move collapsed files
    mv "${OUTPUT_DIR}"/ancient/processed/*.collapsed.gz "${OUTPUT_DIR}/ancient/processed/"
}

# Function to align samples
align_samples() {
    log "Starting alignment process..."
    
    # Index reference if needed
    if [ ! -f "${OUTPUT_DIR}/reference/ref.fasta.bwt" ]; then
        bwa index "${OUTPUT_DIR}/reference/ref.fasta"
    fi
    
    # Align ancient samples in parallel
    for sample in "${ANCIENT_SAMPLES[@]}"; do
        log "Aligning ${sample}..."
        bwa mem -t $THREADS \
            "${OUTPUT_DIR}/reference/ref.fasta" \
            "${OUTPUT_DIR}/ancient/processed/${sample}.trimmed25bp25q.collapsed.gz" \
            > "${OUTPUT_DIR}/ancient/aligned/${sample}.sam" &
    done
    wait
    
    # Convert SAM to BAM and sort
    for sample in "${ANCIENT_SAMPLES[@]}"; do
        samtools view -bS "${OUTPUT_DIR}/ancient/aligned/${sample}.sam" \
            | samtools sort -@ $THREADS -o "${OUTPUT_DIR}/ancient/aligned/${sample}.sorted.bam"
        samtools index "${OUTPUT_DIR}/ancient/aligned/${sample}.sorted.bam"
    done
}

# Function to run quality control
run_qc() {
    log "Running quality control..."
    
    # Run PRINSEQ++ on ancient samples
    for sample in "${ANCIENT_SAMPLES[@]}"; do
        prinseq++ \
            -min_len 25 \
            -min_qual_mean 25 \
            -ns_max_n 1 \
            -min_gc 20 \
            -max_gc 80 \
            -input "${OUTPUT_DIR}/ancient/aligned/${sample}.sorted.bam" \
            -out_name "${OUTPUT_DIR}/ancient/processed/${sample}.qc"
    done
    
    # Generate QC reports
    multiqc "${OUTPUT_DIR}/ancient/processed" -o "${OUTPUT_DIR}/results/qc_report"
}

# Function to run SNP analysis
run_snp_analysis() {
    log "Running SNP analysis..."
    
    # Run snippy on each sample
    for sample in "${ANCIENT_SAMPLES[@]}"; do
        snippy --cpus $THREADS \
            --outdir "${OUTPUT_DIR}/results/snps/${sample}" \
            --ref "${OUTPUT_DIR}/reference/ref.fasta" \
            --se "${OUTPUT_DIR}/ancient/processed/${sample}.qc.fastq"
    done
    
    # Run snippy-core
    snippy-core --ref "${OUTPUT_DIR}/reference/ref.fasta" "${OUTPUT_DIR}/results/snps/*"
}

# Function to build phylogenetic tree
build_tree() {
    log "Building phylogenetic tree..."
    
    # Run IQ-TREE
    iqtree -s "${OUTPUT_DIR}/results/snps/core.aln" \
        -m MFP \
        -bb 1000 \
        -nt AUTO \
        -pre "${OUTPUT_DIR}/results/tree"
}

# Main pipeline execution
main() {
    log "Starting MTB analysis pipeline..."
    
    # Download samples
    download_ancient_samples
    download_modern_samples
    
    # Process ancient samples
    process_ancient_samples
    
    # Align samples
    align_samples
    
    # Run QC
    run_qc
    
    # Run SNP analysis
    run_snp_analysis
    
    # Build phylogenetic tree
    build_tree
    
    log "Pipeline completed successfully!"
}

# Execute main pipeline
main
