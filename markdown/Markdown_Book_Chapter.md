# Book Chapter: Computational workflow for peak calling from scCUT&Tag

> by Md Wahiduzzaman
>
> wahid810@gmail.com; wahid@hnu.edu.cn
>
> School of life sciences, Hunan University
>
> Changsha-410082, China.

## 1 Introduction

&#x20;Single‑cell CUT&Tag (scCUT&Tag) enables genome‑wide mapping of protein‑DNA interactions in individual cells by combining antibody‑targeted tagmentation with single‑cell barcoding [1,2]. Unlike single‑cell ATAC‑seq, which captures open chromatin regions independent of a specific antibody, scCUT&Tag provides targeted, low‑background profiles of specific histone modifications (e.g., H3K27ac, H3K27me3) or transcription factors (e.g., Rad21 and Olig2), offering a versatile approach for investigating regulatory mechanisms across diverse biological systems.

To facilitate consistent and reproducible analysis of scCUT&Tag data across different HMs and TFs, this chapter provides a comprehensive step‑by‑step computational workflow. We present a modular pipeline applicable to widely used scCUT&Tag designs, including human PBMC (hPBMC) and mouse brain (mBrain) datasets. The core analytical workflow comprises four main stages: (i) initial data processing using cellranger‑atac, (ii) generation of input control BAM files, (iii) splitting BAM files by cell barcodes, and (iv) peak calling.

The protocol begins with raw FASTQ files containing cell barcodes, assuming that cell demultiplexing has already been performed by the sequencing facility or via a standard barcode extraction script. Reads are preprocessed and aligned per cell. An input control BAM file is generated for each cell‑type using the remaining cell‑types. Peak calling is then performed using seven established peak callers including DROMPA+ [3], Genrich [4], GoPeaks [5], HOMER [6], MACS2 [7], SEACR [8], SICER2 [9] those are allowing flexibility to compare both input‑guided and input‑free approaches. Count, for each cell, the number of fragments overlapping each peak. Following peak identification, we count the read widths number of each peak per cell across the HMs and TFs.

## 2 Materials

### 2.1 Software requirements

The following software is necessary to run the peak calling workflow introduced in this protocol. The versions shown in parentheses were used during the preparation of this protocol and have been tested accordingly.

·        Bash (version 5.2.2): A Unix shell and command language (https://www.gnu.org/software/bash/).

·        Cell Ranger ATAC (version 2.1.0): A set of analysis pipelines from 10x Genomics for processing and analyzing single-cell Epi ATAC (formerly Single Cell ATAC) data, including read alignment, peak calling, cell calling, and count matrix generation (https://www.10xgenomics.com/support/software/) [10].

·        fastq-dump (version 2.11.3): A tool from the NCBI SRA Toolkit for converting SRA (Sequence Read Archive) files into FASTQ format (https://github.com/ncbi/sra-tools).

·        Samtools (Version 1.13): A suite of programs for processing and analyzing high-throughput sequencing data, including tools for file format conversion, sorting, querying, statistics, and variant calling (https://www.htslib.org/) [11].

·        DROMPA+ (version 1.20.1): A user-friendly peak-calling and visualization tool for multiple ChIP‑seq datasets, designed for quality control, PCR bias filtering, normalization, and differential analysis (nakatolab.iqb.u-tokyo.ac.jp/softwares/)[3].

·        Genrich (version 0.6.1): A peak-caller for genomic enrichment assays such as ChIP‑seq and ATAC‑seq that analyzes alignment files and produces peaks of significant enrichment (https://github.com/jsh58/Genrich) [4].  

·        GoPeaks (version 1.0.0): A peak caller specifically designed for CUT&Tag/CUT&RUN sequencing data, working best with narrow peaks such as H3K4me3 and transcription factors (https://git-lfs.github.com) [5].

·        HOMER (version 5.1): A comprehensive suite for motif discovery and ChIP‑seq analysis, containing the findPeaks program for peak calling and transcript identification (http://homer.ucsd.edu/homer/) [6].

·        MACS2 (version 2.2.9.1): Model‑based Analysis of ChIP‑Seq (MACS) for identifying transcription factor binding sites and enriched regions in ChIP‑seq data, applicable to any DNA enrichment assay (https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html) [7].

·        SEACR (version 1.3): Sparse Enrichment Analysis for CUT&RUN, intended to call peaks and enriched regions from sparse chromatin profiling data where background is dominated by regions with no read coverage (https://github.com/FredHutch/SEACR) [8].

·        SICER2 (version 1.0.3): A redesigned and improved ChIP‑seq broad peak calling tool, specifically for identifying broad histone modification domains (https://github.com/zanglab/SICER2) [9].

·        R (version 4.3.3): A free software environment for statistical computing and graphics. It compiles and runs on a wide variety of UNIX platforms, Windows, and macOS (https://www.r-project.org/).For installation examples refer to **Note 1**.

·        Seurat (version 5.3.0): A R package designed for QC, analysis, and exploration of single-cell RNA-seq data(https://satijalab.org/seurat/) [12].

·        Signac (version 1.15.0): A comprehensive R package for the analysis of single-cell chromatin data that designed for scATAC-seq, single-cell targeted tagmentation methods such as scCUT&Tag, and multimodal datasets (https://stuartlab.org/signac/) [13].

### 2.2 Data preprocessing

The raw sequencing data were obtained from the Sequence Read Archive (SRA) under accession numbers for GSE195725 [14] and GSE157637 [1]. To facilitate downstream processing with Cell Ranger ATAC, the FASTQ files were renamed according to the 10x Genomics naming convention. This was achieved by creating symbolic links from the original SRA-derived files to the required format, which includes sample name, lane, and read type identifiers (I1, R1, I2, R2). The symbolic links were generated to avoid duplicating the large raw data files while ensuring compatibility with the Cell Ranger ATAC pipeline. The cellranger-atac count pipeline was executed with the hg38 (for hPBMC) and  mm10 (for mBrain) reference genome (refdata-cellranger-arc-mm10-2020-A-2.0.0, refdata-cellranger-arc-hg38-2020-A-2.0.0) obtained from 10x Genomics (see **Note 2** ). The pipeline performed read alignment using BWA, barcode processing, cell calling, and transposase‑accessible region identification. The run was configured with 64 cores and 128 GB of memory to ensure efficient processing of the large dataset. The main output was the possorted_bam.bam file, which contains all aligned, cell‑barcoded reads. This BAM file served as the real BAM files.



The scCUT&Tag data preprocessing workflow using Signac [14] begins with loading fragment files into a ChromatinAssay object, followed by quality control wherein cells are filtered based on total fragment counts and nucleosome signal patterns to remove empty droplets, doublets, and technical artifacts; counts are then quantified over fixed genomic bins and normalized using the term frequency–inverse document frequency (TF‑IDF) transformation, which accounts for differences in sequencing depth across cells and emphasizes features with specific biological relevance. Subsequently, linear dimensional reduction is performed via latent semantic indexing (LSI) on the normalized matrix, and the first LSI dimension—typically dominated by technical noise—is discarded while the remaining top dimensions are retained for downstream analysis; these reduced-dimension components serve as input for uniform manifold approximation and projection (UMAP), yielding a two‑dimensional embedding that reveals the underlying cellular heterogeneity. The resulting clusters are annotated by integrating with reference chromatin profiles or matched single‑cell RNA‑seq data via the Weighted Nearest Neighbor framework, which presents the UMAP visualization with cells colored by their assigned identities alongside supplementary quality‑control panels. Figure 1 showing percentage cell and read numbers on of cell-type -specific across the HMs/TFs (**Note 3** ).

### 2.3 BAM file availability

The H3K4me3 hPBMC BAM file is accessible on Zenodo (https://zenodo.org/records/20797529). For the remaining datasets, we retrieved the raw FASTQ files and associated processed data from the Gene Expression Omnibus under accession numbers GSE195725 [14]  (hPBMC) and GSE157637 [1] (mBrain). Using the cell barcode information derived from these processed data, we generated simulated BAM files across all cell types (**Section 3.3**). Furthermore, we extracted cell-type-specific BAM files directly from the real BAM alignments, following the splitting strategy outlined in **Section 3.1**.      

## 3. Methods

### 3.1 Split the BAM file by cell type

To generate cell-type-specific BAM files from scCUT&Tag datasets, we first extracted the cell barcodes associated with each annotated cell type from the processed Seurat object. These barcodes were subsequently used to split the original, unfiltered BAM file, yielding individual BAM files for each cell type.

&#x20;

An example `Bash`workflow for processing the H3K4me3 dataset from human PBMCs (hPBMC) is provided below.



**NOTE:** The following `Bash`code only for manuscript (**H3K4me3 for hPBMC**)&#x20;

````r
```bash
########################################################
##### Splitting bam for H3K4me3 of hPBMC
########################################################

cd ~/project_scHMTF/GSE195725_processed_data/splitbam_realbam/H3K4me3/
#!/bin/bash

INPUT_BAM="H3K4me3.bam"
OUTPUT_DIR="split_celltype_bams_corrected"

mkdir -p "$OUTPUT_DIR"

echo "=== Starting Cell Type BAM Splitting ==="

for barcode_file in *barcodes.txt; do
    # Skip temp files
    if [[ "$barcode_file" == "temp_bam_barcodes.txt" ]]; then
        continue
    fi
    
    celltype=$(echo "$barcode_file" | sed 's/_barcodes.txt//')
    echo "Splitting $celltype..."
    
    # Use awk method - most reliable
    samtools view -h "$INPUT_BAM" | \
    awk -F'\t' -v barcode_file="$barcode_file" '
    BEGIN {
        # Load barcodes into array
        while((getline line < barcode_file) > 0) {
            gsub(/-.*$/, "", line)   # REMOVE -1/-2 suffix
            barcodes[line] = 1
        }
        close(barcode_file)
    }
    /^@/ { 
        # Print all header lines
        print 
        next 
    }
    {
        # Check for CB tag in optional fields
        for(i=12; i<=NF; i++) {
            if($i ~ /^CB:Z:/) {
                split($i, arr, ":")
                cb = arr[3]
                cb2 = cb
                sub(/-.*$/, "", cb2)
                if(cb2 in barcodes) {
                    print
                    break
                }
            }
        }
    }' | samtools view -b - > "${OUTPUT_DIR}/H3K4me3_${celltype}.bam"
    
    # Index the BAM file
    samtools index "${OUTPUT_DIR}/H3K4me3_${celltype}.bam"
    
    # Count reads
    read_count=$(samtools view -c "${OUTPUT_DIR}/H3K4me3_${celltype}.bam")
    echo "$celltype: $read_count reads"
done

echo "=== Splitting Complete ==="
```
````

This approach was applied systematically across all histone modifications and cell types. For the hPBMC dataset, we processed BAM files for H3K27ac, H3K27me3, H3K4me1, H3K4me2, and H3K9me3. For the mouse brain (mBrain) dataset, we similarly processed H3K27ac, H3K27me3, H3K36me3, H3K4me3, Olig2, and Rad21. The complete set of scripts used for splitting BAM files and subsequent analyses is available in our public GitHub repository:

https://github.com/wahidthedream/scCUT-Tag_protocol/





### 3.2 Call peaks from the real BAM files without an input controI using MACS2 and SICER2

&#x20;To identify enriched genomic regions, we performed peak calling on the cell‑type‑specific BAM files for both the hPBMC and mBrain datasets in the absence of an input control. We employed seven distinct peak callers: DROMPA+, Genrich, GoPeaks, HOMER, MACS2, SEACR, and SICER2. For the hPBMC dataset, we processed H3K27me3 and H3K9me3 as broad marks; H3K4me1, H3K4me2, and H3K4me3 as narrow (sharp) marks; and H3K27ac as both broad and narrow. For the mBrain dataset, we processed H3K27me3 and H3K36me3 as broad, and H3K4me3, Olig2, and Rad21 as narrow, with H3K27ac again run in both modes. To ensure comparability across methods, we standardized the key input parameters for each caller according to the peak type, as detailed below:

&#x20;

1. For narrow (sharp) peak calling: DROMPA+ (--pthre_internal 5); Genrich (-a 200, -l 100, -g 100, -p 0.01); GoPeaks (default parameters); HOMER (factor mode); MACS2 (-q 0.05, --keep-dup all); SEACR (relaxed, 0.01, non); and SICER2 (--window_size 50, --gap 100, -p 0.01).

2. For broad peak calling: DROMPA+ (--pthre_internal 4); Genrich (-a 100, -l 500, -g 1000, -p 0.05); GoPeaks (--broad); HOMER (histone mode); MACS2 (--broad, -q 0.05, --keep-dup all); SEACR (stringent, 0.01, non); and SICER2 (--window_size 100, --gap 200, -p 0.05).

 

Following peak calling with each tool, we generated and subsequently explored the resulting peak BED files for downstream analysis.

As an illustration, the following “bash” workflows demonstrate peak calling for the H3K4me3 hPBMC dataset using MACS2 and SICER2, respectively. To execute these scripts, a dedicated Conda environment containing all required dependencies must first be created and activated.



**NOTE:** The following `bash` code for all HMs of hPBMC by using MACS2.

````bash

```bash
#!/bin/bash
###############################################################################################################
### MACS2 Peak Calling with Control (Input BAMs)
### Multi-histone processing: Broad + Narrow marks
### Human PBMC real BAMs
###############################################################################################################

# Define histone marks
BROAD_MARKS=("H3K27ac" "H3K27me3" "H3K9me3")
NARROW_MARKS=("H3K27ac" "H3K4me1" "H3K4me2" "H3K4me3")

# Define cell types
CELLTYPES=("B" "CD4T" "CD8T" "DC" "Mono" "NK" "otherT" "other")

# Base directories
THREADS=32
BASE_DIR="/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam"
OUT_BASE="${BASE_DIR}/HumanPBMC_peakbed/MACS2_peakbed"
GENOME=2.7e9  # Human genome size for MACS2

mkdir -p "${OUT_BASE}"

# -------------------------------------------------------------------------------------------
# Function: Run MACS2 for all cell types under a given histone mark
# -------------------------------------------------------------------------------------------
run_macs2_all_ctrl() {
    local MARK=$1
    local MODE=$2
    local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams_corrected"
    local OUT_DIR="${OUT_BASE}/${MARK}_${MODE}/With_input_corrected"
    mkdir -p "${OUT_DIR}"

    echo "==================== Processing ${MARK} (${MODE}) ===================="

    for CELL in "${CELLTYPES[@]}"; do
        local TREAT="${MARK_DIR}/${MARK}_${CELL}.bam"
        local CTRL="${MARK_DIR}/input_${CELL}.bam"

        echo "========== Running MACS2 for: ${CELL} (${MARK}) =========="

        if [[ ! -f "$TREAT" || ! -f "$CTRL" ]]; then
            echo "Missing treatment or control BAM for ${CELL}. Skipping..."
            continue
        fi

        # Set parameters
        local QVAL=0.05
        local ARGS=" -q $QVAL --keep-dup all"

        if [[ "$MODE" == "broad" ]]; then
            ARGS="$ARGS --broad --broad-cutoff $QVAL"
        fi

        macs2 callpeak \
            -t "$TREAT" \
            -c "$CTRL" \
            -f BAMPE \
            -g "$GENOME" \
            -n "${MARK}_${CELL}" \
            --outdir "$OUT_DIR" \
            $ARGS

        echo "Finished: ${MARK}_${CELL}"
        echo
    done
}

# -------------------------------------------------------------------------------------------
# Run for both broad and narrow marks
# -------------------------------------------------------------------------------------------
for MARK in "${BROAD_MARKS[@]}"; do
    run_macs2_all_ctrl "$MARK" "broad"
done

for MARK in "${NARROW_MARKS[@]}"; do
    run_macs2_all_ctrl "$MARK" "narrow"
done


#!/bin/bash
###############################################################################################################
### MACS2 Peak Calling (Without Control)
### Multi-histone processing: Broad + Narrow marks
### Human PBMC real BAMs
###############################################################################################################

# Define histone marks
BROAD_MARKS=("H3K27ac" "H3K27me3" "H3K9me3")
NARROW_MARKS=("H3K27ac" "H3K4me1" "H3K4me2" "H3K4me3")

# Define cell types
CELLTYPES=("B" "CD4T" "CD8T" "DC" "Mono" "NK" "otherT" "other")

# Base directories
BASE_DIR="/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam"
OUT_BASE="${BASE_DIR}/HumanPBMC_peakbed/MACS2_peakbed"
GENOME=2.7e9  # Human genome size for MACS2

mkdir -p "${OUT_BASE}"

# -------------------------------------------------------------------------------------------
# Function: Run MACS2 for all cell types under a given histone mark (no control)
# -------------------------------------------------------------------------------------------
run_macs2_no_ctrl() {
    local MARK=$1
    local MODE=$2
    local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams_corrected"
    local OUT_DIR="${OUT_BASE}/${MARK}_${MODE}/Without_input_corrected"
    mkdir -p "${OUT_DIR}"

    echo "==================== Processing ${MARK} (${MODE}) ===================="

    for CELL in "${CELLTYPES[@]}"; do
        local TREAT="${MARK_DIR}/${MARK}_${CELL}.bam"

        echo "========== Running MACS2 for: ${CELL} (${MARK}) =========="

        if [[ ! -f "$TREAT" ]]; then
            echo "Missing treatment BAM for ${CELL}. Skipping..."
            continue
        fi

        # Parameters
        local QVAL=0.05
        local ARGS="-q $QVAL --keep-dup all"

        if [[ "$MODE" == "broad" ]]; then
            ARGS="$ARGS --broad --broad-cutoff $QVAL"
        fi

        macs2 callpeak \
            -t "$TREAT" \
            -f BAMPE \
            -g "$GENOME" \
            -n "${MARK}_${CELL}" \
            --outdir "$OUT_DIR" \
            $ARGS

        echo "Finished: ${MARK}_${CELL}"
        echo
    done
}

# -------------------------------------------------------------------------------------------
# Run for both broad and narrow marks
# -------------------------------------------------------------------------------------------
for MARK in "${BROAD_MARKS[@]}"; do
    run_macs2_no_ctrl "$MARK" "broad"
done

for MARK in "${NARROW_MARKS[@]}"; do
    run_macs2_no_ctrl "$MARK" "narrow"
done


#!/bin/bash
###############################################################################################################
### MACS2 Peak Calling with Control (Input BAMs)
### Multi-histone processing: Broad + Narrow marks
### MouseBrain real BAMs
###############################################################################################################

# Define histone marks
BROAD_MARKS=("H3K27ac" "H3K27me3" "H3K36me3")
NARROW_MARKS=("H3K27ac" "H3K4me3" "Olig2" "Rad21")

# Function: define cell types per histone mark
get_celltypes_for_mark() {
  local MARK="$1"
  case "$MARK" in
    H3K27ac)
      CELLTYPES=("Astrocytes" "mOL" "OEC" "OPC" "VLMC")
      ;;
    H3K27me3)
      CELLTYPES=("Astrocytes" "Microglia" "mOL" "Neurons1" "Neurons3" "OEC" "OPC" "VLMC")
      ;;
    H3K36me3)
      CELLTYPES=("Astrocytes" "mOL" "OEC" "OPC")
      ;;
    H3K4me3)
      CELLTYPES=("Astrocytes" "Microglia" "mOL" "Neurons1" "Neurons2" "Neurons3" "OEC" "OPC" "VLMC")
      ;;
    Olig2)
      CELLTYPES=("Astrocytes" "mOL" "OEC" "Unknown")
      ;;
    Rad21)
      CELLTYPES=("Astrocytes" "mOL" "OEC" "Unknown")
      ;;
    *)
      echo "Unknown mark: $MARK"
      CELLTYPES=()
      ;;
  esac
}

# Base directories
BASE_DIR="/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam"
OUT_BASE="${BASE_DIR}/MouseBrain_peakbed/MACS2_peakbed_corrected"
GENOME=1.87e9  # Mouse genome size for MACS2

mkdir -p "${OUT_BASE}"

# -------------------------------------------------------------------------------------------
# Function: Run MACS2 for all cell types under a given histone mark
# -------------------------------------------------------------------------------------------
run_macs2_all_ctrl() {
    local MARK=$1
    local MODE=$2
    local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams"
    local OUT_DIR="${OUT_BASE}/With_input_corrected"
    mkdir -p "${OUT_DIR}"

    get_celltypes_for_mark "$MARK"
    echo "==================== Processing ${MARK} (${MODE}) ===================="

    for CELL in "${CELLTYPES[@]}"; do
        local TREAT="${MARK_DIR}/${MARK}_${CELL}.bam"
        local CTRL="${MARK_DIR}/input_${CELL}.bam"

        echo "========== Running MACS2 for: ${CELL} (${MARK}) =========="

        if [[ ! -f "$TREAT" || ! -f "$CTRL" ]]; then
            echo "Missing treatment or control BAM for ${CELL}. Skipping..."
            continue
        fi

        # Set parameters
        local QVAL=0.05
        local ARGS="-q $QVAL --keep-dup all"

        if [[ "$MODE" == "broad" ]]; then
            ARGS="$ARGS --broad --broad-cutoff $QVAL"
        fi

        macs2 callpeak \
            -t "$TREAT" \
            -c "$CTRL" \
            -f BAMPE \
            -g "$GENOME" \
            -n "${MARK}_${CELL}" \
            --outdir "$OUT_DIR" \
            $ARGS

        echo "Finished: ${MARK}_${CELL}"
        echo
    done
}

# -------------------------------------------------------------------------------------------
# Run for both broad and narrow marks
# -------------------------------------------------------------------------------------------
for MARK in "${BROAD_MARKS[@]}"; do
    run_macs2_all_ctrl "$MARK" "broad"
done 

for MARK in "${NARROW_MARKS[@]}"; do
    run_macs2_all_ctrl "$MARK" "narrow"
done



#!/bin/bash
###############################################################################################################
### MACS2 Peak Calling WITHOUT Control (NO Input BAMs)
### Multi-histone processing: Broad + Narrow marks
### MouseBrain real BAMs
###############################################################################################################

# Define histone marks
BROAD_MARKS=("H3K27ac" "H3K27me3" "H3K36me3")
NARROW_MARKS=("H3K27ac" "H3K4me3" "Olig2" "Rad21")

# Function: define cell types per histone mark
get_celltypes_for_mark() {
  local MARK="$1"
  case "$MARK" in
    H3K27ac)
      CELLTYPES=("Astrocytes" "mOL" "OEC" "OPC" "VLMC")
      ;;
    H3K27me3)
      CELLTYPES=("Astrocytes" "Microglia" "mOL" "Neurons1" "Neurons3" "OEC" "OPC" "VLMC")
      ;;
    H3K36me3)
      CELLTYPES=("Astrocytes" "mOL" "OEC" "OPC")
      ;;
    H3K4me3)
      CELLTYPES=("Astrocytes" "Microglia" "mOL" "Neurons1" "Neurons2" "Neurons3" "OEC" "OPC" "VLMC")
      ;;
    Olig2)
      CELLTYPES=("Astrocytes" "mOL" "OEC" "Unknown")
      ;;
    Rad21)
      CELLTYPES=("Astrocytes" "mOL" "OEC" "Unknown")
      ;;
    *)
      echo " Unknown mark: $MARK"
      CELLTYPES=()
      ;;
  esac
}

# Base directories
BASE_DIR="/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam"
OUT_BASE="${BASE_DIR}/MouseBrain_peakbed/MACS2_peakbed_corrected"
GENOME=1.87e9  # Mouse genome size for MACS2

mkdir -p "${OUT_BASE}"

# -------------------------------------------------------------------------------------------
# Function: Run MACS2 for all cell types under a given histone mark (without control)
# -------------------------------------------------------------------------------------------
run_macs2_no_ctrl() {
    local MARK=$1
    local MODE=$2
    local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams"
    local OUT_DIR="${OUT_BASE}/Without_input_corrected"
    mkdir -p "${OUT_DIR}"

    get_celltypes_for_mark "$MARK"
    echo "==================== Processing ${MARK} (${MODE}) ===================="

    for CELL in "${CELLTYPES[@]}"; do
        local TREAT="${MARK_DIR}/${MARK}_${CELL}.bam"
        echo "========== Running MACS2 for: ${CELL} (${MARK}) =========="

        if [[ ! -f "$TREAT" ]]; then
            echo "Missing treatment BAM for ${CELL}. Skipping..."
            continue
        fi

        # Set parameters
        local QVAL=0.05
        local ARGS="-q $QVAL --keep-dup all"

        if [[ "$MODE" == "broad" ]]; then
            ARGS="$ARGS --broad --broad-cutoff $QVAL"
        fi

        macs2 callpeak \
            -t "$TREAT" \
            -f BAMPE \
            -g "$GENOME" \
            -n "${MARK}_${CELL}" \
            --outdir "$OUT_DIR" \
            $ARGS

        echo "Finished: ${MARK}_${CELL}"
        echo
    done
}

# -------------------------------------------------------------------------------------------
# Run for both broad and narrow marks
# -------------------------------------------------------------------------------------------
for MARK in "${BROAD_MARKS[@]}"; do
    run_macs2_no_ctrl "$MARK" "broad"
done

for MARK in "${NARROW_MARKS[@]}"; do
    run_macs2_no_ctrl "$MARK" "narrow"
done

```
````





**NOTE:** The following `bash` code for all HMs of hPBMC by using SICER2.

````bash
```bash

#!/bin/bash
################################################################################################
### SICER2 Peak Calling WITH Control (Input BAMs)
### Multi-histone processing: Broad + Narrow marks
### Human PBMC real BAMs
#################################################################################################

#------------------------------------------
# Histone marks
#------------------------------------------
BROAD_MARKS=("H3K27ac" "H3K27me3" "H3K9me3")
NARROW_MARKS=("H3K27ac" "H3K4me1" "H3K4me2" "H3K4me3")

#------------------------------------------
# Cell types
#------------------------------------------
CELLTYPES=("B" "CD4T" "CD8T" "DC" "Mono" "NK" "otherT" "other")

#------------------------------------------
# Base paths and parameters
#------------------------------------------
THREADS=32
BASE_DIR="/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam"
OUT_BASE="${BASE_DIR}/HumanPBMC_peakbed/SICER2_peakbed_corrected"
mkdir -p "${OUT_BASE}"

# Human genome parameters
GENOME="hg38"
FDR=0.05

#------------------------------------------
# Function: Run SICER2 with Control
#------------------------------------------
run_sicer2_all_input() {
    local MARK="$1"
    local MODE="$2"  # broad or narrow

    local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams_corrected/"
    local OUT_DIR="${OUT_BASE}/${MARK}_${MODE}/With_input_corrected"
    mkdir -p "$OUT_DIR"

    # Adjust window and gap for broad vs narrow
    if [[ "$MODE" == "broad" ]]; then
        WINDOW=100
        GAP=200
    else
        WINDOW=50
        GAP=100
    fi

    for CELL in "${CELLTYPES[@]}"; do
        echo "========== Running SICER2 for: $MARK (${CELL}) =========="
        CELL_OUT="${OUT_DIR}/${MARK}_${CELL}"
        mkdir -p "$CELL_OUT"

        # Check BAM files
        if [[ ! -f "${MARK_DIR}/${MARK}_${CELL}.bam" ]]; then
            echo "Treatment BAM missing: ${MARK}_${CELL}.bam. Skipping..."
            continue
        fi
        if [[ ! -f "${MARK_DIR}/input_${CELL}.bam" ]]; then
            echo "Control BAM missing: input_${CELL}.bam. Skipping..."
            continue
        fi

        # Run SICER2
        sicer \
          --treatment_file "${MARK_DIR}/${MARK}_${CELL}.bam" \
          --control_file "${MARK_DIR}/input_${CELL}.bam" \
          --species "$GENOME" \
          --fragment_size 150 \
          --window_size $WINDOW \
          --gap_size $GAP \
          --effective_genome_fraction 0.80 \
          --false_discovery_rate $FDR \
          --redundancy_threshold 1 \
          --output_directory "$CELL_OUT" \
          --cpu $THREADS \
          --significant_reads \
          > "${CELL_OUT}/sicer.log" 2>&1

        echo "Finished: ${MARK}_${CELL} (SICER2)"
        echo
    done

    #------------------------------------------
    # Peak summary
    #------------------------------------------
    echo "========= SICER2 Peak Summary for $MARK ($MODE) ========="
    for CELL in "${CELLTYPES[@]}"; do
        DIR="${OUT_DIR}/${MARK}_${CELL}"
        PEAK_FILE=$(find "$DIR" -type f -name "*-FDR0.05-island.bed" | head -n 1)
        if [[ -f "$PEAK_FILE" ]]; then
            COUNT=$(wc -l < "$PEAK_FILE")
            echo "${MARK}_${CELL}: $COUNT peaks"
        else
            echo "${MARK}_${CELL}: 0 peaks (file not found)"
        fi
    done
    echo
}

#------------------------------------------
# Run for Broad and Narrow Marks
#------------------------------------------
for MARK in "${BROAD_MARKS[@]}"; do
    run_sicer2_all_input "$MARK" "broad"
done

for MARK in "${NARROW_MARKS[@]}"; do
    run_sicer2_all_input "$MARK" "narrow"
done



#!/bin/bash
#################################################################################################
### SICER2 Peak Calling WITHOUT Control (No Input BAMs)
### Multi-histone processing: Broad + Narrow marks
### Human PBMC real BAMs
#################################################################################################

#------------------------------------------
# Histone marks
#------------------------------------------
BROAD_MARKS=("H3K27ac" "H3K27me3" "H3K9me3")
NARROW_MARKS=("H3K27ac" "H3K4me1" "H3K4me2" "H3K4me3")

#------------------------------------------
# Cell types
#------------------------------------------
CELLTYPES=("B" "CD4T" "CD8T" "DC" "Mono" "NK" "otherT" "other")

#------------------------------------------
# Base paths and parameters
#------------------------------------------
THREADS=32
BASE_DIR="/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/"
OUT_BASE="${BASE_DIR}/HumanPBMC_peakbed/SICER2_peakbed_corrected"
mkdir -p "${OUT_BASE}"

GENOME="hg38"
FDR=0.05

#------------------------------------------
# Function: Run SICER2 without Control
#------------------------------------------
run_sicer2_no_input() {
  local MARK="$1"
  local MODE="$2"  # broad or narrow

  local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams_corrected"
  local OUT_DIR="${OUT_BASE}/${MARK}_${MODE}/Without_input_corrected"
  mkdir -p "$OUT_DIR"

  # Adjust window and gap for broad vs narrow
  if [[ "$MODE" == "broad" ]]; then
    WINDOW=100; GAP=200
  else
    WINDOW=50; GAP=100
  fi

  for CELL in "${CELLTYPES[@]}"; do
    echo "========== Running SICER2 for: $MARK (${CELL}) =========="

    CELL_OUT="${OUT_DIR}/${MARK}_${CELL}"
    mkdir -p "$CELL_OUT"

    BAM="${MARK_DIR}/${MARK}_${CELL}.bam"
    if [[ ! -f "$BAM" ]]; then
      echo "Missing BAM file: $BAM"
      continue
    fi

    sicer \
      --treatment_file "$BAM" \
      --species "$GENOME" \
      --fragment_size 150 \
      --window_size "$WINDOW" \
      --gap_size "$GAP" \
      --effective_genome_fraction 0.80 \
      --false_discovery_rate "$FDR" \
      --redundancy_threshold 1 \
      --output_directory "$CELL_OUT" \
      --cpu "$THREADS" \
      --e_value 1000 \
      --significant_reads \
      > "${CELL_OUT}/sicer.log" 2>&1

    # Convert .scoreisland to .bed
    SCORE_FILE=$(find "$CELL_OUT" -name "*scoreisland" | head -n 1)
    if [[ -f "$SCORE_FILE" ]]; then
      BED_FILE="${SCORE_FILE%.scoreisland}-FDR${FDR}-island.bed"
      awk 'BEGIN{OFS="\t"} {print $1, $2, $3, ".", $5, "."}' "$SCORE_FILE" > "$BED_FILE"
    fi

    echo "Finished: ${MARK}_${CELL}"
    echo
  done

  echo "========= SICER2 Peak Summary for $MARK ========="
  for CELL in "${CELLTYPES[@]}"; do
    DIR="${OUT_DIR}/${MARK}_${CELL}"
    PEAK_FILE=$(find "$DIR" -type f -name "*-FDR${FDR}-island.bed" | head -n 1)
    if [[ -f "$PEAK_FILE" ]]; then
      count=$(wc -l < "$PEAK_FILE")
      echo "${MARK}_${CELL}: $count peaks"
    else
      echo "${MARK}_${CELL}: 0 peaks"
    fi
  done
  echo
}

#------------------------------------------
# Run SICER2 for Broad and Narrow Marks
#------------------------------------------
for MARK in "${BROAD_MARKS[@]}"; do
  run_sicer2_no_input "$MARK" "broad"
done

for MARK in "${NARROW_MARKS[@]}"; do
  run_sicer2_no_input "$MARK" "narrow"
done


#!/bin/bash
###############################################################################################
### SICER2 Peak Calling WITH Control (with Input BAMs)
### Multi-histone processing: Broad + Narrow marks
### MouseBrain real BAMs
################################################################################################

#------------------------------------------
# Histone marks
#------------------------------------------
BROAD_MARKS=("H3K27ac" "H3K27me3" "H3K36me3")
NARROW_MARKS=("H3K27ac" "H3K4me3" "Olig2" "Rad21")

#------------------------------------------
# Function: define cell types per histone mark
#------------------------------------------
get_celltypes_for_mark() {
  local MARK="$1"
  case "$MARK" in
    H3K27ac) CELLTYPES=("Astrocytes" "mOL" "OEC" "OPC" "VLMC") ;;
    H3K27me3) CELLTYPES=("Astrocytes" "Microglia" "mOL" "Neurons1" "Neurons3" "OEC" "OPC" "VLMC") ;;
    H3K36me3) CELLTYPES=("Astrocytes" "mOL" "OEC" "OPC") ;;
    H3K4me3) CELLTYPES=("Astrocytes" "Microglia" "mOL" "Neurons1" "Neurons2" "Neurons3" "OEC" "OPC" "VLMC") ;;
    Olig2|Rad21) CELLTYPES=("Astrocytes" "mOL" "OEC" "Unknown") ;;
    *) echo "Unknown mark: $MARK"; CELLTYPES=() ;;
  esac
}

#------------------------------------------
# Base paths and parameters
#------------------------------------------
THREADS=32
BASE_DIR="/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam"
OUT_BASE="${BASE_DIR}/MouseBrain_peakbed/SICER2_peakbed"
mkdir -p "${OUT_BASE}"

# Mouse genome
GENOME="mm10"
FDR=0.05

#------------------------------------------
# Function: Run SICER2 with Control
#------------------------------------------
run_sicer2_all_input() {
    local MARK="$1"
    local MODE="$2"  # broad or narrow

    local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams"
    local OUT_DIR="${OUT_BASE}/${MARK}_${MODE}/With_input_corrected"
    mkdir -p "$OUT_DIR"

    # Adjust window and gap for broad vs narrow
    if [[ "$MODE" == "broad" ]]; then
        WINDOW=100
        GAP=200
    else
        WINDOW=50
        GAP=100
    fi

    for CELL in "${CELLTYPES[@]}"; do
        echo "========== Running SICER2 for: $MARK (${CELL}) =========="
        CELL_OUT="${OUT_DIR}/${MARK}_${CELL}"
        mkdir -p "$CELL_OUT"

        # Check BAM files
        if [[ ! -f "${MARK_DIR}/${MARK}_${CELL}.bam" ]]; then
            echo "Treatment BAM missing: ${MARK}_${CELL}.bam. Skipping..."
            continue
        fi
        if [[ ! -f "${MARK_DIR}/input_${CELL}.bam" ]]; then
            echo "Control BAM missing: input_${CELL}.bam. Skipping..."
            continue
        fi

        # Run SICER2
        sicer \
          --treatment_file "${MARK_DIR}/${MARK}_${CELL}.bam" \
          --control_file "${MARK_DIR}/input_${CELL}.bam" \
          --species "$GENOME" \
          --fragment_size 150 \
          --window_size $WINDOW \
          --gap_size $GAP \
          --effective_genome_fraction 0.80 \
          --false_discovery_rate $FDR \
          --redundancy_threshold 1 \
          --output_directory "$CELL_OUT" \
          --cpu $THREADS \
          --significant_reads \
          > "${CELL_OUT}/sicer.log" 2>&1

        echo "Finished: ${MARK}_${CELL} (SICER2)"
        echo
    done

    #------------------------------------------
    # Peak summary
    #------------------------------------------
    echo "========= SICER2 Peak Summary for $MARK ($MODE) ========="
    for CELL in "${CELLTYPES[@]}"; do
        DIR="${OUT_DIR}/${MARK}_${CELL}"
        PEAK_FILE=$(find "$DIR" -type f -name "*-FDR0.05-island.bed" | head -n 1)
        if [[ -f "$PEAK_FILE" ]]; then
            COUNT=$(wc -l < "$PEAK_FILE")
            echo "${MARK}_${CELL}: $COUNT peaks"
        else
            echo "${MARK}_${CELL}: 0 peaks (file not found)"
        fi
    done
    echo
}

#------------------------------------------
# Run SICER2 for Broad Marks
#------------------------------------------
for MARK in "${BROAD_MARKS[@]}"; do
    get_celltypes_for_mark "$MARK"
    run_sicer2_all_input "$MARK" "broad"
done

#------------------------------------------
# Run SICER2 for Narrow Marks
#------------------------------------------
for MARK in "${NARROW_MARKS[@]}"; do
    get_celltypes_for_mark "$MARK"
    run_sicer2_all_input "$MARK" "narrow"
done


#!/bin/bash
##################################################################################################
### SICER2 Peak Calling WITHOUT Control (No Input BAMs)
### Multi-histone processing: Broad + Narrow marks
### MouseBrain real BAMs
##################################################################################################

#------------------------------------------
# Histone marks
#------------------------------------------
BROAD_MARKS=("H3K27ac" "H3K27me3" "H3K36me3")
NARROW_MARKS=("H3K27ac" "H3K4me3" "Olig2" "Rad21")

#------------------------------------------
# Function: define cell types per histone mark
#------------------------------------------
get_celltypes_for_mark() {
  local MARK="$1"
  case "$MARK" in
    H3K27ac) CELLTYPES=("Astrocytes" "mOL" "OEC" "OPC" "VLMC") ;;
    H3K27me3) CELLTYPES=("Astrocytes" "Microglia" "mOL" "Neurons1" "Neurons3" "OEC" "OPC" "VLMC") ;;
    H3K36me3) CELLTYPES=("Astrocytes" "mOL" "OEC" "OPC") ;;
    H3K4me3) CELLTYPES=("Astrocytes" "Microglia" "mOL" "Neurons1" "Neurons2" "Neurons3" "OEC" "OPC" "VLMC") ;;
    Olig2|Rad21) CELLTYPES=("Astrocytes" "mOL" "OEC" "Unknown") ;;
    *) echo "Unknown mark: $MARK"; CELLTYPES=() ;;
  esac
}

#------------------------------------------
# Base paths and parameters
#------------------------------------------
THREADS=32
BASE_DIR="/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam"
OUT_BASE="${BASE_DIR}/MouseBrain_peakbed/SICER2_peakbed"
GENOME="mm10"
FDR=0.05
mkdir -p "$OUT_BASE"

#------------------------------------------
# Function: Run SICER2 without Control
#------------------------------------------
run_sicer2_no_input() {
  local MARK="$1"
  local MODE="$2"  # broad or narrow

  local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams"
  local OUT_DIR="${OUT_BASE}/${MARK}_${MODE}/Without_input_corrected"
  mkdir -p "$OUT_DIR"

  # Adjust window and gap
  if [[ "$MODE" == "broad" ]]; then
    WINDOW=100; GAP=200
  else
    WINDOW=50; GAP=100
  fi

  for CELL in "${CELLTYPES[@]}"; do
    echo "= Running SICER2 for: $MARK (${CELL}) ="

    CELL_OUT="${OUT_DIR}/${MARK}_${CELL}"
    mkdir -p "$CELL_OUT"

    BAM="${MARK_DIR}/${MARK}_${CELL}.bam"
    if [[ ! -f "$BAM" ]]; then
      echo "Missing BAM file: $BAM"
      continue
    fi

    sicer \
      --treatment_file "$BAM" \
      --species "$GENOME" \
      --fragment_size 150 \
      --window_size "$WINDOW" \
      --gap_size "$GAP" \
      --effective_genome_fraction 0.80 \
      --false_discovery_rate "$FDR" \
      --redundancy_threshold 1 \
      --output_directory "$CELL_OUT" \
      --cpu "$THREADS" \
      --e_value 1000 \
      --significant_reads \
      > "${CELL_OUT}/sicer.log" 2>&1

    # Convert .scoreisland to .bed
    SCORE_FILE=$(find "$CELL_OUT" -name "*scoreisland" | head -n 1)
    if [[ -f "$SCORE_FILE" ]]; then
      BED_FILE="${SCORE_FILE%.scoreisland}-FDR${FDR}-island.bed"
      awk 'BEGIN{OFS="\t"} {print $1, $2, $3, ".", $5, "."}' "$SCORE_FILE" > "$BED_FILE"
    fi

    echo "Finished: ${MARK}_${CELL}"
    echo
  done

  echo "========= SICER2 Peak Summary for $MARK ========="
  for CELL in "${CELLTYPES[@]}"; do
    DIR="${OUT_DIR}/${MARK}_${CELL}"
    PEAK_FILE=$(find "$DIR" -type f -name "*-FDR${FDR}-island.bed" | head -n 1)
    if [[ -f "$PEAK_FILE" ]]; then
      count=$(wc -l < "$PEAK_FILE")
      echo "${MARK}_${CELL}: $count peaks"
    else
      echo "${MARK}_${CELL}: 0 peaks"
    fi
  done
  echo
}

#------------------------------------------
# Run SICER2 for Broad and Narrow Marks
#------------------------------------------
for MARK in "${BROAD_MARKS[@]}"; do
  get_celltypes_for_mark "$MARK"
  run_sicer2_no_input "$MARK" "broad"
done

for MARK in "${NARROW_MARKS[@]}"; do
  get_celltypes_for_mark "$MARK"
  run_sicer2_no_input "$MARK" "narrow"
done
```
````



### 3.3 Generate pseudo-input BAM files

&#x20;

To accurately call peaks using scCUT&Tag peak finders such as DROMPA+, MACS2, GoPeaks, Genrich, HOMER, SEACR, and SICER2 a control (input) BAM file is generally required for each sample to model the background read distribution and reduce false positives. However, in datasets derived from sorted or single‑cell populations, a true genomic input library is often unavailable. To circumvent this limitation, we generated pseudo‑input control BAM files for each cell type by pooling all remaining cell types within the same histone modification (HM) or transcription factor (TF) dataset. This strategy provides a representative background sample that accounts for sequencing biases across the experiment.



For each cell type, the procedure follows three main steps: (1) exclude the target cell‑type BAM, (2) merge all other valid cell‑type BAMs, and (3) sort and index the resulting merged file. The entire process is automated using the Bash script provided below, which demonstrates the workflow for the H3K4me3 mark in the hPBMC dataset. The script includes robust file‑validation checks (“samtools”) to automatically skip any corrupted BAM files, ensuring that the final merged outputs are of high quality.





````bash
```bash
###=====================================================
###    Making Input bam
###=====================================================
cd ~/project_scHMTF/GSE195725_processed_data/splitbam_realbam/H3K4me3/split_celltype_bams_corrected

#!/bin/bash

# Function to create input for a specific cell type
create_input() {
    local exclude_cell="$1"
    local input_name="input_${exclude_cell}.bam"
    local temp_merged="input_merged_${exclude_cell}.bam"
    
    echo "Creating input for $exclude_cell (excluding H3K4me3_${exclude_cell}.bam)"
    
    # Get all valid BAM files except the one to exclude
    input_files=()
    for bam in H3K4me3_*.bam; do
        # Skip the excluded cell type and any corrupted files
        if [[ "$bam" == "H3K4me3_${exclude_cell}.bam" ]] || [[ "$bam" == "H3K4me3_temp_bam.bam" ]]; then
            continue
        fi
        
        # Check if BAM file is valid
        if samtools quickcheck "$bam" 2>/dev/null; then
            input_files+=("$bam")
            echo "  Including: $bam"
        else
            echo "  Skipping corrupted file: $bam"
        fi
    done
    
    # Check if we have any valid files to merge
    if [ ${#input_files[@]} -eq 0 ]; then
        echo "  No valid BAM files found for input creation"
        return 1
    fi
    
    echo "  Merging ${#input_files[@]} files..."
    
    # Merge all other cell types
    samtools merge "$temp_merged" "${input_files[@]}"
    samtools sort -o "$input_name" "$temp_merged"
    samtools index "$input_name"
    rm -f "$temp_merged"
    
    read_count=$(samtools view -c "$input_name" 2>/dev/null || echo "0")
    echo "Created: $input_name ($read_count reads)"
}

# First, check and remove the problematic temp_bam file
if [[ -f "H3K4me3_temp_bam.bam" ]]; then
    echo "Removing problematic file: H3K4me3_temp_bam.bam"
    rm -f H3K4me3_temp_bam.bam
fi

# List available BAM files
echo "Available BAM files:"
for bam in H3K4me3_*.bam; do
    if samtools quickcheck "$bam" 2>/dev/null; then
        reads=$(samtools view -c "$bam" 2>/dev/null || echo "corrupted")
        echo "$bam ($reads reads)"
    else
        echo "$bam (corrupted)"
    fi
done

echo ""

# Create input for each cell type
create_input "B"
create_input "CD4T" 
create_input "CD8T"
create_input "DC"
create_input "Mono"
create_input "NK"
create_input "other"
create_input "otherT"

echo "=== All input files created ==="
ls -lh input_*.bam 2>/dev/null || echo "No input files were created"
```
````



### 3.4 Call peaks using the pseudo-input BAM files

&#x20;

To identify enriched genomic regions, we performed peak calling on the cell‑type‑specific BAM files for both the hPBMC and mBrain datasets using pseudo input control BAM files. We employed seven distinct peak callers for the hPBMC dataset, we explored H3K27me3 and H3K9me3 as broad marks; H3K4me1, H3K4me2, and H3K4me3 as narrow (sharp) marks; and H3K27ac as both broad and narrow. For the mBrain dataset, we extracted H3K27me3 and H3K36me3 as broad, and H3K4me3, Olig2, and Rad21 as narrow, with H3K27ac again run in both modes. To ensure comparability across methods, we standardized the key input parameters for each caller according to the peak type, as detailed below:



1. For narrow (sharp) peak calling: DROMPA+ (--pthre_internal 5, --pthre_enrich 4); Genrich (-a 200, -l 100, -g 100, -p 0.01); GoPeaks (default parameters); HOMER (factor mode); MACS2 (-q 0.05, --keep-dup all); SEACR (relaxed, nrom); and SICER2 (--window_size 50, --gap 100, FDR 0.05).

2. For broad peak calling: DROMPA+ (--pthre_internal 4, --pthre_enrich 3); Genrich (-a 100, -l 500, -g 1000, -p 0.05); GoPeaks (--broad); HOMER (histone mode); MACS2 (--broad, -q 0.05, --keep-dup all); SEACR (stringent, nrom); and SICER2 (--window_size 100, --gap 200).

 

Following peak calling with each tool, we generated and subsequently explored the resulting peak BED files for downstream analysis.

As an illustration, the following “bash” workflows demonstrate peak calling using pseudo input BAM files for the H3K4me3 hPBMC dataset using MACS2 and SICER2, respectively. To execute these scripts, a dedicated Conda environment containing all required dependencies must first be created and activated.





**NOTE:** The following `bash` code for H3K4me3 of hPBMC by using MACS2.

````bash
```bash
#!/bin/bash
##################################################
### MACS2 Peak Calling with pseudo input control 
### Multi-histone processing: Broad + Narrow marks
### Human PBMC real BAMs
##################################################

# Define histone marks
BROAD_MARKS=("")
NARROW_MARKS=("H3K4me3")

# Define cell types
CELLTYPES=("B" "CD4T" "CD8T" "DC" "Mono" "NK" "otherT" "other")

# Base directories
THREADS=32
BASE_DIR="/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam"
OUT_BASE="${BASE_DIR}/HumanPBMC_peakbed/MACS2_peakbed"
GENOME=2.7e9  # Human genome size for MACS2

mkdir -p "${OUT_BASE}"

# -----------------------------------------------------
# Function: Run MACS2 for all cell types
#  under a given histone mark
# -----------------------------------------------------
run_macs2_all_ctrl() {
    local MARK=$1
    local MODE=$2
    local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams_corrected"
    local OUT_DIR="${OUT_BASE}/${MARK}_${MODE}/With_input_corrected"
    mkdir -p "${OUT_DIR}"

    echo "== Processing ${MARK} (${MODE}) =="

    for CELL in "${CELLTYPES[@]}"; do
        local TREAT="${MARK_DIR}/${MARK}_${CELL}.bam"
        local CTRL="${MARK_DIR}/input_${CELL}.bam"

        echo "==== Running MACS2 for: ${CELL} (${MARK}) ====="

        if [[ ! -f "$TREAT" || ! -f "$CTRL" ]]; then
            echo "Missing treatment or control BAM for ${CELL}. Skipping..."
            continue
        fi

        # Set parameters
        local QVAL=0.05
        local ARGS=" -q $QVAL --keep-dup all"

        if [[ "$MODE" == "broad" ]]; then
            ARGS="$ARGS --broad --broad-cutoff $QVAL"
        fi

        macs2 callpeak \
            -t "$TREAT" \
            -c "$CTRL" \
            -f BAMPE \
            -g "$GENOME" \
            -n "${MARK}_${CELL}" \
            --outdir "$OUT_DIR" \
            $ARGS

        echo "Finished: ${MARK}_${CELL}"
        echo
    done
}

#-----------------------------------------------------
# Run for both broad and narrow marks
#-----------------------------------------------------
for MARK in "${BROAD_MARKS[@]}"; do
    run_macs2_all_ctrl "$MARK" "broad"
done

for MARK in "${NARROW_MARKS[@]}"; do
    run_macs2_all_ctrl "$MARK" "narrow"
done
```
````



**NOTE:** The following `bash` code for H3K4me3 of hPBMC by using SICER2.



````bash
```bash 

#!/bin/bash
########################################################
### SICER2 Peak Calling using pseudo input control 
### Multi-histone processing: Broad + Narrow marks
### Human PBMC real BAMs
######################################################### 
#------------------------------------------
# Histone marks
#------------------------------------------
BROAD_MARKS=("")
NARROW_MARKS=("H3K4me3")

#------------------------------------------
# Cell types
#------------------------------------------
CELLTYPES=("B" "CD4T" "CD8T" "DC" "Mono" "NK" "otherT" "other")

#------------------------------------------
# Base paths and parameters
#------------------------------------------
THREADS=32
BASE_DIR="/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam"
OUT_BASE="${BASE_DIR}/HumanPBMC_peakbed/SICER2_peakbed_corrected"
mkdir -p "${OUT_BASE}"

# Human genome parameters
GENOME="hg38"
FDR=0.05

#------------------------------------------
# Function: Run SICER2 with Control
#------------------------------------------
run_sicer2_all_input() {
    local MARK="$1"
    local MODE="$2"  # broad or narrow

    local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams_corrected/"
    local OUT_DIR="${OUT_BASE}/${MARK}_${MODE}/With_input_corrected"
    mkdir -p "$OUT_DIR"

    # Adjust window and gap for broad vs narrow
    if [[ "$MODE" == "broad" ]]; then
        WINDOW=100
        GAP=200
    else
        WINDOW=50
        GAP=100
    fi

    for CELL in "${CELLTYPES[@]}"; do
        echo "=== Running SICER2 for: $MARK (${CELL}) ==="
        CELL_OUT="${OUT_DIR}/${MARK}_${CELL}"
        mkdir -p "$CELL_OUT"

        # Check BAM files
        if [[ ! -f "${MARK_DIR}/${MARK}_${CELL}.bam" ]]; then
            echo "Treatment BAM missing: ${MARK}_${CELL}.bam. Skipping..."
            continue
        fi
        if [[ ! -f "${MARK_DIR}/input_${CELL}.bam" ]]; then
            echo "Control BAM missing: input_${CELL}.bam. Skipping..."
            continue
        fi

        # Run SICER2
        sicer \
          --treatment_file "${MARK_DIR}/${MARK}_${CELL}.bam" \
          --control_file "${MARK_DIR}/input_${CELL}.bam" \
          --species "$GENOME" \
          --fragment_size 150 \
          --window_size $WINDOW \
          --gap_size $GAP \
          --effective_genome_fraction 0.80 \
          --false_discovery_rate $FDR \
          --redundancy_threshold 1 \
          --output_directory "$CELL_OUT" \
          --cpu $THREADS \
          --significant_reads \
          > "${CELL_OUT}/sicer.log" 2>&1

        echo "Finished: ${MARK}_${CELL} (SICER2)"
        echo
    done

    #------------------------------------------
    # Peak summary
    #------------------------------------------
    echo "========= SICER2 Peak Summary for $MARK ($MODE) ========="
    for CELL in "${CELLTYPES[@]}"; do
        DIR="${OUT_DIR}/${MARK}_${CELL}"
        PEAK_FILE=$(find "$DIR" -type f -name "*-FDR0.05-island.bed" | head -n 1)
        if [[ -f "$PEAK_FILE" ]]; then
            COUNT=$(wc -l < "$PEAK_FILE")
            echo "${MARK}_${CELL}: $COUNT peaks"
        else
            echo "${MARK}_${CELL}: 0 peaks (file not found)"
        fi
    done
    echo
}

#------------------------------------------
# Run for Broad and Narrow Marks
#------------------------------------------
for MARK in "${BROAD_MARKS[@]}"; do
    run_sicer2_all_input "$MARK" "broad"
done

for MARK in "${NARROW_MARKS[@]}"; do
    run_sicer2_all_input "$MARK" "narrow"
done
```
````



### 3.5 Generate simulated BAM files

&#x20;

To generate cell‑type‑specific simulated treatment and input control samples for the hPBMC and mBrain datasets, we implemented a custom computational pipeline using R (version 4.3.3) with the Signac, Seurat, GenomeInfoDb, and data.table packages, alongside Samtools (version 1.13) for binary conversion (see **Note 4**  for full dependency details). The procedure is described below for the hPBMC dataset as a representative example; the same workflow was applied to the mBrain data.



First, for each of the six HMs under investigation (H3K27ac, H3K27me3, H3K4me1, H3K4me2, H3K4me3, and H3K9me3) for hPBMC and (H3K27ac, H3K27me3, H3K36me3, H3K4me3, Olig2, and Rad21) for mBrain, the corresponding chromatin assay object was loaded into the R environment. A fragment object was instantiated from the cell‑specific fragment file (TSV.GZ format) using the “CreateFragmentObject” function, retaining only the cells present in the assay. A new “ChromatinAssay” was then constructed by combining the tile‑level count matrix with the fragment object, using a manually curated “Seqinfo” object derived from the “hg38” and “mm10” genome reference (chromosome sizes obtained from the UCSC table). Cell identities were assigned according to the “predicted.celltype.l” for hPBMC and “celltype” for mBrain annotation metadata field (**Note 5** ).

Within each HMs mark, the pipeline iterated over all distinct cell types. For every cell type, two paired‑end SAM files were generated:

(i) a treatment file containing fragments from cells belonging to that type, and

(ii) a pseudo input control file containing fragments from all remaining (off‑target) cells of the same HM/TF.

This separation was achieved by subsetting the fragment data frame according to cell barcodes. Each subset was converted to SAM format using a custom function that writes paired‑end reads (read length 50 bp, template length derived from fragment coordinates), with the standard SAM header prepended. The iteration over cell types was parallelised across multiple cores to reduce processing time (the number of cores is configurable based on available hardware).



````R
```r
########################################################
#### Making BAM file for all HMs  for Human PBMC dataset
#### CellTypes wise
########################################################
###-----------------------------------------------------
## Required packages
###-----------------------------------------------------
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(data.table)
###-----------------------------------------------------
## Make Sure the histones name
###-----------------------------------------------------
histones<-c("H3K27ac", "H3K27me3", "H3K4me1", "H3K4me2", "H3K4me3", "H3K9me3")
###-----------------------------------------------------
## Setting the "for loop" 
###-----------------------------------------------------
for(hist in histones){

histone <- readRDS(paste0("/home/wahid/project_scHMTF/GSE195725_processed_data/", hist, ".rds"))
# --- 2. Read and create Fragment object ---
fragments <- CreateFragmentObject(
  path = paste0("/home/wahid/project_scHMTF/GSE195725_processed_data/", hist, "_fragments.tsv.gz"),
  cells = colnames(histone)
)
# --- 3. Create new ChromatinAssay and reassign it to object ---
library(GenomeInfoDb)
# Read your genome file (chromosome sizes)
genome_file <- "/home/wahid/project_scHMTF/GSE195725_processed_data/ref/hg38.genome"
chrom_sizes <- read.table(genome_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(chrom_sizes) <- c("seqnames", "seqlengths")
# Create header
header <- c("@HD\tVN:1.6\tSO:unsorted")
header <- c(header, paste0("@SQ\tSN:", chrom_sizes$seqnames, "\tLN:", chrom_sizes$seqlengths))
# Create Seqinfo manually
hg38_seqinfo <- Seqinfo(seqnames = chrom_sizes$seqnames,
                        seqlengths = chrom_sizes$seqlengths,
                        genome = "hg38")

                        # Now pass it to CreateChromatinAssay
chrom_assay <- CreateChromatinAssay(
  counts = GetAssayData(histone[["tiles"]], slot = "counts"),
  fragments = fragments,
  sep = c("-", "-"),
  genome = hg38_seqinfo  # Custom genome
)
histone[["tiles"]] <- chrom_assay

# Path to original fragment file
frag_path <- Fragments(histone)[[1]]@path

# Read compressed fragment file
frags <- fread(
  cmd = paste("zcat", frag_path),
  header = FALSE,
  colClasses = c("character", "integer", "integer", "character", "integer")
)

colnames(frags) <- c("chr", "start", "end", "cell", "score")
# === 2. Output directory ===
outdir <- paste0("/home/wahid/project_scHMTF/GSE195725_processed_data/BAM/celltypeswisebam_l1/",hist)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
# === 4. Paired-end SAM writer ===
write_pe_sam <- function(df, outfile) {
  sam_lines <- lapply(1:nrow(df), function(i) {
    row <- df[i]
    read_id <- paste0(row$cell, "_", i)
    chr <- row$chr
    start1 <- format(as.integer(row$start) + 1, scientific = FALSE)
    start2 <- format(as.integer(row$end) - 50, scientific = FALSE)
    if (as.integer(start2) < 1) start2 <- 1
    flag1 <- 99
    flag2 <- 147
    r1 <- paste(read_id, flag1, chr, start1, 255, "50M", "=", start2, as.integer(start2) - as.integer(start1), "*", "*", sep = "\t")
    r2 <- paste(read_id, flag2, chr, start2, 255, "50M", "=", start1, as.integer(start1) - as.integer(start2), "*", "*", sep = "\t")
    c(r1, r2)
  })
  writeLines(c(header, unlist(sam_lines)), con = outfile)
}
#Setting Idents instead of celltype.l1 table(H3K27me3$predicted.celltype.l1)
Idents(histone) <- histone$predicted.celltype.l1
# === 5. Iterate over all clusters ===
all_clusters <- unique(Idents(histone))
library(parallel)
# Delete old corrupted SAMs if needed
unlink(list.files(outdir, pattern = "\\.sam$", full.names = TRUE))
# Re-run parallel SAM generation
library(parallel)
mclapply(all_clusters, function(cluster) {
  message("Processing cluster: ", cluster)
  treatment_cells <- colnames(histone)[Idents(histone) == cluster]
  input_control_cells   <- colnames(histone)[Idents(histone) != cluster]
  frag_treat <- frags[cell %in% treatment_cells]
  frag_ctrl  <- frags[cell %in% input_control_cells]
  sam_treat <- file.path(outdir, paste0("treatment_", cluster, ".sam"))
  sam_ctrl  <- file.path(outdir, paste0("input_control_", cluster, ".sam"))
  write_pe_sam(frag_treat, sam_treat)
  write_pe_sam(frag_ctrl, sam_ctrl)
}, mc.cores = 64)
}
```
````

Following SAM generation, all .sam files were converted to the compressed binary BAM format using the samtools view command with the -S -b flags on `bash`.&#x20;

````bash
```bash
###################################################
### Sam to bam in bash
###################################################
for file in *.sam; do
  # Remove extension and handle spaces using quotes
  base=$(basename "$file" .sam)
  # Convert SAM to BAM
  samtools view -S -b "$file" > "${base}.bam"
done
```
````

After that, pseudo input control BAM files were generated by merging the remaining BAM files from all cell types of a specific cell type. This was achieved by merging, for a given target cell type, all the BAM files corresponding to the input control samples of all other cell types from the same HM.



### 3.6 Call peaks from the simulated BAM files (without an input control)

&#x20;

To identify enriched genomic regions, we performed peak calling on the cell‑type‑specific simulated BAM files for both the hPBMC and mBrain datasets under without input control. We employed seven distinct peak callers for the hPBMC dataset, we processed H3K27me3 and H3K9me3 as broad marks; H3K4me1, H3K4me2, and H3K4me3 as narrow (sharp) marks; and H3K27ac as both broad and narrow. For the mBrain dataset, we processed H3K27me3 and H3K36me3 as broad, and H3K4me3, Olig2, and Rad21 as narrow, with H3K27ac again run in both modes. To ensure comparability across methods, we standardized the key input parameters for each caller according to the peak type, as detailed below:



1. For narrow (sharp) peak calling: DROMPA+ (--pthre_internal 5); Genrich (-a 200, -l 100, -g 100, -p 0.01); GoPeaks (default parameters); HOMER (factor mode); MACS2 (-q 0.05, --keep-dup all); SEACR (relaxed, 0.01, non); and SICER2 (--window_size 50, --gap 100, -p 0.01).

2. For broad peak calling: DROMPA+ (--pthre_internal 4); Genrich (-a 100, -l 500, -g 1000, -p 0.05); GoPeaks (--broad); HOMER (histone mode); MACS2 (--broad, -q 0.05, --keep-dup all); SEACR (stringent, 0.01, non); and SICER2 (--window_size 100, --gap 200, -p 0.05).

 

Following peak calling with each tool, we generated and subsequently explored the resulting peak BED files for downstream analysis.

As an illustration, the following “bash” workflows demonstrate peak calling for the H3K4me3 hPBMC dataset using MACS2 and SICER2, respectively. To execute these scripts, a dedicated Conda environment containing all required dependencies must first be created and activated.





**NOTE:** The following `bash` code for H3K4me3 of hPBMC simulated BAM by using MACS2.

````bash
```bash
#!/bin/bash
########################################################
### Multi-histone processing: Broad + Narrow marks
### Human PBMC real BAMs
########################################################

# Define histone marks
BROAD_MARKS=("")
NARROW_MARKS=("H3K4me3")

# Define cell types
CELLTYPES=("B" "CD4T" "CD8T" "DC" "Mono" "NK" "otherT" "other")

# Base directories
BASE_DIR="simulated BAM directories"
OUT_BASE="${BASE_DIR}/output_peakbed of simulated BAM/MACS2_peakbed"
GENOME=2.7e9  # Human genome size for MACS2

mkdir -p "${OUT_BASE}"

# ----------------------------------------------------
# Function: Run MACS2 for all cell types 
#under a given histone mark (no control)
# ----------------------------------------------------
run_macs2_no_ctrl() {
    local MARK=$1
    local MODE=$2
    local MARK_DIR="${BASE_DIR}/${MARK}/simulated_BAM"
    local OUT_DIR="${OUT_BASE}/${MARK}_${MODE}/Without_input_corrected"
    mkdir -p "${OUT_DIR}"

    echo "== Processing ${MARK} (${MODE}) ======"

    for CELL in "${CELLTYPES[@]}"; do
        local TREAT="${MARK_DIR}/${MARK}_${CELL}.bam"

        echo "=== Running MACS2 for: ${CELL} (${MARK}) =="

        if [[ ! -f "$TREAT" ]]; then
            echo "Missing treatment BAM for ${CELL}. Skipping..."
            continue
        fi

        # Parameters
        local QVAL=0.05
        local ARGS="-q $QVAL --keep-dup all"

        if [[ "$MODE" == "broad" ]]; then
            ARGS="$ARGS --broad --broad-cutoff $QVAL"
        fi

        macs2 callpeak \
            -t "$TREAT" \
            -f BAMPE \
            -g "$GENOME" \
            -n "${MARK}_${CELL}" \
            --outdir "$OUT_DIR" \
            $ARGS

        echo "Finished: ${MARK}_${CELL}"
        echo
    done
}

# ----------------------------------------------------
# Run for both broad and narrow marks
# ----------------------------------------------------
for MARK in "${BROAD_MARKS[@]}"; do
    run_macs2_no_ctrl "$MARK" "broad"
done

for MARK in "${NARROW_MARKS[@]}"; do
    run_macs2_no_ctrl "$MARK" "narrow"
done
```
````

**NOTE:** The following `bash` code for H3K4me3 of hPBMC simulated BAM by using SICER2.

````bash
```bash

#!/bin/bash
#########################################################
### SICER2 Peak Calling WITHOUT Control (No Input BAMs)
### Multi-histone processing: Broad + Narrow marks
### Human PBMC real BAMs
#########################################################

#------------------------------------------
# Histone marks
#------------------------------------------
BROAD_MARKS=("")
NARROW_MARKS=("H3K4me3")

#------------------------------------------
# Cell types
#------------------------------------------
CELLTYPES=("B" "CD4T" "CD8T" "DC" "Mono" "NK" "otherT" "other")

#------------------------------------------
# Base paths and parameters
#------------------------------------------
THREADS=32
BASE_DIR="simulated BAM directories"
OUT_BASE="${BASE_DIR}/output_peakbed of simulated BAM/SICER2_peakbed"
mkdir -p "${OUT_BASE}"

GENOME="hg38"
FDR=0.05

#------------------------------------------
# Function: Run SICER2 without Control
#------------------------------------------
run_sicer2_no_input() {
  local MARK="$1"
  local MODE="$2"  # broad or narrow

  local MARK_DIR="${BASE_DIR}/${MARK}/simulated_BAM"
  local OUT_DIR="${OUT_BASE}/${MARK}_${MODE}/Without_input_corrected"
  mkdir -p "$OUT_DIR"

  # Adjust window and gap for broad vs narrow
  if [[ "$MODE" == "broad" ]]; then
    WINDOW=100; GAP=200
  else
    WINDOW=50; GAP=100
  fi

  for CELL in "${CELLTYPES[@]}"; do
    echo "= Running SICER2 for: $MARK (${CELL}) ="

    CELL_OUT="${OUT_DIR}/${MARK}_${CELL}"
    mkdir -p "$CELL_OUT"

    BAM="${MARK_DIR}/${MARK}_${CELL}.bam"
    if [[ ! -f "$BAM" ]]; then
      echo "Missing BAM file: $BAM"
      continue
    fi

    sicer \
      --treatment_file "$BAM" \
      --species "$GENOME" \
      --fragment_size 150 \
      --window_size "$WINDOW" \
      --gap_size "$GAP" \
      --effective_genome_fraction 0.80 \
      --false_discovery_rate "$FDR" \
      --redundancy_threshold 1 \
      --output_directory "$CELL_OUT" \
      --cpu "$THREADS" \
      --e_value 1000 \
      --significant_reads \
      > "${CELL_OUT}/sicer.log" 2>&1

    # Convert .scoreisland to .bed
    SCORE_FILE=$(find "$CELL_OUT" -name "*scoreisland" | head -n 1)
    if [[ -f "$SCORE_FILE" ]]; then
      BED_FILE="${SCORE_FILE%.scoreisland}-FDR${FDR}-island.bed"
      awk 'BEGIN{OFS="\t"} {print $1, $2, $3, ".", $5, "."}' "$SCORE_FILE" > "$BED_FILE"
    fi

    echo "Finished: ${MARK}_${CELL}"
    echo
  done

  echo "========= SICER2 Peak Summary for $MARK ========="
  for CELL in "${CELLTYPES[@]}"; do
    DIR="${OUT_DIR}/${MARK}_${CELL}"
    PEAK_FILE=$(find "$DIR" -type f -name "*-FDR${FDR}-island.bed" | head -n 1)
    if [[ -f "$PEAK_FILE" ]]; then
      count=$(wc -l < "$PEAK_FILE")
      echo "${MARK}_${CELL}: $count peaks"
    else
      echo "${MARK}_${CELL}: 0 peaks"
    fi
  done
  echo
}

#------------------------------------------
# Run SICER2 for Broad and Narrow Marks
#------------------------------------------
for MARK in "${BROAD_MARKS[@]}"; do
  run_sicer2_no_input "$MARK" "broad"
done

for MARK in "${NARROW_MARKS[@]}"; do
  run_sicer2_no_input "$MARK" "narrow"
done
```
````

### 3.7 Summarize the peak numbers and peak widths

&#x20;

In this section, we describe the summarized of cell‑type‑specific peak count matrices for human PBMC and mouse brain datasets, under both with‑input and without‑input control strategies, following peak calling with the seven tools. Below is an example “bash” script that processes the human PBMC dataset with input control. For each combination of histone mark (or transcription factor), peak calling method, and cell type, the script calculates: (1) Number of peaks (total peak count per sample) and (2) Mean peak width (average length of peaks in base pairs).



````bash
```bash
#!/bin/bash
##############################################################
# Calculation of Number of Peaks and Peak Widths for HumanPBMC 
##############################################################
Histones=("H3K27ac-b" "H3K27ac-s" "H3K27me3" "H3K4me1" "H3K4me2" "H3K4me3" "H3K9me3")
Methods=("DROMPAplus" "Genrich" "GoPeaks" "HOMER" "MACS2" "SEACR" "SICER2")
# Input directory containing BED files
PEAK_DIR="peak bed directory of hPBMC with input control"
# Output CSV file
output_file="Output summary of hPBMC with input control as csv"
# Write header
echo "Histone/TF,Method,Sample,Number_of_Peaks,Mean_Peak_Width" > "$output_file"
# Loop through each Histone/TF
for histone in "${Histones[@]}"; do
    # Loop through each Method
    for method in "${Methods[@]}"; do
        # Find all BED files matching the method and histone
        for bed in "$PEAK_DIR"/${method}_${histone}*.bed; do
            # Skip if file does not exist
            [[ -f "$bed" ]] || continue
            # Extract sample name
            sample=$(basename "$bed" .bed)          
            # Calculate number of peaks
            num_peaks=$(wc -l < "$bed")
            # Calculate mean peak width
            mean_width=$(awk '{sum += $3 - $2} END {if (NR>0) print int(sum / NR); else print 0}' "$bed")
            # Append to output CSV
            echo "${histone},${method},${sample},${num_peaks},${mean_width}" >> "$output_file"
        done
    done
done
echo "Peak summary calculation completed. Output saved to $output_file"
```

````

We compared the performance of seven peak calling tools including DROMPA+, Genrich, GoPeaks, HOMER, MACS2, SEACR, and SICER2 using 10 HMs and 2 HTs scCUT&Tag datasets derived from hPBMC and mbrain. For each dataset, we evaluated both real and simulated data under two experimental conditions: with pseudo input control and without pseudo input control. As for example H3K4me3 scCUT&Tag datasets derived from human PBMC CD8T cells and mouse brain mOL.

&#x20;

Performance was assessed based on the total number of peaks called, which are presented as log10‑transformed values in **Figure 2**. In human PBMC, for both real and simulated datasets, DROMPA+ and Genrich yielded the highest peak counts when an input control was available, whereas HOMER and MACS2 produced the lowest counts. In the absence of an input control, SEACR, DROMPA+, and Genrich demonstrated the most consistent performance, exhibiting only a moderate reduction in peak numbers. For mouse brain real data, GoPeaks and HOMER generated the highest peak counts in the presence of an input control, while SEACR recorded the lowest; however, without input control, SEACR produced the highest peak count, followed closely by SICER2, which remained among the top performers. For mbrain simulated data, DROMPA+ and Genrich again achieved the highest peak numbers when input was available, with HOMER showing the lowest counts; in the absence of input control, GoPeaks resulted in the lowest peak count among all tools. Collectively, these findings indicate that the choice of an optimal peak caller depends on multiple factors, including species, data type (real versus simulated), and the availability of an input control. 



For genomic track visualization, we used the online version of the Integrative Genomics Viewer (IGV) [16]. After obtaining peak BED files from both real and simulated datasets of human PBMC and mouse brain under conditions with and without input controls, we selected the H3K4me3 peak BED files and corresponding bigWig files for CD8T cells (human PBMC) and mOL cells (mouse brain) (see **Note 6**). All files were uploaded to the IGV web portal to visually compare the performance of seven peak callers across real and simulated data, as well as with and without input controls (**Figure 3**). The following genomic regions were examined: human PBMC, Chr1:52,570,837–53,103,680; mouse brain, Chr1:52,434,885–52,967,728.



IGV inspection revealed marked differences in peak calling performance among the seven tools, driven by data type (real vs. simulated) and input control inclusion. Without input controls, several callers, especially MACS2, HOMER, and SICER2 tended to produce dense, low‑confidence interval calls. Including matched input controls suppressed this background noise, yielding sharper, localized peaks that aligned with robust signals in the master track.

Real biological datasets displayed sharp, well‑resolved enrichment peaks. In contrast, simulated tracks produced broad, uniform fragment calls across nearly all tools. While DROMPA+, GoPeaks, and Gennich maintained stable, stringent peak calls under both real and simulated conditions, other tools showed heightened sensitivity to background fluctuations, particularly in the absence of control tracks.



### 4. Notes

&#x20;

1. Software prerequisites can be installed using On Linux, the recommended approach is to install R from the official Comprehensive R Archive Network (CRAN) repository to obtain the most recent stable version. For Ubuntu/Debian systems, this involves adding the CRAN repository and its GPG key, then installing the r-base and r-base-dev packages via bioconda.

2. For the “cellranger-atac” count pipeline, reference genomes must be downloaded from 10x Genomics. The human (hg38/GRCh38) and mouse (mm10) reference genomes used in this study were refdata-cellranger-arc-GRCh38-2020-A-2.0.0 and refdata-cellranger-arc-mm10-2020-A-2.0.0, respectively. These can be obtained from the official 10x Genomics support website or via direct download links provided by the vendor. The files are large (approximately 14 GB for the human genome) and should be downloaded over a stable internet connection. After download, the compressed archives must be extracted using tar -xzvf, generating a directory containing the reference genome structure. When running cellranger-atac count or cellranger-arc count, the path to the extracted reference directory is specified using the --reference parameter.

3. For each histone modification (HM) and transcription factor (TF) dataset, the number of cells per cell type was determined. Using the cell-type-specific percentages, the pooled BAM files were split into cell-type-specific BAM files. Read counts were then extracted from these split BAM files and log10-transformed to improve data distribution and visualization. The resulting log10-transformed read counts are shown in **Figure 1**.

4. Conda Environment and Dependency Specifications To ensure full reproducibility of the simulated BAM file generation pipeline, all software dependencies were managed using a dedicated Conda environment.

   name: simulate-bam

   channels:

   - conda-forge

   - bioconda

   - defaults

   - bcftools=1.16

   - bedtools=2.31.0

5. Cell type identities were assigned using two distinct annotation strategies based on available metadata: for the hPBMC dataset, labels were derived from the predicted.celltype.l field generated by Signac label transfer, which distinguishes major immune lineages (B, CD4T, CD8T, DC, Mono, NK, otherT and other) to ensure sufficient cell counts per type across all histone modifications; for the mBrain dataset, the celltype metadata provided with the original processed data was used directly, containing expert-curated annotations (Astrocytes, OEC, OPC, Microglia, mOL, Neurons1-3, and VLMC). Prior to fragment subsetting, a validation step excluded any cells lacking valid annotations or marked as “Unknown” to avoid contaminating the pseudo-input controls, and cell-type counts were tabulated per histone mark—with all major types meeting the minimum threshold of 50 cells demonstrating the pipeline’s adaptability to diverse single-cell datasets through configurable metadata field specification.

6. To validate and visually inspect the peak calls generated by the various peak callers, we employed the Integrative Genomics Viewer (IGV), a high-performance, interactive tool for the visual exploration of genomic datasets. The IGV web application was accessed via the official website at https://igv.org/app/ [16]. Upon launching the application, the appropriate reference genome was selected for each dataset: the human reference genome GRCh38 (hg38) was chosen for the human peripheral blood mononuclear cell (hPBMC) samples, while the mouse reference genome mm10 was selected for the mouse brain (mBrain) samples. Subsequently, the peak BED files—representing the genomic intervals identified as enriched regions—were loaded into IGV for all conditions examined. This included peak files derived from both the real and simulated datasets, generated across all peak caller tools evaluated in this study. For each of these conditions, visualizations were performed both with and without the inclusion of a pseudo-input control, enabling a comparative assessment of peak calling performance and the impact of background signal correction.



### References

&#x20;

1. Bartosovic M, Kabbe M, Castelo-Branco G (2021) Single-cell CUT&Tag profiles histone modifications and transcription factors in complex tissues. Nature biotechnology 39 (7):825-835. doi:10.1038/s41587-021-00869-9

2. Grosselin K, Durand A, Marsolier J, Poitou A, Marangoni E, Nemati F, Dahmani A, Lameiras S, Reyal F, Frenoy O, Pousse Y, Reichen M, Woolfe A, Brenan C, Griffiths AD, Vallot C, Gérard A (2019) High-throughput single-cell ChIP-seq identifies heterogeneity of chromatin states in breast cancer. Nature genetics 51 (6):1060-1066. doi:10.1038/s41588-019-0424-9

3. Wu J, Wahiduzzaman M, Yin P, Sun P, Chen H, Ding Y, Wang J (2026) Advances in scCUT&Tag and computational analysis for single-cell gene regulatory element mapping. Briefings in bioinformatics 27 (1). doi:10.1093/bib/bbag015

4. Nakato R, Sakata T (2021) Methods for ChIP-seq analysis: A practical workflow and advanced applications. Methods 187:44-53. doi:10.1016/j.ymeth.2020.03.005

5. Smith JP, Sheffield NC (2020) Analytical Approaches for ATAC-seq Data Analysis. Curr Protoc Hum Genet 106 (1):e101. doi:10.1002/cphg.101

6. Yashar WM, Kong G, VanCampen J, Curtiss BM, Coleman DJ, Carbone L, Yardimci GG, Maxson JE, Braun TP (2022) GoPeaks: histone modification peak calling for CUT&Tag. Genome Biol 23 (1):144. doi:10.1186/s13059-022-02707-w

7. Heinz S, Benner C, Spann N, Bertolino E, Lin YC, Laslo P, Cheng JX, Murre C, Singh H, Glass CK (2010) Simple combinations of lineage-determining transcription factors prime cis-regulatory elements required for macrophage and B cell identities. Mol Cell 38 (4):576-589. doi:10.1016/j.molcel.2010.05.004

8. Zhang Y, Liu T, Meyer CA, Eeckhoute J, Johnson DS, Bernstein BE, Nusbaum C, Myers RM, Brown M, Li W, Liu XS (2008) Model-based analysis of ChIP-Seq (MACS). Genome Biol 9 (9):R137. doi:10.1186/gb-2008-9-9-r137

9. Meers MP, Tenenbaum D, Henikoff S (2019) Peak calling by Sparse Enrichment Analysis for CUT&RUN chromatin profiling. Epigenetics Chromatin 12 (1):42. doi:10.1186/s13072-019-0287-4

10. Eder T, Grebien F (2022) Comprehensive assessment of differential ChIP-seq tools guides optimal algorithm selection. Genome Biol 23 (1):119. doi:10.1186/s13059-022-02686-y

11. Frankish A, Carbonell-Sala S, Diekhans M, Jungreis I, Loveland JE, Mudge JM, Sisu C, Wright JC, Arnan C, Barnes I, Banerjee A, Bennett R, Berry A, Bignell A, Boix C, Calvet F, Cerdán-Vélez D, Cunningham F, Davidson C, Donaldson S, Dursun C, Fatima R, Giorgetti S, Giron CG, Gonzalez JM, Hardy M, Harrison PW, Hourlier T, Hollis Z, Hunt T, James B, Jiang Y, Johnson R, Kay M, Lagarde J, Martin FJ, Gómez LM, Nair S, Ni P, Pozo F, Ramalingam V, Ruffier M, Schmitt BM, Schreiber JM, Steed E, Suner MM, Sumathipala D, Sycheva I, Uszczynska-Ratajczak B, Wass E, Yang YT, Yates A, Zafrulla Z, Choudhary JS, Gerstein M, Guigo R, Hubbard TJP, Kellis M, Kundaje A, Paten B, Tress ML, Flicek P (2023) GENCODE: reference annotation for the human and mouse genomes in 2023. Nucleic acids research 51 (D1):D942-d949. doi:10.1093/nar/gkac1071

12. Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, Whitwham A, Keane T, McCarthy SA, Davies RM, Li H (2021) Twelve years of SAMtools and BCFtools. GigaScience 10 (2). doi:10.1093/gigascience/giab008

13. Hao Y, Stuart T, Kowalski MH, Choudhary S, Hoffman P, Hartman A, Srivastava A, Molla G, Madad S, Fernandez-Granda C, Satija R (2024) Dictionary learning for integrative, multimodal and scalable single-cell analysis. Nature biotechnology 42 (2):293-304. doi:10.1038/s41587-023-01767-y

14. Stuart T, Srivastava A, Madad S, Lareau CA, Satija R (2021) Single-cell chromatin state analysis with Signac. Nat Methods 18 (11):1333-1341. doi:10.1038/s41592-021-01282-5

15. Zhang B, Srivastava A, Mimitou E, Stuart T, Raimondi I, Hao Y, Smibert P, Satija R (2022) Characterizing cellular heterogeneity in chromatin state with scCUT&Tag-pro. Nature biotechnology 40 (8):1220-1230. doi:10.1038/s41587-022-01250-0

16. Robinson JT, Thorvaldsdottir H, Turner D, Mesirov JP (2023) igv.js: an embeddable JavaScript implementation of the Integrative Genomics Viewer (IGV). Bioinformatics (Oxford, England) 39 (1). doi:10.1093/bioinformatics/btac830







