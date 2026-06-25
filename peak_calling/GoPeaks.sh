#!/bin/bash
###############################################################################################################
### GoPeaks Peak Calling with Control (Input BAMs)
### Multi-histone processing: Broad + Narrow marks
### Human PBMC real BAMs
###############################################################################################################

BROAD_MARKS=("H3K27ac" "H3K27me3" "H3K9me3")
NARROW_MARKS=("H3K27ac" "H3K4me1" "H3K4me2" "H3K4me3")
CELLTYPES=("B" "CD4T" "CD8T" "DC" "Mono" "NK" "otherT" "other")

THREADS=32
BASE_DIR="/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam"
OUT_BASE="${BASE_DIR}/HumanPBMC_peakbed/GoPeaks_peakbed"
CHROMSIZE="/home/wahid/project_scHMTF/GSE195725_processed_data/ref/hg38.chrom.sizes"

mkdir -p "${OUT_BASE}"

###--------------------------------------------------
### Function: Run GoPeaks with Control (Input BAM)
###--------------------------------------------------
run_gopeaks_with_ctrl() {
    local MARK=$1
    local MODE=$2
    local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams"
    local OUT_DIR="${OUT_BASE}/${MARK}_${MODE}/With_input"

    mkdir -p "${OUT_DIR}"

    echo "#######--------------------------------------------------######"
    echo "### Running GoPeaks for ${MARK} (${MODE} peaks)"
    echo "#######--------------------------------------------------######"

    for CELL in "${CELLTYPES[@]}"; do
        echo "========== Processing: ${MARK}_${CELL} =========="

        local BAM="${MARK_DIR}/${MARK}_${CELL}.bam"
        local CTRL="${MARK_DIR}/input_${CELL}.bam"
        local PREFIX="${OUT_DIR}/${MARK}_${CELL}_${MODE}"

        if [[ -f "$BAM" && -f "$CTRL" ]]; then
            if [[ "$MODE" == "broad" ]]; then
                gopeaks \
                    -b "$BAM" \
                    -c "$CTRL" \
                    -s "$CHROMSIZE" \
                    -o "$PREFIX" \
                    --broad \
                    && echo "✅ Done: ${MARK}_${CELL}" \
                    || echo "❌ Failed: ${MARK}_${CELL}"
            else
                gopeaks \
                    -b "$BAM" \
                    -c "$CTRL" \
                    -s "$CHROMSIZE" \
                    -o "$PREFIX" \
                    && echo "✅ Done: ${MARK}_${CELL}" \
                    || echo "❌ Failed: ${MARK}_${CELL}"
            fi
        else
            echo "⚠️ Missing BAM or Control file for ${MARK}_${CELL}"
        fi
    done
}

###--------------------------------------------------
### Run for broad histone marks
###--------------------------------------------------
for MARK in "${BROAD_MARKS[@]}"; do
    run_gopeaks_with_ctrl "$MARK" "broad"
done

###--------------------------------------------------
### Run for narrow histone marks
###--------------------------------------------------
for MARK in "${NARROW_MARKS[@]}"; do
    run_gopeaks_with_ctrl "$MARK" "narrow"
done


#!/bin/bash
###========================================================================================================####
### GoPeaks Peak Calling without Control (No Input BAMs)
### Multi-histone processing: Broad + Narrow marks
### Human PBMC real BAMs
###========================================================================================================####

BROAD_MARKS=("H3K27ac" "H3K27me3" "H3K9me3")
NARROW_MARKS=("H3K27ac" "H3K4me1" "H3K4me2" "H3K4me3")
CELLTYPES=("B" "CD4T" "CD8T" "DC" "Mono" "NK" "otherT" "other")

BASE_DIR="/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam"
OUT_BASE="${BASE_DIR}/HumanPBMC_peakbed/GoPeaks_peakbed"
CHROMSIZE="/home/wahid/project_scHMTF/GSE195725_processed_data/ref/hg38.chrom.sizes"

mkdir -p "${OUT_BASE}"

###--------------------------------------------------
### Function: Run GoPeaks without Control
###--------------------------------------------------
run_gopeaks_no_ctrl() {
    local MARK=$1
    local MODE=$2
    local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams"
    local OUT_DIR="${OUT_BASE}/${MARK}_${MODE}/Without_input"

    mkdir -p "${OUT_DIR}"

    echo "#########--------------------------------------------------#######"
    echo "### Running GoPeaks for ${MARK} (${MODE} peaks; no control)"
    echo "#########--------------------------------------------------########"

    for CELL in "${CELLTYPES[@]}"; do
        echo "========== Processing: ${MARK}_${CELL} =========="

        local BAM="${MARK_DIR}/${MARK}_${CELL}.bam"
        local PREFIX="${OUT_DIR}/${MARK}_${CELL}_${MODE}"

        if [[ -f "$BAM" ]]; then
            if [[ "$MODE" == "broad" ]]; then
                gopeaks \
                    -b "$BAM" \
                    -s "$CHROMSIZE" \
                    -o "$PREFIX" \
                    --broad \
                    && echo "✅ Done: ${MARK}_${CELL}" \
                    || echo "❌ Failed: ${MARK}_${CELL}"
            else
                gopeaks \
                    -b "$BAM" \
                    -s "$CHROMSIZE" \
                    -o "$PREFIX" \
                    && echo "✅ Done: ${MARK}_${CELL}" \
                    || echo "❌ Failed: ${MARK}_${CELL}"
            fi
        else
            echo "⚠️ Missing BAM file for ${MARK}_${CELL}"
        fi
    done
}

###--------------------------------------------------
### Run for broad histone marks
###--------------------------------------------------
for MARK in "${BROAD_MARKS[@]}"; do
    run_gopeaks_no_ctrl "$MARK" "broad"
done

###--------------------------------------------------
### Run for narrow histone marks
###--------------------------------------------------
for MARK in "${NARROW_MARKS[@]}"; do
    run_gopeaks_no_ctrl "$MARK" "narrow"
done



#!/bin/bash
###############################################################################################################
### GoPeaks Peak Calling with Control (Input BAMs)
### Multi-histone processing: Broad + Narrow marks
### Mouse Brain Real BAMs
###############################################################################################################

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
      echo "⚠️ Unknown mark: $MARK"
      CELLTYPES=()
      ;;
  esac
}

THREADS=32
BASE_DIR="/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam"
OUT_BASE="${BASE_DIR}/MouseBrain_peakbed/GoPeaks_peakbed"
CHROMSIZE="/home/wahid/project_scHMTF/GSE157637_processed_data/ref/mm10.chrom.sizes"

mkdir -p "${OUT_BASE}"

###--------------------------------------------------
### Function: Run GoPeaks with Control (Input BAM)
###--------------------------------------------------
run_gopeaks_with_ctrl() {
    local MARK=$1
    local MODE=$2
    local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams"
    local OUT_DIR="${OUT_BASE}/${MARK}_${MODE}/With_input2"

    mkdir -p "${OUT_DIR}"

    echo "#######--------------------------------------------------#####"
    echo "### Running GoPeaks for ${MARK} (${MODE} peaks)"
    echo "#######--------------------------------------------------#####"

    for CELL in "${CELLTYPES[@]}"; do
        echo "========== Processing: ${MARK}_${CELL} =========="

        local BAM="${MARK_DIR}/${MARK}_${CELL}.bam"
        local CTRL="${MARK_DIR}/input_${CELL}.bam"
        local PREFIX="${OUT_DIR}/${MARK}_${CELL}_${MODE}"

        if [[ -f "$BAM" && -f "$CTRL" ]]; then
            if [[ "$MODE" == "broad" ]]; then
                gopeaks \
                    -b "$BAM" \
                    -c "$CTRL" \
                    -s "$CHROMSIZE" \
                    -o "$PREFIX" \
                    --broad \
                    && echo "✅ Done: ${MARK}_${CELL}" \
                    || echo "❌ Failed: ${MARK}_${CELL}"
            else
                gopeaks \
                    -b "$BAM" \
                    -c "$CTRL" \
                    -s "$CHROMSIZE" \
                    -o "$PREFIX" \
                    && echo "✅ Done: ${MARK}_${CELL}" \
                    || echo "❌ Failed: ${MARK}_${CELL}"
            fi
        else
            echo "⚠️ Missing BAM or Control file for ${MARK}_${CELL}"
        fi
    done
}

###--------------------------------------------------
### Run for broad histone marks
###--------------------------------------------------
for MARK in "${BROAD_MARKS[@]}"; do
    get_celltypes_for_mark "$MARK"
    run_gopeaks_with_ctrl "$MARK" "broad"
done

###--------------------------------------------------
### Run for narrow histone marks
###--------------------------------------------------
for MARK in "${NARROW_MARKS[@]}"; do
    get_celltypes_for_mark "$MARK"
    run_gopeaks_with_ctrl "$MARK" "narrow"
done


#!/bin/bash
###############################################################################################################
### GoPeaks Peak Calling WITHOUT Control (Input BAMs)
### Multi-histone processing: Broad + Narrow marks
### Mouse Brain Real BAMs
###############################################################################################################

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
      echo "⚠️ Unknown mark: $MARK"
      CELLTYPES=()
      ;;
  esac
}

BASE_DIR="/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam"
OUT_BASE="${BASE_DIR}/MouseBrain_peakbed/GoPeaks_peakbed"
CHROMSIZE="/home/wahid/project_scHMTF/GSE157637_processed_data/ref/mm10.chrom.sizes"

mkdir -p "${OUT_BASE}"

###--------------------------------------------------
### Function: Run GoPeaks WITHOUT Control
###--------------------------------------------------
run_gopeaks_no_ctrl() {
    local MARK=$1
    local MODE=$2
    local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams"
    local OUT_DIR="${OUT_BASE}/${MARK}_${MODE}/Without_input"

    mkdir -p "${OUT_DIR}"

    echo "#######--------------------------------------------------####"
    echo "### Running GoPeaks for ${MARK} (${MODE} peaks) WITHOUT control"
    echo "#######--------------------------------------------------####"

    for CELL in "${CELLTYPES[@]}"; do
        echo "========== Processing: ${MARK}_${CELL} =========="

        local BAM="${MARK_DIR}/${MARK}_${CELL}.bam"
        local PREFIX="${OUT_DIR}/${MARK}_${CELL}_${MODE}"

        if [[ -f "$BAM" ]]; then
            if [[ "$MODE" == "broad" ]]; then
                gopeaks \
                    -b "$BAM" \
                    -s "$CHROMSIZE" \
                    -o "$PREFIX" \
                    --broad \
                    && echo "✅ Done: ${MARK}_${CELL}" \
                    || echo "❌ Failed: ${MARK}_${CELL}"
            else
                gopeaks \
                    -b "$BAM" \
                    -s "$CHROMSIZE" \
                    -o "$PREFIX" \
                    && echo "✅ Done: ${MARK}_${CELL}" \
                    || echo "❌ Failed: ${MARK}_${CELL}"
            fi
        else
            echo "⚠️ Missing BAM file for ${MARK}_${CELL}"
        fi
    done
}

###--------------------------------------------------
### Run for broad histone marks
###--------------------------------------------------
for MARK in "${BROAD_MARKS[@]}"; do
    get_celltypes_for_mark "$MARK"
    run_gopeaks_no_ctrl "$MARK" "broad"
done

###--------------------------------------------------
### Run for narrow histone marks
###--------------------------------------------------
for MARK in "${NARROW_MARKS[@]}"; do
    get_celltypes_for_mark "$MARK"
    run_gopeaks_no_ctrl "$MARK" "narrow"
done
