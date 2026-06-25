#!/bin/bash
###############################################################################################################
### HOMER Peak Calling WITH Control (Input BAMs)
### Multi-histone processing: Broad + Narrow marks
### Human PBMC real BAMs
###############################################################################################################

BROAD_MARKS=("H3K27ac" "H3K27me3" "H3K9me3")
NARROW_MARKS=("H3K27ac" "H3K4me1" "H3K4me2" "H3K4me3")
CELLTYPES=("B" "CD4T" "CD8T" "DC" "Mono" "NK" "otherT" "other")

BASE_DIR="/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam"
OUT_BASE="${BASE_DIR}/HumanPBMC_peakbed/HOMER_peakbed"
TAG_BASE="${BASE_DIR}/tag_directories"

mkdir -p "${OUT_BASE}" "${TAG_BASE}"

run_homer_all_ctrl() {
    local MARK=$1
    local MODE=$2
    local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams"
    local OUT_DIR="${OUT_BASE}/${MARK}_${MODE}/With_input"

    mkdir -p "${OUT_DIR}"

    for CELL in "${CELLTYPES[@]}"; do
        echo "========== Running HOMER for: $CELL (${MARK}) =========="

        TREATMENT_BAM="${MARK_DIR}/${MARK}_${CELL}.bam"
        CONTROL_BAM="${MARK_DIR}/input_${CELL}.bam"

        if [[ ! -f "$TREATMENT_BAM" || ! -f "$CONTROL_BAM" ]]; then
            echo "❌ Missing treatment or control BAM for $CELL. Skipping..."
            continue
        fi

        TREATMENT_TAG="${TAG_BASE}/${MARK}_${CELL}_treatment_tagdir"
        CONTROL_TAG="${TAG_BASE}/${MARK}_${CELL}_control_tagdir"

        echo "🔹 Creating tag directories..."
        makeTagDirectory "$TREATMENT_TAG" "$TREATMENT_BAM"
        makeTagDirectory "$CONTROL_TAG" "$CONTROL_BAM"

        echo "🔹 Running HOMER findPeaks..."
        PEAK_TXT="${OUT_DIR}/${CELL}_${MARK}_peaks.txt"

        # Use broad style for broad marks, factor for narrow
        if [[ " ${BROAD_MARKS[*]} " =~ " ${MARK} " ]]; then
            STYLE="histone"
        else
            STYLE="factor"
        fi

        findPeaks "$TREATMENT_TAG" -style "$STYLE" -i "$CONTROL_TAG" -o "$PEAK_TXT"

        if [[ -s "$PEAK_TXT" ]]; then
            echo "🔹 Converting peaks.txt to standard BED6..."
            PEAK_BED="${OUT_DIR}/${CELL}_${MARK}_peaks.bed"
            awk 'BEGIN{OFS="\t"} NR>1 {print $2, $3, $4, $1, $8, "."}' "$PEAK_TXT" > "$PEAK_BED"
            echo "✅ BED saved: $PEAK_BED"
        else
            echo "⚠️ No peaks found (empty $PEAK_TXT)"
        fi

        echo
    done

    echo "========= HOMER Peak Summary ========="
    for PEAK_BED in "${OUT_DIR}"/*_peaks.bed; do
        [[ -f "$PEAK_BED" ]] || continue
        COUNT=$(wc -l < "$PEAK_BED")
        echo "$(basename "$PEAK_BED"): $COUNT peaks"
    done
}

# Run HOMER for all marks (broad + narrow)
for MARK in "${BROAD_MARKS[@]}"; do
    run_homer_all_ctrl "$MARK" "broad"
done

for MARK in "${NARROW_MARKS[@]}"; do
    run_homer_all_ctrl "$MARK" "narrow"
done


#!/bin/bash
###############################################################################################################
### HOMER Peak Calling WITHOUT Control (Input BAMs)
### Multi-histone processing: Broad + Narrow marks
### Human PBMC real BAMs
###############################################################################################################

BROAD_MARKS=("H3K27ac" "H3K27me3" "H3K9me3")
NARROW_MARKS=("H3K27ac" "H3K4me1" "H3K4me2" "H3K4me3")
CELLTYPES=("B" "CD4T" "CD8T" "DC" "Mono" "NK" "otherT" "other")

BASE_DIR="/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam"
OUT_BASE="${BASE_DIR}/HumanPBMC_peakbed/HOMER_peakbed"
TAG_BASE="${BASE_DIR}/tag_directories2"

mkdir -p "${OUT_BASE}" "${TAG_BASE}"

run_homer_all_noctrl() {
    local MARK=$1
    local MODE=$2
    local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams"
    local OUT_DIR="${OUT_BASE}/${MARK}_${MODE}/Without_input"

    mkdir -p "${OUT_DIR}"

    for CELL in "${CELLTYPES[@]}"; do
        echo "========== Running HOMER for: $CELL (${MARK}) =========="

        TREATMENT_BAM="${MARK_DIR}/${MARK}_${CELL}.bam"

        if [[ ! -f "$TREATMENT_BAM" ]]; then
            echo "❌ Missing treatment BAM for $CELL. Skipping..."
            continue
        fi

        TREATMENT_TAG="${TAG_BASE}/${MARK}_${CELL}_treatment_tagdir"

        echo "🔹 Creating tag directory..."
        makeTagDirectory "$TREATMENT_TAG" "$TREATMENT_BAM"

        echo "🔹 Running HOMER findPeaks..."
        PEAK_TXT="${OUT_DIR}/${CELL}_${MARK}_peaks.txt"

        # Use broad style for broad marks, factor for narrow
        if [[ " ${BROAD_MARKS[*]} " =~ " ${MARK} " ]]; then
            STYLE="histone"
        else
            STYLE="factor"
        fi

        findPeaks "$TREATMENT_TAG" -style "$STYLE" -o "$PEAK_TXT"

        if [[ -s "$PEAK_TXT" ]]; then
            echo "🔹 Converting peaks.txt to standard BED6..."
            PEAK_BED="${OUT_DIR}/${CELL}_${MARK}_peaks.bed"
            awk 'BEGIN{OFS="\t"} NR>1 {print $2, $3, $4, $1, $8, "."}' "$PEAK_TXT" > "$PEAK_BED"
            echo "✅ BED saved: $PEAK_BED"
        else
            echo "⚠️ No peaks found (empty $PEAK_TXT)"
        fi

        echo
    done

    echo "========= HOMER Peak Summary ========="
    for PEAK_BED in "${OUT_DIR}"/*_peaks.bed; do
        [[ -f "$PEAK_BED" ]] || continue
        COUNT=$(wc -l < "$PEAK_BED")
        echo "$(basename "$PEAK_BED"): $COUNT peaks"
    done
}

# Run HOMER for all marks (broad + narrow)
for MARK in "${BROAD_MARKS[@]}"; do
    run_homer_all_noctrl "$MARK" "broad"
done

for MARK in "${NARROW_MARKS[@]}"; do
    run_homer_all_noctrl "$MARK" "narrow"
done





#!/bin/bash
###############################################################################################################
### HOMER Peak Calling WITH Control (Input BAMs)
### Multi-histone processing: Broad + Narrow marks
### MouseBrain real BAMs
###############################################################################################################

# Define broad and narrow marks
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

# Base directories
BASE_DIR="/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam"
OUT_BASE="${BASE_DIR}/MouseBrain_peakbed/HOMER_peakbed"
TAG_BASE="${BASE_DIR}/tag_directories"

mkdir -p "${OUT_BASE}" "${TAG_BASE}"

# Function: run HOMER for a given mark
run_homer_all_ctrl() {
    local MARK="$1"
    local MODE="$2"
    local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams"
    local OUT_DIR="${OUT_BASE}/${MARK}_${MODE}/With_input"

    mkdir -p "${OUT_DIR}"

    # Assign cell types dynamically
    get_celltypes_for_mark "$MARK"

    for CELL in "${CELLTYPES[@]}"; do
        echo "========== Running HOMER for: $CELL (${MARK}) =========="

        TREATMENT_BAM="${MARK_DIR}/${MARK}_${CELL}.bam"
        CONTROL_BAM="${MARK_DIR}/input_${CELL}.bam"

        if [[ ! -f "$TREATMENT_BAM" || ! -f "$CONTROL_BAM" ]]; then
            echo "❌ Missing treatment or control BAM for $CELL. Skipping..."
            continue
        fi

        TREATMENT_TAG="${TAG_BASE}/${MARK}_${CELL}_treatment_tagdir"
        CONTROL_TAG="${TAG_BASE}/${MARK}_${CELL}_control_tagdir"

        echo "🔹 Creating tag directories..."
        makeTagDirectory "$TREATMENT_TAG" "$TREATMENT_BAM"
        makeTagDirectory "$CONTROL_TAG" "$CONTROL_BAM"

        echo "🔹 Running HOMER findPeaks..."
        PEAK_TXT="${OUT_DIR}/${CELL}_${MARK}_peaks.txt"

        # Use broad style for broad marks, factor for narrow
        if [[ " ${BROAD_MARKS[*]} " =~ " ${MARK} " ]]; then
            STYLE="histone"
        else
            STYLE="factor"
        fi

        findPeaks "$TREATMENT_TAG" -style "$STYLE" -i "$CONTROL_TAG" -o "$PEAK_TXT"

        if [[ -s "$PEAK_TXT" ]]; then
            echo "🔹 Converting peaks.txt to standard BED6..."
            PEAK_BED="${OUT_DIR}/${CELL}_${MARK}_peaks.bed"
            awk 'BEGIN{OFS="\t"} NR>1 {print $2, $3, $4, $1, $8, "."}' "$PEAK_TXT" > "$PEAK_BED"
            echo "✅ BED saved: $PEAK_BED"
        else
            echo "⚠️ No peaks found (empty $PEAK_TXT)"
        fi

        echo
    done

    echo "========= HOMER Peak Summary ========="
    for PEAK_BED in "${OUT_DIR}"/*_peaks.bed; do
        [[ -f "$PEAK_BED" ]] || continue
        COUNT=$(wc -l < "$PEAK_BED")
        echo "$(basename "$PEAK_BED"): $COUNT peaks"
    done
}

# Run HOMER for all broad marks
for MARK in "${BROAD_MARKS[@]}"; do
    run_homer_all_ctrl "$MARK" "broad"
done

# Run HOMER for all narrow marks
for MARK in "${NARROW_MARKS[@]}"; do
    run_homer_all_ctrl "$MARK" "narrow"
done



#!/bin/bash
###############################################################################################################
### HOMER Peak Calling WITHOUT Control
### Multi-histone processing: Broad + Narrow marks
### MouseBrain real BAMs
###############################################################################################################

# Define broad and narrow marks
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

# Base directories
BASE_DIR="/home/wahid/project_scHMTF/GSE157637_processed_data/splitbam_realbam"
OUT_BASE="${BASE_DIR}/MouseBrain_peakbed/HOMER_peakbed"
TAG_BASE="${BASE_DIR}/tag_directories2"

mkdir -p "${OUT_BASE}" "${TAG_BASE}"

# Function: run HOMER for a given mark (without control)
run_homer_no_ctrl() {
    local MARK="$1"
    local MODE="$2"
    local MARK_DIR="${BASE_DIR}/${MARK}/split_celltype_bams"
    local OUT_DIR="${OUT_BASE}/${MARK}_${MODE}/Without_input"

    mkdir -p "${OUT_DIR}"

    # Assign cell types dynamically
    get_celltypes_for_mark "$MARK"

    for CELL in "${CELLTYPES[@]}"; do
        echo "========== Running HOMER for: $CELL (${MARK}) =========="

        TREATMENT_BAM="${MARK_DIR}/${MARK}_${CELL}.bam"

        if [[ ! -f "$TREATMENT_BAM" ]]; then
            echo "❌ Missing treatment BAM for $CELL. Skipping..."
            continue
        fi

        TREATMENT_TAG="${TAG_BASE}/${MARK}_${CELL}_treatment_tagdir"

        echo "🔹 Creating tag directory..."
        makeTagDirectory "$TREATMENT_TAG" "$TREATMENT_BAM"

        echo "🔹 Running HOMER findPeaks..."
        PEAK_TXT="${OUT_DIR}/${CELL}_${MARK}_peaks.txt"

        # Use broad style for broad marks, factor for narrow
        if [[ " ${BROAD_MARKS[*]} " =~ " ${MARK} " ]]; then
            STYLE="histone"
        else
            STYLE="factor"
        fi

        findPeaks "$TREATMENT_TAG" -style "$STYLE" -o "$PEAK_TXT"

        if [[ -s "$PEAK_TXT" ]]; then
            echo "🔹 Converting peaks.txt to standard BED6..."
            PEAK_BED="${OUT_DIR}/${CELL}_${MARK}_peaks.bed"
            awk 'BEGIN{OFS="\t"} NR>1 {print $2, $3, $4, $1, $8, "."}' "$PEAK_TXT" > "$PEAK_BED"
            echo "✅ BED saved: $PEAK_BED"
        else
            echo "⚠️ No peaks found (empty $PEAK_TXT)"
        fi

        echo
    done

    echo "========= HOMER Peak Summary ========="
    for PEAK_BED in "${OUT_DIR}"/*_peaks.bed; do
        [[ -f "$PEAK_BED" ]] || continue
        COUNT=$(wc -l < "$PEAK_BED")
        echo "$(basename "$PEAK_BED"): $COUNT peaks"
    done
}

# Run HOMER for all broad marks
for MARK in "${BROAD_MARKS[@]}"; do
    run_homer_no_ctrl "$MARK" "broad"
done

# Run HOMER for all narrow marks
for MARK in "${NARROW_MARKS[@]}"; do
    run_homer_no_ctrl "$MARK" "narrow"
done


