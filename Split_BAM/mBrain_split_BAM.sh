
```bash
########################################################
##### Splitting bam for H3K4me3 of mBrain
########################################################
cd ~/project_scHMTF/GSE157637_processed_data/splitbam_realbam/H3K4me3  
#!/bin/bash

INPUT_BAM="H3K4me3.bam"
OUTPUT_DIR="split_celltype_bams"

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
                cb_value = arr[3]
                if(cb_value in barcodes) {
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
