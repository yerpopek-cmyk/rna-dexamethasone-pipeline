#!/usr/bin/env bash
set -euo pipefail

# Navigate to project root directory
cd "$(dirname "$0")/.."

mkdir -p data

SAMPLES=(
    "SRR1039508" "SRR1039509" "SRR1039512" "SRR1039513"
    "SRR1039516" "SRR1039517" "SRR1039520" "SRR1039521"
)

# Helper function to get EBI HTTP URL for a given accession and read number
get_ebi_url() {
    local srr=$1
    local read_num=$2
    local first6="${srr:0:6}"
    local last_char="${srr: -1}"
    local suffix="00${last_char}"
    echo "https://ftp.sra.ebi.ac.uk/vol1/fastq/${first6}/${suffix}/${srr}/${srr}_${read_num}.fastq.gz"
}

MAX_JOBS=4
pids=()

for SRR in "${SAMPLES[@]}"; do
    for READ in 1 2; do
        FILE="data/${SRR}_${READ}.fastq.gz"
        if [[ -f "$FILE" ]]; then
            echo "$FILE already exists. Skipping."
            continue
        fi
        
        URL=$(get_ebi_url "$SRR" "$READ")
        
        # Check active pids and wait if we reached max jobs
        while [ ${#pids[@]} -ge $MAX_JOBS ]; do
            # Filter out finished pids
            temp_pids=()
            for pid in "${pids[@]}"; do
                if kill -0 "$pid" 2>/dev/null; then
                    temp_pids+=("$pid")
                fi
            done
            pids=("${temp_pids[@]}")
            if [ ${#pids[@]} -ge $MAX_JOBS ]; then
                sleep 1
            fi
        done
        
        echo "Downloading $FILE in parallel..."
        wget -q -c "$URL" -O "$FILE" &
        pids+=($!)
    done
done

# Wait for remaining background processes
if [ ${#pids[@]} -gt 0 ]; then
    echo "Waiting for remaining downloads to finish..."
    wait "${pids[@]}"
fi

echo "All downloads completed successfully!"

