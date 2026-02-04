#!/bin/bash
# Test data generation script for pysec2pri
# Generates all export formats for all datasources using --format all

source .venv/bin/activate

OUTPUT_DIR="test_data_generation"
DATABASES="chebi hgnc hmdb ncbi uniprot"

# Archived versions for applicable datasources
declare -A ARCHIVED_VERSIONS
ARCHIVED_VERSIONS[chebi]="245"
ARCHIVED_VERSIONS[hgnc]="2023-07-01"
ARCHIVED_VERSIONS[uniprot]="2024_01"

# Create output directory
mkdir -p "$OUTPUT_DIR"

echo "=== Testing pysec2pri data generation ==="
echo "Output directory: $OUTPUT_DIR"
echo ""

# Function to run pysec2pri for a datasource with --format all
run_datasource() {
    local db=$1
    local version=$2
    local output_dir=$3
    local db_dir="${output_dir}/${db}"
    
    if [ -n "$version" ]; then
        db_dir="${db_dir}_${version}"
        echo "Processing $db version $version..."
        pysec2pri "$db" --format all --version "$version" -o "$db_dir/${db}_all.tsv"
    else
        echo "Processing $db (latest)..."
        pysec2pri "$db" --format all -o "$db_dir/${db}_all.tsv"
    fi
    
    echo "  Done: $db_dir"
}

# Export function for xargs
export -f run_datasource
export OUTPUT_DIR

echo "=== Processing latest versions (parallel) ==="
# Run all databases with --format all using xargs
echo "$DATABASES" | tr ' ' '\n' | xargs -P 4 -I {} bash -c 'run_datasource "$1" "" "$2"' _ {} "$OUTPUT_DIR"

echo ""
echo "=== Testing archived versions ==="
# Run archived versions for applicable datasources
for db in chebi hgnc uniprot; do
    version="${ARCHIVED_VERSIONS[$db]}"
    if [ -n "$version" ]; then
        run_datasource "$db" "$version" "$OUTPUT_DIR"
    fi
done

echo ""
echo "=== All data generation complete ==="
echo "Results in: $OUTPUT_DIR"
ls -la "$OUTPUT_DIR"
