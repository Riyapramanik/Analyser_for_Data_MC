#!/bin/bash
# create_file_lists.sh - Script to create file lists for each directory

# Define your base directories
DIRECTORIES=(
    "/eos/user/p/pravesh/run3/2023/MC/WW_TuneCP5_13p6TeV_pythia8/WW/250821_151512/0000/"
)

# Create output directory for file lists
mkdir -p filelists_WW

# Function to create file list for each directory
create_file_list() {
    local dir=$1
    local index=$2
    
    echo "Creating file list for directory $index: $dir"
    
    # Extract a meaningful name for the file list
    local dirname=$(basename $(dirname $(dirname $dir)))
    local subdir=$(basename $dir)
    local filename="filelists_WW/filelist_${dirname}_${subdir}_job${index}.txt"
    
    # Find all ROOT files and write to file list
    find $dir -name "*.root" -type f > $filename
    
    # Count files
    local file_count=$(wc -l < $filename)
    echo "Found $file_count ROOT files in $dir"
    echo "File list saved as: $filename"
    
    echo $filename
}

# Create file lists for each directory
echo "Creating file lists for Condor jobs..."
job_index=0

for dir in "${DIRECTORIES[@]}"; do
    if [ -d "$dir" ]; then
        create_file_list "$dir" $job_index
        ((job_index++))
    else
        echo "Warning: Directory does not exist: $dir"
    fi
done

