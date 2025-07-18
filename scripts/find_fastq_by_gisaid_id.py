#!/usr/bin/env python3
import argparse
import configparser
import sys
import shutil
from pathlib import Path
 
import pandas as pd
import re
 
"""
================================================================================
Find and Copy BAM Files by GISAID ID
================================================================================
 
Purpose:
    Given a list of GISAID EPI_ISL accession numbers, this script traces them
    back to their corresponding BAM mapping files and copies them to a
    specified output directory for further analysis.
 
Workflow:
    1. Reads a list of target GISAID IDs from an input file.
    2. Reads the sample sheet that links GISAID IDs to internal sample aliases.
    3. Reads the master metadata Excel file to get sequencing run info.
    4. Merges this information to find the SequenceRun and Barcode for each target ID.
    5. Reconstructs and verifies the raw BAM file path.
    6. Copies all found BAM files to the user-specified output directory.
 
Example Usage:
    python scripts/find_and_copy_bam_by_gisaid_id.py \
        --gisaid-ids /path/to/your/gisaid_ids.csv \
        --sample-sheet input_tsv/sam_with_ids.tsv \
        --master-excel /path/to/your/master_metadata.xls \
        --config config/sars.ini \
        --output-dir bams_for_review
"""
 
def standardize_barcode(barcode: str) -> str:
    """Standardizes a barcode string to the format 'BC' followed by two digits."""
    if pd.isna(barcode):
        return ""
    match = re.search(r'\d+', str(barcode))
    if not match:
        return ""
    return f"BC{int(match.group(0)):02d}"
 
def find_bam_path(row, base_path: Path) -> Path | None:
    """
    Constructs and verifies the existence of a BAM file path for a given row.
    Checks two possible directory structures.
    """
    if pd.isna(row.SequenceRun) or not row.StandardizedBarcode:
        return None
 
    barcode_num_only = row.StandardizedBarcode.replace("BC", "")
    # Path structure 1: e.g., .../mapped/barcode16_mapped.bam
    path1 = base_path / str(row.SequenceRun) / "mapped" / f"barcode{barcode_num_only}_mapped.bam"
    # Path structure 2: e.g., .../result/mapped/BC16_mapped.bam
    path2 = base_path / str(row.SequenceRun) / "result" / "mapped" / f"{row.StandardizedBarcode}_mapped.bam"
 
    if path2.is_file():
        return path2
    if path1.is_file():
        return path1
    return None
 
def main():
    parser = argparse.ArgumentParser(
        description="Find and copy BAM mapping files corresponding to a list of GISAID EPI_ISL IDs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--gisaid-ids", required=True, help="Path to the input file containing GISAID IDs (one per line).")
    parser.add_argument("--sample-sheet", required=True, help="Path to the sam.tsv file that contains the 'gisaid accession id' column.")
    parser.add_argument("--master-excel", required=True, help="Path to the master metadata .xls or .xlsx file.")
    parser.add_argument("--config", required=True, help="Path to the configuration INI file.")
    parser.add_argument("--output-dir", required=True, help="Directory where the found BAM files will be copied.")
    args = parser.parse_args()
 
    # 1. Setup and Load Config
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    config = configparser.ConfigParser()
    config.read(args.config)
    try:
        sequence_data_root = Path(config['PATHS']['sequence_data_root'])
    except KeyError:
        print(f" FATAL: 'sequence_data_root' key not found under [PATHS] in '{args.config}'.", file=sys.stderr)
        sys.exit(1)
 
    # 2. Load and Prepare Data
    print(f"Reading target GISAID IDs from '{args.gisaid_ids}'...")
    try:
        target_ids_df = pd.read_csv(args.gisaid_ids, header=None)
        target_ids = set(target_ids_df.iloc[:, 0].str.strip().dropna().unique())
        print(f"Found {len(target_ids)} unique target IDs.")
    except FileNotFoundError:
        print(f" FATAL: GISAID ID file not found at '{args.gisaid_ids}'.", file=sys.stderr)
        sys.exit(1)
 
    print(f"Reading sample sheet from '{args.sample_sheet}'...")
    try:
        sam_df = pd.read_csv(args.sample_sheet, sep='\t')
        target_samples_df = sam_df[sam_df['gisaid accession id'].isin(target_ids)].copy()
        target_samples_df.rename(columns={'gisaid accession id': 'gisaid_accession_id'}, inplace=True)
        print(f"Found {len(target_samples_df)} matching samples in the sample sheet.")
    except (FileNotFoundError, KeyError):
        print(f" FATAL: Sample sheet not found or is missing the 'gisaid accession id' column at '{args.sample_sheet}'.", file=sys.stderr)
        sys.exit(1)
    print(f"Reading master metadata from '{args.master_excel}'...")
    try:
        master_df = pd.read_excel(args.master_excel, engine='openpyxl' if args.master_excel.endswith('xlsx') else 'xlrd')
        master_df['StandardizedBarcode'] = master_df['Barcode'].apply(standardize_barcode)
    except FileNotFoundError:
        print(f" FATAL: Master Excel file not found at '{args.master_excel}'.", file=sys.stderr)
        sys.exit(1)
    # 3. Merge DataFrames to link GISAID ID to Run Info
    print("Linking GISAID IDs to sequencing run information...")
    merged_df = pd.merge(
        target_samples_df[['alias', 'gisaid_accession_id']],
        master_df[['SampleID', 'SequenceRun', 'StandardizedBarcode']],
        left_on='alias',
        right_on='SampleID',
        how='inner'
    )
    # 4. Find and Copy BAM Files
    print(f"Searching for BAM files under '{sequence_data_root}' and copying to '{output_dir}'...")
    found_count = 0
    for row in merged_df.itertuples():
        source_path = find_bam_path(row, sequence_data_root)
        if source_path:
            # Create a more descriptive filename for the copied file
            dest_filename = f"{row.gisaid_accession_id}_{source_path.name}"
            dest_path = output_dir / dest_filename
            try:
                shutil.copy2(source_path, dest_path)
                print(f"  [COPIED]   {row.gisaid_accession_id} -> {dest_path}")
                found_count += 1
            except Exception as e:
                print(f"  [ERROR]    Failed to copy {source_path}: {e}", file=sys.stderr)
        else:
            print(f"  [MISSING]  {row.gisaid_accession_id}: Could not find BAM for Run {row.SequenceRun}, Barcode {row.StandardizedBarcode}", file=sys.stderr)
 
    # 5. Final Summary
    print("\n--- Summary ---")
    print(f"Successfully found and copied {found_count} out of {len(merged_df)} target BAM files.")
    print(f"Files are located in: {output_dir.resolve()}")
    print("✔ Done.")
 
if __name__ == "__main__":
    main()