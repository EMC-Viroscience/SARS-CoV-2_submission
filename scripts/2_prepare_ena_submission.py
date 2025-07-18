#!/usr/bin/env python3
import argparse
import configparser
import gzip
import hashlib
import os
import shutil
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import re
from tqdm import tqdm

"""
================================================================================
ENA Submission Preparation and FASTQ Processing Script
================================================================================

Usage:
    python scripts/2_prepare_ena_submission.py \
        --input-excel "/path/to/your/master_sheet.xls" \
        --output-dir "input_tsv" \
        --config "config/sars.ini" \
        --fastq-output-dir "raw_fastq" \
        --workers 8

Description:
    This script orchestrates the entire pre-submission process for ENA.

    Part 1: FASTQ Processing
    1. Reads metadata from the input Excel file.
    2. For each sample, finds the correct raw FASTQ file from two possible locations.
    3. Copies and compresses each found FASTQ file into the --fastq-output-dir.
    4. Names the new file using the scheme: {SequenceRun}_{StandardizedBarcode}.fastq.gz.
    5. Calculates the MD5 checksum of the *newly created gzipped file*.

    Part 2: TSV Generation
    1. Generates 'sam.tsv' with 11 mandatory columns for sample registration.
    2. Generates 'exp.tsv' and 'run.tsv', using the new FASTQ filenames and their calculated checksums.
    3. Saves all TSV files to the specified --output-dir.

Dependencies:
    pip install pandas xlrd openpyxl tqdm
"""

# Helper Functions

def standardize_barcode(barcode: str) -> str:
    """Standardizes a barcode string to the format 'BC' followed by two digits."""
    if pd.isna(barcode):
        return ""
    match = re.search(r'\d+', str(barcode))
    if not match:
        return ""
    return f"BC{int(match.group(0)):02d}"

def compress_and_get_md5(source_path: Path, dest_path: Path) -> tuple[str, str]:
    """
    Reads a source file, compresses it with gzip, writes it to the destination,
    and then calculates and returns the MD5 checksum of the *compressed* file.
    """
    try:
        # Step 1: Compress the file
        with open(source_path, 'rb') as f_in, gzip.open(dest_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

        # Step 2: Calculate MD5 on the newly created gzipped file
        md5 = hashlib.md5()
        with open(dest_path, 'rb') as f_gz:
            # Read in chunks to handle large files efficiently
            while chunk := f_gz.read(8192):
                md5.update(chunk)
        
        return dest_path.name, md5.hexdigest()

    except Exception as e:
        # If any part of the process fails, ensure the incomplete destination file is removed.
        if dest_path.exists():
            dest_path.unlink()
        raise IOError(f"Failed to process {source_path}: {e}") from e

# Main Logic Functions

def process_all_fastq_files(df: pd.DataFrame, fastq_output_dir: Path, max_workers: int, sequence_data_root: str) -> pd.DataFrame:
    """Finds, compresses, copies, and checksums all FASTQ files in parallel."""
    print(f"\n--- Part 1: Processing FASTQ Files ---")
    
    base_path = Path(sequence_data_root)
    tasks = []
    
    for row in df.itertuples():
        if pd.isna(row.SequenceRun) or not row.StandardizedBarcode:
            continue
        
        barcode_num_only = row.StandardizedBarcode.replace("BC", "")
        path1 = base_path / str(row.SequenceRun) / "filtered" / f"barcode{barcode_num_only}_filtered.fastq"
        path2 = base_path / str(row.SequenceRun) / "result" / "filtered" / f"{row.StandardizedBarcode}_filtered.fastq"

        source_path = path2 if path2.is_file() else (path1 if path1.is_file() else None)
        
        if source_path:
            dest_name = f"{row.SequenceRun}_{row.StandardizedBarcode}.fastq.gz"
            dest_path = fastq_output_dir / dest_name
            tasks.append({'source': source_path, 'dest': dest_path, 'index': row.Index})
        else:
            print(f"    WARNING: No FASTQ file found for Run: {row.SequenceRun}, Barcode: {row.Barcode}", file=sys.stderr)

    if not tasks:
        print("\nNo valid FASTQ files found to process. Cannot generate exp.tsv or run.tsv.", file=sys.stderr)
        return pd.DataFrame()

    print(f"Found {len(tasks)} FASTQ files. Compressing and copying with {max_workers} worker(s)...")
    
    results = []
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_task = {
            executor.submit(compress_and_get_md5, task['source'], task['dest']): task 
            for task in tasks
        }
        
        progress = tqdm(as_completed(future_to_task), total=len(tasks), desc="Processing FASTQ")
        for future in progress:
            task = future_to_task[future]
            try:
                filename, checksum = future.result()
                results.append({'index': task['index'], 'fastq_filename': filename, 'md5': checksum})
            except Exception as exc:
                tqdm.write(f"  ? FAILED to process {task['source']}: {exc}", file=sys.stderr)

    return pd.DataFrame(results).set_index('index')

def generate_ena_sheets(df: pd.DataFrame, output_dir: Path):
    """Generates sam.tsv, exp.tsv, and run.tsv from the processed DataFrame."""
    print("\n--- Part 2: Generating ENA Submission Sheets ---")
    
    # sam.tsv
    column_map = {
        'SampleID': 'alias', 'SamplingDate': 'collection date',
        'PatientCountry': 'geographic location (country and/or sea)',
        'Host': 'host common name', 'PatientCode': 'host subject id',
        'SequenceID_Internal': 'isolate'
    }
    sam_df = df.rename(columns=column_map)
    
    sam_df['host scientific name'] = np.where(sam_df['host common name'] == 'Human', 'Homo sapiens', 'not provided')
    sam_df['host health state'] = 'not provided'
    sam_df['host sex'] = 'not provided'
    sam_df['collector name'] = collector_name
    sam_df['host subject id'] = 'not provided'
    sam_df['collecting institution'] = 'not provided'

    dates = pd.to_datetime(sam_df['collection date'], errors='coerce')
    sam_df['collection date'] = np.where(
        dates.notna(),
        dates.dt.strftime('%Y-%m-%d'),
        'not provided'
    )
    
    sam_df = sam_df[list(column_map.values()) + ['host scientific name', 'host health state', 'host sex', 'collecting institution', 'collector name']]
    sam_df.to_csv(output_dir / 'sam.tsv', sep='\t', index=False, na_rep='not provided')
    print(f" 'sam.tsv' created with {len(sam_df)} records.")

    # exp.tsv & run.tsv
    if 'fastq_filename' not in df.columns:
        print("No FASTQ data available, skipping exp.tsv and run.tsv.", file=sys.stderr)
        return
        
    # Filter out rows where FASTQ processing failed
    processed_df = df.dropna(subset=['fastq_filename', 'md5'])
    
    # Create a unique experiment alias
    processed_df['expAlias'] = processed_df['SequenceRun'].astype(str) + '_' + processed_df['StandardizedBarcode'].astype(str)

    # exp.tsv
    exp_df = pd.DataFrame({
        'samAlias': processed_df['SampleID'],
        'expAlias': processed_df['expAlias'],
        'filename': processed_df['fastq_filename'],
        'md5': processed_df['md5']
    })
    exp_df.to_csv(output_dir / 'exp.tsv', sep='\t', index=False)
    print(f" 'exp.tsv' created with {len(exp_df)} records.")

    # run.tsv
    run_df = pd.DataFrame({
        'runAlias': processed_df['SampleID'],
        'expAlias': processed_df['expAlias'],
        'filename': processed_df['fastq_filename'],
        'md5': processed_df['md5']
    })
    run_df.to_csv(output_dir / 'run.tsv', sep='\t', index=False)
    print(f" 'run.tsv' created with {len(run_df)} records.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Prepare ENA submission files (TSVs and gzipped FASTQs) from a master metadata Excel sheet.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input-excel", required=True, help="Path to the master '.xls' or '.xlsx' metadata file.")
    parser.add_argument("--output-dir", default="input_tsv", help="Directory to save the generated TSV files.")
    parser.add_argument("--config", required=True, help="Path to the configuration INI file.")
    parser.add_argument("--fastq-output-dir", default="raw_fastq", help="Directory to save the compressed FASTQ files.")
    parser.add_argument("--workers", type=int, default=4, help="Number of parallel processes for FASTQ compression.")
    
    args = parser.parse_args()

    # --- Main Orchestration ---
    output_path = Path(args.output_dir)
    fastq_path = Path(args.fastq_output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    fastq_path.mkdir(parents=True, exist_ok=True)

    # --- Read Config ---
    config = configparser.ConfigParser()
    config.read(args.config)
    try:
        sequence_data_root = config['PATHS']['sequence_data_root']
        collector_name = config['ENA_METADATA']['collector_name']
        print(f"Using sequence data root from config: '{sequence_data_root}'")
    except KeyError:
        print(f" FATAL: 'sequence_data_root' or 'collector_name' not found in '{args.config}'.", file=sys.stderr)
        sys.exit(1)

    # Read Excel

    try:
        print(f"Reading main metadata file: '{args.input_excel}'")
        main_df = pd.read_excel(args.input_excel, engine='openpyxl' if args.input_excel.endswith('xlsx') else 'xlrd')
        main_df['StandardizedBarcode'] = main_df['Barcode'].apply(standardize_barcode)
    except Exception as e:
        print(f"? FATAL: Could not read Excel file: {e}. Aborting.", file=sys.stderr)
        sys.exit(1)

    # Step 1: Process FASTQs and get back a dataframe with results
    fastq_results_df = process_all_fastq_files(main_df, fastq_path, args.workers, sequence_data_root)

    # Step 2: Merge FASTQ results back into the main dataframe
    if not fastq_results_df.empty:
        main_df = main_df.merge(fastq_results_df, left_index=True, right_index=True, how='left')

    # Step 3: Generate all TSV files from the combined dataframe
    generate_ena_sheets(main_df, output_path)

    print("\n All tasks complete.")