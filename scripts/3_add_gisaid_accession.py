import pandas as pd
import argparse
from pathlib import Path
import json
import re
import sys

def extract_unique_key(text: str) -> str:
    """
    Extracts the unique identifier (e.g., 'ZH-EMC-8906') from a string.
    """
    if pd.isna(text):
        return None
    # Updated regex to be more general for different lab codes
    match = re.search(r'[A-Z]{2,}-[A-Z]{2,}-\d+', text)
    if match:
        return match.group(0)
    return None

def parse_gisaid_log(log_path: Path) -> dict:
    """
    Parses a covCLI log file to extract a mapping of unique keys to GISAID EPI ISL IDs.
    Handles both successful new submissions and errors for already existing samples.
    """
    print(f"Parsing GISAID log file: '{log_path}'...")
    accession_map = {}
    
    with open(log_path, 'r') as f:
        for line in f:
            try:
                log_entry = json.loads(line.strip())
            except json.JSONDecodeError:
                continue
                
            code = log_entry.get('code')
            msg = log_entry.get('msg', '')
            
            if code == 'epi_isl_id':
                parts = msg.split(';')
                if len(parts) >= 2:
                    virus_name, accession_id = parts[0].strip(), parts[1].strip()
                    key = extract_unique_key(virus_name)
                    if key and accession_id:
                        accession_map[key] = accession_id

            elif code == 'validation_error' and 'already exists' in msg:
                virus_name_match = re.search(r"existing_virus_name: (.*?);", msg)
                accession_match = re.search(r"EPI_ISL_\d+", msg)
                
                if virus_name_match and accession_match:
                    virus_name = virus_name_match.group(1).strip()
                    key = extract_unique_key(virus_name)
                    accession_id = accession_match.group(0)
                    if key and accession_id:
                        accession_map[key] = accession_id

    if not accession_map:
        print("    - Warning: No successful or existing accession IDs found in the log file.")
    else:
        print(f" Successfully parsed {len(accession_map)} total accession IDs from the log.")
        
    return accession_map

def sync_ena_sheets_with_gisaid(input_dir: str, log_file: str, keep_all: bool):
    """
    Adds GISAID accession IDs to sam.tsv and filters all three ENA TSV files
    (sam, exp, run) to keep them synchronized.
    
    Args:
        input_dir (str): Path to the directory containing sam.tsv, exp.tsv, and run.tsv.
        log_file (str): Path to the covCLI log file.
        keep_all (bool): If True, keeps all samples, filling missing IDs with 'not provided'.
                         If False, only samples with a found GISAID ID are kept.
    """
    input_path = Path(input_dir)
    sam_path = input_path / "sam.tsv"
    exp_path = input_path / "exp.tsv"
    run_path = input_path / "run.tsv"
    
    accession_map = parse_gisaid_log(Path(log_file))
    
    # 1. Read all three TSV files
    print(f"\nReading ENA sheets from directory: '{input_path}'...")
    try:
        sam_df = pd.read_csv(sam_path, sep='\t')
        exp_df = pd.read_csv(exp_path, sep='\t')
        run_df = pd.read_csv(run_path, sep='\t')
        original_count = len(sam_df)
    except FileNotFoundError as e:
        print(f" Error: Missing required TSV file: {e.filename}", file=sys.stderr)
        return
        
    # 2. Add GISAID IDs to the sample sheet
    print("Extracting unique keys from the 'isolate' column...")
    sam_df['join_key'] = sam_df['isolate'].apply(extract_unique_key)
    
    print("Mapping GISAID accession IDs to samples...")
    sam_df['gisaid accession id'] = sam_df['join_key'].map(accession_map)
    sam_df = sam_df.drop(columns=['join_key'])

    # 3. Filter all three DataFrames if needed
    if keep_all:
        print("`--keep-all` specified: Keeping all samples and experiments.")
        sam_df['gisaid accession id'] = sam_df['gisaid accession id'].fillna('not provided')
    else:
        print("Default behavior: Filtering all sheets to keep only samples with a GISAID accession ID.")
        
        # Identify which samples to keep based on successful mapping
        samples_to_keep = sam_df[sam_df['gisaid accession id'].notna()].copy()
        
        # Get the list of sample aliases to keep
        aliases_to_keep = set(samples_to_keep['alias'])
        
        # Filter all three dataframes
        sam_df_filtered = samples_to_keep
        exp_df_filtered = exp_df[exp_df['samAlias'].isin(aliases_to_keep)].copy()
        run_df_filtered = run_df[run_df['runAlias'].isin(aliases_to_keep)].copy()
        
        print(f"  - Samples: Removed {original_count - len(sam_df_filtered)} records. Kept {len(sam_df_filtered)}.")
        print(f"  - Experiments: Removed {len(exp_df) - len(exp_df_filtered)} records. Kept {len(exp_df_filtered)}.")
        print(f"  - Runs: Removed {len(run_df) - len(run_df_filtered)} records. Kept {len(run_df_filtered)}.")
        
        # Replace original dataframes with filtered ones
        sam_df, exp_df, run_df = sam_df_filtered, exp_df_filtered, run_df_filtered

    # 4. Overwrite the original TSV files
    print(f"\nOverwriting TSV files in '{input_path}'...")
    sam_df.to_csv(sam_path, sep='\t', index=False)
    exp_df.to_csv(exp_path, sep='\t', index=False)
    run_df.to_csv(run_path, sep='\t', index=False)
    
    num_with_id = (sam_df['gisaid accession id'] != 'not provided').sum()
    print(f" Success! Files are synchronized. Final sample count: {len(sam_df)}, with {num_with_id} GISAID IDs.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Adds GISAID accession IDs to sam.tsv and filters sam.tsv, exp.tsv, and run.tsv to keep them synchronized.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--tsv-dir",
        default="input_tsv",
        help="Path to the directory containing sam.tsv, exp.tsv, and run.tsv."
    )
    
    parser.add_argument(
        "--log-file",
        required=True,
        help="Path to the covCLI log file containing submission results."
    )
    
    parser.add_argument(
        "--keep-all",
        action='store_true',
        help="If specified, keep all samples in the output, even those without an accession ID (fills with 'not provided')."
    )

    args = parser.parse_args()
    
    sync_ena_sheets_with_gisaid(args.tsv_dir, args.log_file, args.keep_all)