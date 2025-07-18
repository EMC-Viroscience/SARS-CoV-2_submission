import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import configparser
import sys

def generate_gisaid_files(input_file: str, output_dir: str, config: configparser.ConfigParser):
    """
    Reads an Excel file to generate a GISAID-formatted metadata.csv
    and an all_sequences.fasta file with custom formatting.

    Args:
        input_file (str): Path to the input Excel file.
        output_dir (str): Directory to save the output files.
        config (configparser.ConfigParser): Parsed configuration object.
    """
    # 1. Read config values
    try:
        cfg_gisaid = config['GISAID']
        submitter_email = cfg_gisaid['submitter']
        author_list_string = cfg_gisaid['author_list_string']
    except KeyError as e:
        print(f" Error: Missing key {e} in [GISAID] section of your config file.", file=sys.stderr)
        sys.exit(1)

    # 2. Read and process input data
    print(f"Reading data from '{input_file}'...")
    try:
        df = pd.read_excel(input_file)
    except FileNotFoundError:
        print(f" Error: Input file not found at '{input_file}'")
        return
    except Exception as e:
        print(f" Error reading Excel file: {e}")
        return

    # 3. Create metadata.csv
    print("Generating 'metadata.csv'...")

    # Start of formatting logic
    dates = pd.to_datetime(df['SamplingDate'], errors='coerce')
    year_str = np.where(dates.notna(), dates.dt.strftime('%Y'), '')
    virus_name_with_year = df['SequenceID_gisaid'].astype(str) + '/' + year_str
    formatted_location = df['ContinentCountryProv'].str.replace('\\', ' / ', regex=False)
    # End of formatting logic

    # Create the GISAID DataFrame using the newly formatted Series
    df_gisaid = pd.DataFrame({
        'covv_virus_name': virus_name_with_year,
        'covv_collection_date': dates.dt.strftime('%Y-%m-%d'),
        'covv_location': formatted_location,
        'covv_host': df['Host'],
        'covv_sampling_strategy': df['SampleQuestion'],
        'covv_specimen': df['SampleType'],
        'covv_coverage': df['Coverage']
    })

    # Add hardcoded and config-driven values
    df_gisaid['submitter'] = submitter_email  # From config
    df_gisaid['covv_authors'] = author_list_string # From config
    
    df_gisaid['fn'] = "all_sequences.fasta"
    df_gisaid['covv_type'] = "betacoronavirus"
    df_gisaid['covv_passage'] = "Original"
    df_gisaid['covv_add_location'] = "unknown"
    df_gisaid['covv_add_host_info'] = "unknown"
    df_gisaid['covv_gender'] = "unknown"
    df_gisaid['covv_patient_age'] = "unknown"
    df_gisaid['covv_patient_status'] = "unknown"
    df_gisaid['covv_outbreak'] = "unknown"
    df_gisaid['covv_last_vaccinated'] = "unknown"
    df_gisaid['covv_treatment'] = "unknown"
    df_gisaid['covv_seq_technology'] = "Nanopore"
    df_gisaid['covv_assembly_method'] = "unknown"
    df_gisaid['covv_orig_lab_addr'] = "Dr. Molewaterplein 3015GE Rotterdam Netherlands"
    df_gisaid['covv_provider_sample_id'] = "unknown"
    df_gisaid['covv_subm_lab_addr'] = "Dr. Molewaterplein 3015GE Rotterdam Netherlands"
    df_gisaid['covv_consortium'] = "Dutch National COVID-19 Response Team"
    df_gisaid['covv_orig_lab'] = "Dutch COVID-19 response team"
    df_gisaid['covv_subm_lab'] = "Erasmus Medical Center"
    df_gisaid['covv_subm_sample_id'] = "unknown"

    # Fill any remaining missing values (NaN/NaT) with "unknown"
    df_gisaid = df_gisaid.fillna("unknown")

    # Reorder columns to the final required format
    final_column_order = [
        'submitter', 'fn', 'covv_virus_name', 'covv_type', 'covv_passage',
        'covv_collection_date', 'covv_location', 'covv_add_location', 'covv_host',
        'covv_add_host_info', 'covv_sampling_strategy', 'covv_gender', 'covv_patient_age',
        'covv_patient_status', 'covv_specimen', 'covv_outbreak', 'covv_last_vaccinated',
        'covv_treatment', 'covv_seq_technology', 'covv_assembly_method', 'covv_coverage',
        'covv_orig_lab', 'covv_orig_lab_addr', 'covv_provider_sample_id',
        'covv_subm_lab', 'covv_subm_lab_addr', 'covv_subm_sample_id',
        'covv_consortium', 'covv_authors'
    ]
    df_gisaid = df_gisaid[final_column_order]

    # Prepare output directory and save the CSV
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    metadata_output_file = output_path / "metadata.csv"
    df_gisaid.to_csv(metadata_output_file, index=False)
    print(f"? 'metadata.csv' successfully created at '{metadata_output_file}'")


    # 4. Create all_sequences.fasta
    print("Generating 'all_sequences.fasta'...")
    fasta_output_file = output_path / "all_sequences.fasta"
    sequences_written = 0
    
    df_gisaid_for_fasta = df_gisaid.copy()
    df_gisaid_for_fasta['Sequence'] = df['Sequence']

    with open(fasta_output_file, 'w') as f_out:
        for row in df_gisaid_for_fasta.itertuples():
            header = row.covv_virus_name
            sequence = row.Sequence
            
            if pd.notna(header) and pd.notna(sequence):
                f_out.write(f">{header}\n")
                f_out.write(f"{sequence}\n")
                sequences_written += 1
            else:
                print(f"    - Skipping row {row.Index}: missing 'SequenceID_gisaid' or 'Sequence'.")
    
    print(f"? 'all_sequences.fasta' successfully created at '{fasta_output_file}'")
    print(f"    Total sequences written: {sequences_written}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Prepare GISAID submission files (metadata.csv and all_sequences.fasta) from an Excel sheet.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--input-file",
        required=True,
        help="Path to the input Excel metadata file."
    )
    
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Path to the output directory where submission files will be stored."
    )

    parser.add_argument(
        "--config",
        required=True,
        help="Path to the configuration INI file."
    )

    args = parser.parse_args()
    
    # --- Read Config and pass it to the function ---
    config = configparser.ConfigParser()
    config.read(args.config)

    generate_gisaid_files(args.input_file, args.output_dir, config)