#!/usr/bin/env python3
import sys
import argparse
import configparser
from lxml import etree

"""
Usage: create_run_xml.py [-h] --input INPUT --output OUTPUT --config CONFIG

Converts a `run.tsv` file into an ENA-compliant `run.xml` file using a configuration file.

This script takes a tab-separated file with run metadata (linking samples to
experiments and data files) and transforms it into the XML format required
for an ENA RUN object submission.
"""

def tsv_to_run_xml(tsv_input_file, xml_output_file):
    """
    Converts a TSV file of run data into an ENA-formatted XML file.
    """

    # 1. Get config values
    try:
        center_name = config['ENA_METADATA']['center_name']
    except KeyError:
        print(f" Error: Missing key 'center_name' in [ENA_METADATA] section of your config file.", file=sys.stderr)
        sys.exit(1)

    root = etree.Element("RUN_SET")
    line_counter = 0
    skipped_lines = 0
    written_lines = 0

    try:
        with open(tsv_input_file, 'r') as f:
            for line in f:
                line_counter += 1
                line = line.strip()

                if not line:
                    continue

                # Skip header row (case-insensitive check)
                if "runalias" in line.lower() or "samalias" in line.lower():
                    skipped_lines += 1
                    continue

                columns = line.split('\t')
                if len(columns) != 4:
                    print(f"Warning: Line {line_counter} is malformed (expected 4 columns, found {len(columns)}). Skipping.", file=sys.stderr)
                    skipped_lines += 1
                    continue

                run_alias, exp_alias, gz_file, md5 = columns

                if not all([run_alias, exp_alias, gz_file, md5]):
                    print(f"Warning: Line {line_counter} contains empty required values. Skipping.", file=sys.stderr)
                    skipped_lines += 1
                    continue

                # Create RUN element and its children
                run = etree.SubElement(root, "RUN", {
                    "alias": run_alias,
                    "center_name": center_name
                })

                experiment_ref = etree.SubElement(run, "EXPERIMENT_REF")
                experiment_ref.attrib["refname"] = exp_alias

                data_block = etree.SubElement(run, "DATA_BLOCK")
                files = etree.SubElement(data_block, "FILES")

                etree.SubElement(files, "FILE", {
                    "filename": gz_file,
                    "filetype": "fastq",
                    "checksum_method": "MD5",
                    "checksum": md5
                })
                
                written_lines += 1

    except FileNotFoundError:
        print(f"Error: The input file '{tsv_input_file}' was not found.", file=sys.stderr)
        sys.exit(1)

    # Format and write the XML to the output file
    tree = etree.ElementTree(root)
    with open(xml_output_file, 'wb') as file:
        tree.write(file, pretty_print=True, xml_declaration=True, encoding="UTF-8")

    return written_lines

# Main execution block
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Convert a `run.tsv` file into an ENA-compliant `run.xml` file using a configuration file.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to the input TSV file containing run metadata."
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path for the output XML file."
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Path to the configuration INI file."
    )

    args = parser.parse_args()

    # Read Config and pass it to the function
    config = configparser.ConfigParser()
    config.read(args.config)

    try:
        # Pass the parsed arguments to the conversion function
        count = tsv_to_run_xml(args.input, args.output)
        print(f"Successfully created {count} run objects in '{args.output}'.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)