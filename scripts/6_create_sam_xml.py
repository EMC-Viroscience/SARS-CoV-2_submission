#!/usr/bin/env python3
import sys
import argparse
import configparser
import pandas as pd
from lxml import etree

"""
Usage:
    python scripts/create_sam_xml.py \
        --input "input_tsv/sam.tsv" \
        --output "output_xml/sam.xml" \
        --config "config/ena_config.ini"

Description:
    Converts a standardized 11-column TSV file into an ENA-compliant SAMPLE XML file.
    This script is data-driven, mapping TSV columns directly to ENA sample attributes.
    Static metadata (like taxon ID, description) is loaded from an external config file.
"""

def create_sample_attribute(parent, tag, value, units=None):
    """Helper function to create a standard SAMPLE_ATTRIBUTE XML element."""
    sample_attr = etree.SubElement(parent, "SAMPLE_ATTRIBUTE")
    etree.SubElement(sample_attr, "TAG").text = str(tag)
    etree.SubElement(sample_attr, "VALUE").text = str(value)
    if units:
        etree.SubElement(sample_attr, "UNITS").text = str(units)

def tsv_to_sample_xml(tsv_input_file, xml_output_file, config_file):
    """
    Converts a TSV file of sample data into an ENA-formatted XML file using
    parameters from a configuration file.
    """
    # 1. Load Configuration
    config = configparser.ConfigParser()
    config.read(config_file)
    try:
        cfg = config['ENA_METADATA']
    except KeyError:
        print(f"Error: [ENA_METADATA] section not found in '{config_file}'", file=sys.stderr)
        sys.exit(1)

    # 2. Load and Validate Input TSV
    try:
        df = pd.read_csv(tsv_input_file, sep='\t', dtype=str)
        # The 11 mandatory ENA attribute fields from your new TSV
        expected_columns = [
            "alias", "collection date", "geographic location (country and/or sea)",
            "host common name", "host subject id", "collecting institution", "isolate",
            "host scientific name", "host health state", "host sex", "collector name", "gisaid accession id"
        ]
        if not all(col in df.columns for col in expected_columns):
            missing = set(expected_columns) - set(df.columns)
            print(f"Error: Input TSV is missing required columns: {', '.join(missing)}", file=sys.stderr)
            sys.exit(1)
    except FileNotFoundError:
        print(f"Error: The input file '{tsv_input_file}' was not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading or parsing TSV file: {e}", file=sys.stderr)
        sys.exit(1)

    # 3. Build XML tree
    root = etree.Element("SAMPLE_SET")
    
    for index, row in df.iterrows():
        # Create SAMPLE element
        sample = etree.SubElement(root, "SAMPLE", {
            "alias": row["alias"],
            "center_name": cfg.get('center_name', 'DEFAULT_CENTER') # Use a fallback
        })

        # Add Static Metadata from Config
        etree.SubElement(sample, "TITLE").text = cfg.get('title', 'Default Sample Title')
        
        sample_name = etree.SubElement(sample, "SAMPLE_NAME")
        etree.SubElement(sample_name, "TAXON_ID").text = cfg.get('taxon_id')
        etree.SubElement(sample_name, "SCIENTIFIC_NAME").text = cfg.get('scientific_name')
        etree.SubElement(sample_name, "COMMON_NAME").text = cfg.get('common_name')
        
        etree.SubElement(sample, "DESCRIPTION").text = cfg.get('description')

        # Add Dynamic Sample Attributes from TSV
        sample_attributes = etree.SubElement(sample, "SAMPLE_ATTRIBUTES")
        
        # Dynamically create attributes from the 11 specified columns
        for col_name in expected_columns:
            # The alias is an attribute of the SAMPLE tag, not a SAMPLE_ATTRIBUTE
            if col_name != "alias":
                create_sample_attribute(sample_attributes, col_name, row[col_name])
        
        # Add Hardcoded Sample Attributes
        create_sample_attribute(sample_attributes, "ENA-CHECKLIST", "ERC000033")

    # 4. Write XML to file
    tree = etree.ElementTree(root)
    with open(xml_output_file, 'wb') as file:
        tree.write(file, pretty_print=True, xml_declaration=True, encoding="UTF-8")

    return len(df)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Convert a standardized sample TSV into an ENA-compliant XML file using a config file.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to the input TSV file containing sample metadata."
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

    try:
        count = tsv_to_sample_xml(args.input, args.output, args.config)
        print(f"Successfully created {count} sample objects in '{args.output}'.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)