#!/usr/bin/env python3
import sys
import argparse
import configparser
from lxml import etree

"""
Usage: create_exp_xml.py [-h] --input INPUT --output OUTPUT --config CONFIG

Converts an `exp.tsv` file into an ENA-compliant `exp.xml` file using a configuration file.

This script takes a tab-separated file with experiment metadata and transforms it
into the XML format required for an ENA EXPERIMENT object submission.

Note: The corresponding SAMPLE objects must be successfully submitted to ENA
before submitting the experiment and run metadata.
"""

def create_experiment_attribute(parent, tag, value):
    """
    Helper function to create a standard EXPERIMENT_ATTRIBUTE XML element.
    """
    experiment_attr = etree.SubElement(parent, "EXPERIMENT_ATTRIBUTE")
    etree.SubElement(experiment_attr, "TAG").text = tag
    etree.SubElement(experiment_attr, "VALUE").text = value

def tsv_to_experiment_xml(tsv_input_file, xml_output_file, config):
    """
    Converts a TSV file of experiment data into an ENA-formatted XML file.
    """
    # 1. Get config values
    try:
        cfg_exp = config['ENA_EXPERIMENT_DESIGN']
        title = cfg_exp['title']
        study_accession = cfg_exp['accession']
        platform = cfg_exp['platform']
        instrument_model = cfg_exp['instrument_model']
        library_strategy = cfg_exp['library_strategy']
        library_source = cfg_exp['library_source']
        library_selection = cfg_exp['library_selection']
        library_layout = cfg_exp['library_layout']
        
        cfg_meta = config['ENA_METADATA']
        center_name = cfg_meta['center_name']
    except KeyError as e:
        print(f"Error: Missing key {e} in [ENA_EXPERIMENT_DESIGN] or [ENA_METADATA] section of your config file.", file=sys.stderr)
        sys.exit(1)

    # 2. Process TSV and build XML
    root = etree.Element("EXPERIMENT_SET")
    line_counter = 0
    skipped_lines = 0
    written_lines = 0

    try:
        with open(tsv_input_file, 'r') as f:
            # We read the config once, but loop through the TSV to create an experiment for each line
            for line in f:
                line_counter += 1
                line = line.strip()

                if not line:
                    continue

                if "expalias" in line.lower():
                    skipped_lines += 1
                    continue

                columns = line.split('\t')
                if len(columns) != 4:
                    print(f"Warning: Line {line_counter} is malformed (expected 4 columns, found {len(columns)}). Skipping.", file=sys.stderr)
                    skipped_lines += 1
                    continue

                sam_alias, exp_alias, _, _ = columns # Only need first two columns here

                if not all([sam_alias, exp_alias]):
                    print(f"Warning: Line {line_counter} contains empty required values (samAlias or expAlias). Skipping.", file=sys.stderr)
                    skipped_lines += 1
                    continue

                # Create EXPERIMENT element and its children
                experiment = etree.SubElement(root, "EXPERIMENT", {
                    "alias": exp_alias,
                    "center_name": center_name
                })

                etree.SubElement(experiment, "TITLE").text = title

                study_ref = etree.SubElement(experiment, "STUDY_REF")
                study_ref.attrib["accession"] = study_accession

                design = etree.SubElement(experiment, "DESIGN")
                etree.SubElement(design, "DESIGN_DESCRIPTION").text = ""
                sample_descriptor = etree.SubElement(design, "SAMPLE_DESCRIPTOR")
                sample_descriptor.attrib["refname"] = sam_alias

                library_descriptor = etree.SubElement(design, "LIBRARY_DESCRIPTOR")
                etree.SubElement(library_descriptor, "LIBRARY_NAME").text = "" # Library Name is usually sample-specific, can be left blank
                etree.SubElement(library_descriptor, "LIBRARY_STRATEGY").text = library_strategy
                etree.SubElement(library_descriptor, "LIBRARY_SOURCE").text = library_source
                etree.SubElement(library_descriptor, "LIBRARY_SELECTION").text = library_selection

                # This is the key change for library_layout
                layout_element = etree.SubElement(library_descriptor, "LIBRARY_LAYOUT")
                etree.SubElement(layout_element, library_layout) # Use the config value as the TAG name

                # This is the key change for platform
                platform_element = etree.SubElement(experiment, "PLATFORM")
                # Use the config value as the TAG name for the platform-specific block
                platform_block = etree.SubElement(platform_element, platform)
                etree.SubElement(platform_block, "INSTRUMENT_MODEL").text = instrument_model

                experiment_attributes = etree.SubElement(experiment, "EXPERIMENT_ATTRIBUTES")
                create_experiment_attribute(experiment_attributes, "library preparation date", "not collected")

                written_lines += 1

    except FileNotFoundError:
        print(f"Error: The input file '{tsv_input_file}' was not found.", file=sys.stderr)
        sys.exit(1)

    # 3. Write XML to file
    tree = etree.ElementTree(root)
    with open(xml_output_file, 'wb') as file:
        tree.write(file, pretty_print=True, xml_declaration=True, encoding="UTF-8")

    return written_lines

# Main execution block
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Convert an `exp.tsv` file into an ENA-compliant `exp.xml` file using a configuration file.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to the input TSV file containing experiment metadata."
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

    # --- Read Config and pass it to the function ---
    config = configparser.ConfigParser()
    config.read(args.config)

    try:
        count = tsv_to_experiment_xml(args.input, args.output, config)
        print(f"Successfully created {count} experiment objects in '{args.output}'.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)