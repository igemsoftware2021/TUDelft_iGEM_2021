import csv
import re


def read_ngs_references(path: str) -> dict:
    """Function reads the csv file containing the reference sequences with their
    respective names and returns a dictionary with the sequence coupled to the name."""
    ngs_references = {}
    with open(path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            ngs_references[row["sequence"]] = row["name"]
    return ngs_references


def create_ngs_references_patterns(ngs_references_dict: dict) -> dict:
    patterns_dict = {}
    for seq in ngs_references_dict:
        patterns_dict[re.compile(seq)] = ngs_references_dict[seq]
    return patterns_dict


def reference_seq(sequence: str, ngs_references_dict: dict):
    """Function returns the name of the reference sequence or the string 'NULL'"""
    for seq_pattern in ngs_references_dict:
        match = seq_pattern.search(sequence)
        if match:
            # True
            return ngs_references_dict[seq_pattern]
    # False
    return "NULL"


def determine_clvd_prefix(sequence: str, clvd_prefix_info: dict = {"seq": "ACAAAACAAAAC", "name": "Z"}, unclvd_prefix_info: dict = {"seq": "AAACAAACAAA", "name": "W"}, additional_prefix=None) -> tuple:
    """
    Function returns whether a sequence matches with a cleaved or uncleaved prefix sequences and returns
    the name of the prefix.\n
    args:\n
    sequence: (str) sequence you want to check.\n
    clvd_prefix_info: (dict) dictionary containing the sequence and name of the cleaved prefix.\n
    unclvd_prefix_info: (dict) dictionary containing the sequence and name of the uncleaved prefix.\n
    additional_prefix: (dict) dictionary containing the sequence and name of an optional additional prefix.\n
    \n
    return values:\n
    clvd_prefix_info: (int) value that tells you whether it is the cleaved prefix. 1 = Yes, 0 = No, 2 = Neither clvd or unclvd prefix.\n
    prefix_name: (str) name of the found prefix.
    """

    prefix_match = re.search(clvd_prefix_info["seq"], sequence)
    if bool(prefix_match):
        prefix = clvd_prefix_info["seq"]
        prefix_name = clvd_prefix_info["name"]
        clvd_prefix_info = 1
    else:
        prefix_match = re.search(unclvd_prefix_info["seq"], sequence)
        if bool(prefix_match):
            prefix = unclvd_prefix_info["seq"]
            prefix_name = unclvd_prefix_info["name"]
            clvd_prefix_info = 0
        else:
            prefix = None
            prefix_name = "NULL"
            clvd_prefix_info = 2  # TODO decide whether to implement additional sequence argument
    return (clvd_prefix_info, prefix_name, prefix)


def retrieve_barcode(sequence: str, prefix: str) -> str:
    if prefix is not None:
        prefix_match = re.search(prefix, sequence)
        return sequence[:prefix_match.start()]
    else:
        return "NULL"


def cleanup_sequence(sequence: str, prefix: str, suffix: str) -> str:
    if prefix is not None and suffix is not None:
        prefix_match = re.search(prefix, sequence)
        suffix_match = re.search(suffix, sequence)
        if not bool(suffix_match):
            return "NULL"
        return sequence[prefix_match.end():suffix_match.start()]
    return "NULL"
