import csv
import re


def complement_reverse_seq(sequence: str):
    """Function creates a complementary sequence and reverses the sequence."""
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    # Create a list containing the complement bases
    complement_seq = [complement[base] for base in sequence]
    # Reverse the list
    complement_seq = complement_seq[::-1]
    complement_seq = ''.join(complement_seq)
    return complement_seq


def read_ngs_references(path: str) -> dict:
    """Function reads the csv file containing the reference sequences with their
    respective names and returns a dictionary with the sequence coupled to the name."""
    ngs_references = {}
    with open(path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            ngs_references[row["sequence"]] = row["name"]
    return ngs_references


def complement_reverse_references(ngs_references_dict: dict) -> dict:
    patterns_dict = {}
    for seq in ngs_references_dict:
        seq = complement_reverse_seq(seq)
        seq = clean_reference_sequence(seq, "AAAAAGAAA")
        patterns_dict[seq] = ngs_references_dict[seq]
    return patterns_dict


def cleanup_reference_sequences_dict(ngs_references_dict: dict):
    cleaned_ngs_references = {}
    for reference_seq in ngs_references_dict:
        reversed_seq = complement_reverse_seq(reference_seq)
        cleaned_seq = clean_reference_sequence(reversed_seq, "AAAAAGAAA")
        print(cleaned_seq)
        cleaned_ngs_references[cleaned_seq] = ngs_references_dict[reference_seq]
    return cleaned_ngs_references


def create_ngs_references_patterns(ngs_references_dict: dict) -> dict:
    patterns_dict = {}
    for seq in ngs_references_dict:
        seq = clean_reference_sequence(seq, "AAAAAGAAA")
        print(seq)
        patterns_dict[re.compile(seq)] = ngs_references_dict[seq]
    return patterns_dict


def clean_reference_sequence(sequence, suffix: str):
    prefix_seq, _ = determine_prefix(sequence)
    cleaned_seq = cleanup_sequence(sequence, prefix_seq, suffix)
    return cleaned_seq


def reference_seq(sequence: str, ngs_references_dict: dict):
    """Function returns the name of the reference sequence or the string 'NULL'"""
    for reference_seq in ngs_references_dict:
        match = re.search(reference_seq, sequence)
        if match:
            # True
            return ngs_references_dict[reference_seq]
    # False
    return "NULL"


def determine_prefix(sequence: str, prefix_info={"ACAAAACAAAAC": "Z", "AAACAAACAAA": "W", "CTTTTCCGTATATCTCGCCAG": "A"}):
    for prefix_seq in prefix_info:
        prefix_match = re.search(prefix_seq, sequence)
        if bool(prefix_match):
            prefix_name = prefix_info[prefix_seq]
            return prefix_seq, prefix_name
    return None, None


def determine_clvd_prefix(prefix_name, clvd_prefix_name="Z", unclvd_prefix_name="W"):
    if prefix_name.lower() == clvd_prefix_name.lower():
        clvd_prefix = 1
    elif prefix_name.lower() == unclvd_prefix_name.lower():
        clvd_prefix = 0
    else:
        clvd_prefix = 2
    return clvd_prefix


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
