import yaml
import regex


def retrieve_dataset_name(yaml_file: str):
    """Function that retrieves the name of the dataset."""
    with open(yaml_file, "r") as rf:
        try:
            yaml_info = yaml.safe_load(rf)

            dataset_name = yaml_info["dataset"]

        except yaml.YAMLError as exc:
            print(exc)
    return dataset_name


def retrieve_minimum_reads(yaml_file: str) -> int:
    """Function that retrieves the minimum amount of reads needed before allowing it in the database."""
    with open(yaml_file, "r") as rf:
        try:
            yaml_info = yaml.safe_load(rf)

            minimum_reads = int(yaml_info["minimum_reads"])

        except yaml.YAMLError as exc:
            print(exc)
    return minimum_reads


def retrieve_compiled_patterns(yaml_file: str, pattern: str = "prefix", ligand_present: bool = True) -> dict:
    """Function reads a yaml file and retrieves either prefixes or suffixes for a certain condition for
    the ligand. If the ligand present is set to True, then it will retrieve the prefixes or suffixes
    for which there was ligand present and otherwise for which no ligand was present.\n
    \n
    args:\n
    yaml_file: (str) path to yaml file.\n
    pattern: (str) 'prefix' or 'suffix', determines whether you want to retrieve the prefix or suffix patterns.\n
    ligand_present: (bool) indicates for what ligand condition you want to retrieve the patterns.\n
    \n
    returns:\n
    compiled_patterns: (dict) dictionary where the key, value pair is as follows: regex compiled pattern, name of the sequence.
    """
    # First do some instance checks
    if not isinstance(pattern, str):
        raise TypeError(
            "pattern keyword argument should be of type string")
    if not isinstance(ligand_present, bool):
        raise TypeError(
            "ligand_present keyword argument should be of type boolean")

    if pattern == "prefix":
        pattern_key = "prefixes"
    elif pattern == "suffix":
        pattern_key = "suffixes"
    else:
        raise ValueError(
            "pattern keyword argument should be 'prefix' or 'suffix'")

    if ligand_present:
        ligand_key = "ligand_positive"
    else:
        ligand_key = "ligand_negative"

    with open(yaml_file, "r") as rf:
        try:
            yaml_info = yaml.safe_load(rf)

            # Retrieve the cleave states for a certain ligand condition
            # This can be either: cleaved/uncleaved/other
            cleave_states = yaml_info[ligand_key][pattern_key]

            # Create regex patterns from the prefix sequences
            compiled_patterns = {}
            for cleave_state in cleave_states:
                sequence = yaml_info[ligand_key][pattern_key][cleave_state]["sequence"]
                sequence_name = yaml_info[ligand_key][pattern_key][cleave_state]["name"]
                max_error = yaml_info[ligand_key][pattern_key][cleave_state]["max_error"]
                compiled_patterns[regex.compile(
                    fr"(?e)({sequence}){{e<={max_error}}}")] = sequence_name

        except yaml.YAMLError as exc:
            print(exc)
    return compiled_patterns


def retrieve_compiled_reference_patterns(yaml_file: str, pattern: str = "prefix") -> dict:
    """Function reads a yaml file and retrieves either prefixes or suffixes for a certain condition for
    the ligand. If the ligand present is set to True, then it will retrieve the prefixes or suffixes
    for which there was ligand present and otherwise for which no ligand was present.\n
    \n
    args:\n
    yaml_file: (str) path to yaml file.\n
    pattern: (str) 'prefix' or 'suffix', determines whether you want to retrieve the prefix or suffix patterns.\n
    ligand_present: (bool) indicates for what ligand condition you want to retrieve the patterns.\n
    \n
    returns:\n
    compiled_patterns: (dict) dictionary where the key, value pair is as follows: regex compiled pattern, name of the sequence.
    """
    # First do some instance checks
    if not isinstance(pattern, str):
        raise TypeError(
            "pattern keyword argument should be of type string")

    if pattern == "prefix":
        pattern_key = "prefixes"
    elif pattern == "suffix":
        pattern_key = "suffixes"
    else:
        raise ValueError(
            "pattern keyword argument should be 'prefix' or 'suffix'")
    with open(yaml_file, "r") as rf:
        try:
            yaml_info = yaml.safe_load(rf)

            # Retrieve the cleave states for a certain ligand condition
            # This can be either: cleaved/uncleaved/other
            cleave_states = yaml_info["references"][pattern_key]

            # Create regex patterns from the prefix sequences
            compiled_patterns = {}
            for cleave_state in cleave_states:
                sequence = yaml_info["references"][pattern_key][cleave_state]["sequence"]
                sequence_name = yaml_info["references"][pattern_key][cleave_state]["name"]
                max_error = yaml_info["references"][pattern_key][cleave_state]["max_error"]
                compiled_patterns[regex.compile(
                    fr"(?e)({sequence}){{e<={max_error}}}")] = sequence_name

        except yaml.YAMLError as exc:
            print(exc)
    return compiled_patterns


def retrieve_prefix_name(yaml_file: str, cleaved: bool = True, ligand_present: bool = True) -> str:

    # First do instance checks
    if not isinstance(cleaved, bool):
        raise TypeError("cleaved keyword argument should be of type boolean")

    if not isinstance(ligand_present, bool):
        raise TypeError(
            "ligand_present keyword argument should be of type boolean")

    if ligand_present:
        ligand_key = "ligand_positive"
    else:
        ligand_key = "ligand_negative"

    if cleaved:
        cleave_key = "cleaved"
    else:
        cleave_key = "uncleaved"

    with open(yaml_file, "r") as rf:
        try:
            yaml_info = yaml.safe_load(rf)

            # Store the name of the prefix
            prefix_name = yaml_info[ligand_key]["prefixes"][cleave_key]["name"]

        except yaml.YAMLError as exc:
            print(exc)
    return prefix_name


def retrieve_compiled_info_patterns(yaml_file: str) -> tuple:
    with open(yaml_file, "r") as rf:
        try:
            yaml_info = yaml.safe_load(rf)

            driver_round_pattern = regex.compile(
                yaml_info["info_patterns"]["driver_round"])
            selection_pattern = regex.compile(
                yaml_info["info_patterns"]["selection"])
            ligand_present_pattern = regex.compile(
                yaml_info["info_patterns"]["ligand_present"])

        except yaml.YAMLError as exc:
            print(exc)
    return driver_round_pattern, selection_pattern, ligand_present_pattern
