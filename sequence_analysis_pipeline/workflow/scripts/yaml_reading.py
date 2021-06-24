import yaml
import regex


def retrieve_compiled_patterns(yaml_file: str, pattern: str = "prefix", ligand_present: bool = True):
    """Function reads a yaml file and retrieves either prefixes or suffixes for a certain condition for
    the ligand. If the ligand present is set to True, then it will retrieve the prefixes or suffixes
    for which there was ligand present and otherwise for which no ligand was present.\n
    \n
    args:\n

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


def retrieve_cleaved_prefix_name(yaml_file: str, cleaved: bool = True, ligand_present: bool = True):

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


def retrieve_compiled_info_patterns(yaml_file: str):
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
