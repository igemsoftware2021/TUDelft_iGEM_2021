import yaml
import re

with open("sequence_analysis_pipeline/config.yaml", "r") as rf:
    try:
        yaml_info = yaml.safe_load(rf)
        driver_round_pattern = yaml_info["info_patterns"]["driver_round"]
        selection_pat = yaml_info["info_patterns"]["selection"]
        ligand_present_pat = yaml_info["info_patterns"]["ligand_present"]
        print(driver_round_pattern)
    except yaml.YAMLError as exc:
        print(exc)

test = "test"

print(driver_round_pattern)
test_pattern = r"D([0-9]+)"
filename = "S1_D80_L0_read_count.txt"
print(re.search(test_pattern, filename).group())
print(re.search(driver_round_pattern, filename).group())

print(type(re.search(selection_pat, filename).group()))
print(re.search(ligand_present_pat, filename).group())
