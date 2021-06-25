import regex

driver_round_pattern = regex.compile(r"(?<=D)([0-9]+)")

test_string = "T1_D80_L0"

print(driver_round_pattern.search(test_string).group())
