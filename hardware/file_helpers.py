import csv
import time
import numpy as np
import re


def write_temperature_csv(time, temperature, path=f"temperature-{time.time()}.csv"):
    if len(time) != len(temperature):
        raise ValueError(
            "The time list and temperature list are not of the same length")


def write_absorbance_csv(time, absorbance, path=f"absorbance-{time.time()}.csv"):
    if len(time) != len(absorbance):
        raise ValueError(
            "The variable time and variable absorbance should contain the same amount of lists")

    num_channels = len(absorbance)

    max_len = max(len(absorbance_list) for absorbance_list in absorbance)

    # The absorbance is a list of lists, where every list contains the values of one sensor
    with open(path, mode="w") as csv_wf:
        dict_to_write = dict()

        fieldnames = []
        for i in range(num_channels):
            fieldnames.append(f"time-ch{i}")
            fieldnames.append(f"absorbance-ch{i}")

            dict_to_write[f"time-ch{i}"] = None
            dict_to_write[f"absorbance-ch{i}"] = None

        writer = csv.DictWriter(csv_wf, fieldnames=fieldnames)

        # Write the fieldnames to the top of the file
        writer.writeheader()

        for j in range(max_len):
            for i in range(num_channels):
                if j < len(absorbance[i]):
                    dict_to_write[f"time-ch{i}"] = time[i][j]
                    dict_to_write[f"absorbance-ch{i}"] = absorbance[i][j]
                else:
                    dict_to_write[f"time-ch{i}"] = "-"
                    dict_to_write[f"absorbance-ch{i}"] = "-"

            # Write the information to the corresponding csv file
            writer.writerow(dict_to_write)


def read_absorbance_csv(path):
    """Function reads an absorbance file and returns a list of timepoints and a list of absorbance points"""
    with open(path, mode="r") as csv_rf:

        reader = csv.DictReader(path)
        headers = reader.fieldnames

        # Find the total number of channels of which there is absorbance in the csv file
        num_channels = 0
        for fieldname in headers:
            # This regex matches an integer
            result = re.match(r"[\d]+", fieldname)
            number = result.group(0)
            if number > num_channels:
                num_channels = number

        # Create the lists in which to store the timepoints and absorbance
        time = []
        absorbance = []

        # Add the same number of lists as there are channels
        for _ in range(num_channels):
            time.append([])
            absorbance.append([])

        for row in reader:
            for i in range(num_channels):
                if row[f"time-ch{i}"] != "-":
                    time[i].append(float(row[f"time-ch{i}"]))
                    absorbance[i].append(float(row[f"absorbance-ch{i}"]))

    return time, absorbance
