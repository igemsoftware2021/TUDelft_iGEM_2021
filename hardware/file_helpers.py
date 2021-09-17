import csv
import time
import re


def write_temperature_csv(time, temperature, path=f"temperature-{time.time()}.csv"):
    """Function to write temperature data to a csv file"""
    if len(time) != len(temperature):
        raise ValueError(
            "The time list and temperature list are not of the same length")

    with open(path, mode="w") as csv_wf:
        fieldnames = ["time", "temperature"]

        writer = csv.DictWriter(csv_wf, fieldnames=fieldnames)

        # Write the fieldnames to the top of the file
        writer.writeheader()

        for i in range(len(time)):
            writer.writerow({"time": time[i], "temperature": temperature[i]})


def read_temperature_csv(path):
    """Function to read the temperature data from a csv file"""
    with open(path, "r") as csv_rf:
        reader = csv.DictReader(csv_rf)
        fieldnames = reader.fieldnames

        if "time" not in fieldnames and "temperature" not in fieldnames:
            raise ValueError(
                f"Check the fieldnames of {path}, time and temperature are missing")

        time = []
        temperature = []

        for row in reader:
            time.append(float(row["time"]))
            temperature.append(float(row["temperature"]))
    return time, temperature


def write_absorbance_csv(time: list, absorbance: list, path=f"absorbance-{int(time.time())}.csv"):
    if len(time) != len(absorbance):
        raise ValueError(
            "The variable time and variable absorbance should contain the same amount of lists")

    num_channels = len(absorbance)

    max_len = max(len(absorbance_list) for absorbance_list in absorbance)

    # The absorbance is a list of lists, where every list contains the values of one sensor
    with open(path, mode="w", encoding="utf-8") as csv_wf:
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
    with open(file=path, mode="r", encoding="utf-8") as csv_rf:

        reader = csv.DictReader(csv_rf)

        headers = reader.fieldnames

        # Find the total number of channels of which there is absorbance in the csv file
        max_chan_number = 0
        for fieldname in headers:

            # This regex matches an integer in the name, use re.search
            # not re.match because re.match only checks from the first character
            # onwards
            result = re.search(r"[\d]+", fieldname)
            chan_number = int(result.group(0))
            if chan_number > max_chan_number:
                max_chan_number = chan_number

        num_channels = max_chan_number + 1  # Do +1, because the chan_number starts at 0

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
