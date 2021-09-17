import csv
import time
import re


def write_temperature_csv(time: list, temperature: list, path: str = f"temperature-{time.time()}.csv"):
    """Function to write time and temperature data to a csv file at a specific location.

    Parameters
    ----------
    time: list
        The timepoints at which the temperature measurements were taken. Timepoint at
        index 1 corresponds to the temperature value at index 1 in the temperature list
    temperature: list
        The temperature values stored at certain timepoints. Temperature value at index 1
        corresponds to the timepoint at index 1 in the time list
    path: str
        The path to where the csv file should be written (default f"temperature-{time.time()}.csv")

    Raises
    ------
    ValueError
        If the time list and temperature are not of the same length
    """
    if len(time) != len(temperature):
        raise ValueError(
            "The time list and temperature list are not of the same length")

    with open(file=path, mode="w", encoding="utf-8") as csv_wf:
        fieldnames = ["time", "temperature"]

        writer = csv.DictWriter(csv_wf, fieldnames=fieldnames)

        # Write the fieldnames to the top of the file
        writer.writeheader()

        for i in range(len(time)):
            writer.writerow({"time": time[i], "temperature": temperature[i]})


def read_temperature_csv(path: str):
    """Function to read time and temperature data from a csv file at a specific location.

    Parameters
    ----------
    path: str
        The path to the csv file from which the data should be read.

    Returns
    -------
    time: list
        The timepoints at which the temperature measurements were taken. Timepoint at
        index 1 corresponds to the temperature value at index 1 in the temperature list.
    temperature: list
        The temperature values stored at certain timepoints. Temperature value at index 1
        corresponds to the timepoint at index 1 in the time list.

    Raises
    ------
    ValueError
        If the csv file does not contain the proper fieldnames, namely
        "time" and "temperature"
    """
    with open(file=path, mode="r", encoding="utf-8") as csv_rf:
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


def write_absorbance_csv(time: list, absorbance: list, path: str = f"absorbance-{int(time.time())}.csv"):
    if len(time) != len(absorbance):
        raise ValueError(
            "The variable time and variable absorbance should contain the same amount of lists")

    num_channels = len(absorbance)

    max_len = max(len(absorbance_list) for absorbance_list in absorbance)

    # The absorbance is a list of lists, where every list contains the values of one sensor
    with open(file=path, mode="w", encoding="utf-8") as csv_wf:
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


def read_absorbance_csv(path: str):
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
