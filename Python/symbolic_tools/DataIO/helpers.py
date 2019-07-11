import os
import re
import numpy as np


def find_oldest_iteration(folder, extension='.pvti'):
    iterations = []
    for root, dirs, files in os.walk(folder):
        for file in files:
            if file.endswith(extension):
                # it_counter = re.findall(r"(\d{8})", file)
                it_counter = re.findall(r"(\d{8})", file)[0]
                iterations.append(it_counter)

    # int_iterations = [int(i) for i in iterations]
    # oldest = max(int_iterations)
    # idx = int_iterations.index(oldest)
    oldest = max(iterations)
    return oldest


def calc_mse(anal, num):
    return np.sum((anal - num) * (anal - num)) / len(anal)


def calc_L2(anal, num):
    # Eq. 4.57
    return np.sqrt(np.sum((anal - num) * (anal - num)) / np.sum(anal * anal))
