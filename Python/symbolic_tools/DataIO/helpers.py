import os
import re
import numpy as np


def get_r_from_xy(x, y, x0=0, y0=0):
    r = np.sqrt(pow(x0 - x, 2) + pow(y0 - y, 2))
    return r


def eat_dots_for_texmaker(value):
    s_value = str(value)
    s_value = re.sub(r"\.", 'o', s_value)
    return s_value


def strip_folder_name(some_folder):
    some_folder = some_folder.rstrip('0')  # remove all trailing zeros
    some_folder = some_folder.rstrip('.')
    return some_folder


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
    if not iterations:
        raise FileNotFoundError(f'Check the path: \n {folder}')

    oldest = max(iterations)
    return oldest


def get_vti_from_iteration(folder, iteration, extension='.vti'):
    pattern = f"VTK_P00_{iteration}{extension}"
    matched_files = []
    for root, dirs, files in os.walk(folder):
        for file in files:
            if file.endswith(pattern):
                matched_files.append(file)

    if not matched_files:
        raise FileNotFoundError(f'Check the path: \n {folder}')

    if len(matched_files) > 1:
        raise Exception("More than 1 matching file")

    return matched_files[0]


def calc_mse(anal, num):
    return np.sum((anal - num) * (anal - num)) / len(anal)


def calc_L2(anal, num):
    # Eq. 4.57
    return np.sqrt(np.sum((anal - num) * (anal - num)) / np.sum(anal * anal))
