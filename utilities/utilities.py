import csv
from typing import List
import re
import os


def sorted_human(iterable: List[str]) -> List[str]:
    """ Sort the given iterable in the way that humans expect."""
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(iterable, key=alphanum_key)


def is_dir_exists(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)


def write_dict_to_csv(dict, csv_path):
    with open(csv_path, 'w') as f:  # You will need 'wb' mode in Python 2.x
        w = csv.DictWriter(f, dict.keys())
        w.writeheader()
        w.writerow(dict)
