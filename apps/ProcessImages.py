import multiprocessing
from syntreesizer import Images
import re
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool
from sys import argv

def filter_train(x: Path):
    return "gt" not in str(x)

def filter_valid(x: Path):
    return "gt" in str(x)

if len(argv) != 2:
    print(f"ERROR: Expected one argument, got {len(argv) - 1}.\n")
    print("Usage:\n"
          "python3 ProcessImages.py <dir>\n\n"
          "<dir>\tAbsolute path to directory containing exported images from CityEngine (both training and ground "
          "truth passes)")
    exit(1)
else:
    # "C:/Users/flori/Desktop/Uni/Bachelorarbeit/city-engine/images/experiment_1"
    base_path = Path(argv[1])

files = base_path.glob("*.png")
grouped_files = defaultdict(list)
normal_captures, gt_captures = list(), list()

for file in files:
    key = re.split(r"(\d+)", file.stem)
    grouped_files[key[1]].append(file)

for key, values in grouped_files.items():
    if len(values) != 2:
        print("Non matching images omitted\n")
        continue
    normal_captures.append(*filter(filter_valid, values))
    gt_captures.append(*filter(filter_valid, values))

del grouped_files

if __name__ == "__main__":
    multiprocessing.freeze_support()

    with Pool() as pool:
        results = pool.imap_unordered(Images.Image.contained_processing_chain, zip(normal_captures, gt_captures))

        for _ in results:
            pass
