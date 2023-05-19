import multiprocessing
from syntreesizer import images
import re
from pathlib import Path
from collections import defaultdict
from multiprocessing import Pool
from sys import argv

if len(argv) != 2:
    print(f"ERROR: Expected one argument, got {len(argv) - 1}.\n")
    print("Usage:\n"
          "python3 -c ProcessImages.py <dir>\n\n"
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
    normal_captures.append(values[0])
    gt_captures.append(values[1])

del grouped_files

if __name__ == "__main__":
    multiprocessing.freeze_support()

    with Pool() as pool:
        results = pool.imap_unordered(images.Image.contained_processing_chain, zip(normal_captures, gt_captures))

        for _ in results:
            pass
