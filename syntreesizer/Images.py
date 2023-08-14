from pathlib import Path
from typing import Optional, List

import numpy
import numpy as np
import png


class Image:
    def __init__(self, fp: str) -> None:
        self.fp: Path = Path(fp)
        self.png: png.Reader = png.Reader(fp)
        self.direct = None
        self.processed_data: Optional[np.ndarray] = None

    def prepare_segmentation(self, out: Optional[str] = None):
        direct_representation = self.png.asDirect()

        self.processed_data = np.vstack(map(np.uint8, direct_representation[2]))

        self.processed_data[self.processed_data == 255] = 0

        self.processed_data[(0 < self.processed_data) & (self.processed_data < 255)] = 1

        self.processed_data = self.processed_data[:, ::3].copy()

        if out:
            w = png.Writer(direct_representation[0], direct_representation[1], greyscale=True, bitdepth=1)

            with open(out, "wb") as f:
                w.write(f, self.processed_data)

    def subset_band(self, band_index: int):
        self.direct = self.png.asDirect()

        image_2d = np.vstack(map(np.uint8, self.direct[2]))

        image_3d = np.reshape(image_2d, (self.direct[1], self.direct[0], self.direct[3]["planes"]))

        self.processed_data = image_3d[:, :, band_index].copy()

    def ground_truth_from_diff(self, other: "Image", out: Optional[str] = None) -> Optional[np.ndarray]:
        out_data = self.processed_data - other.processed_data

        out_data[out_data > 0] = 1

        if out:
            w = png.Writer(self.png.width, self.png.height, greyscale=True, bitdepth=1)

            with open(out, "wb") as f:
                w.write(f, out_data)
        else:
            return out_data

    @staticmethod
    def contained_processing_chain(vals: List[Path]):
        im = Image(str(vals[0]))
        gt = Image(str(vals[1]))

        im.subset_band(1)
        gt.subset_band(1)

        im.ground_truth_from_diff(gt, str(vals[1].parent / vals[1].name.replace("diff_", "")))
