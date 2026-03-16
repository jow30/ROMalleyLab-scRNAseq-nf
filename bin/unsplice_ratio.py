import sys
import numpy as np
import loompy
import pandas as pd

def compute_unsplice_ratio_to_txt(loom_file, output_txt):
    with loompy.connect(loom_file) as ds:
        # Extract spliced and unspliced matrices
        spliced = ds.layer['spliced'][:]
        unspliced = ds.layer['unspliced'][:]

        # Get cell barcodes (column attributes)
        if 'CellID' in ds.ca:
            cell_barcodes = ds.ca['CellID']
        elif 'cell_names' in ds.ca:
            cell_barcodes = ds.ca['cell_names']
        else:
            cell_barcodes = np.array([f"cell_{i}" for i in range(spliced.shape[1])])

        # Compute totals per cell
        total_unspliced = np.sum(unspliced, axis=0)
        total_spliced = np.sum(spliced, axis=0)
        total = total_unspliced + total_spliced

        # Compute ratio safely
        ratio = np.divide(
            total_unspliced,
            total,
            out=np.zeros_like(total_unspliced, dtype=float),
            where=total != 0
        )

        # Make a DataFrame
        df = pd.DataFrame({
            "cell_barcode": cell_barcodes,
            "total_unspliced": total_unspliced,
            "total_spliced": total_spliced,
            "unsplice_ratio": ratio
        })

        # Save as tab-separated file
        df.to_csv(output_txt, sep="\t", index=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python unsplice_ratio.py <loom_file> <output_txt>")
        sys.exit(1)

    compute_unsplice_ratio_to_txt(sys.argv[1], sys.argv[2])


