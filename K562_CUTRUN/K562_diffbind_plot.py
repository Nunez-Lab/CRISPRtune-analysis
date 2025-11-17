import polars as pl
import matplotlib.pyplot as plt
import matplotlib as mpl
import altair as alt
import numpy as np

import mygene
import pyperclip

import glob
import os
import sys
import math

mpl.rcParams.update({
    "pdf.fonttype": 42,           # TrueType = editable text
    "font.family": "Liberation Sans", # or Arial, Helvetica
    "text.usetex": False          # Absolutely must be False
})

CD55_ID = "ENSG00000196352"

# INPUT_DIR = "/Users/p.coli/Desktop/CRISPRtune/RPE_CUTRUN/diffbind/diffbind_broad/DiffBind_broad/edgeR/"
# OUTPUT_DIR = "/Users/p.coli/Desktop/CRISPRtune/RPE_CUTRUN/Diff_bind_volcanoplots/"

INPUT_DIR = sys.argv[1]
OUTPUT_DIR = sys.argv[2]

DESeq_dict = {}

order = ["H3K9me3","H3K27ac","H3K4me3"] #hardcoded
print(f"Note: This script hardcoded for {order}")

for i in order:
    for file in glob.glob(f"{INPUT_DIR}edgeR_{i}*_NT_CD55.csv"):
        plot_name = (os.path.basename(file)).replace("edgeR_", "").replace(".csv", "")
        df = pl.read_csv(file, null_values=["NA"])
        df = df.with_columns(
            ((pl.col("Fold") * -1).alias("neg_log2_Fold")),
            ((pl.col("FDR").log10() * -1).alias("neg_log10_fdr"))
        )
        if i != "H3K9me3":
            df = df.with_columns(
                pl.when(
                    (pl.col("seqnames") == "chr1") & 
                    # (pl.col("start") > 207_321_000) &
                    # (pl.col("start") < 207_340_000)
                    (pl.col("start") > 207_321_678-15_000) &
                    (pl.col("end") < 207_321_678+15_000)
                ).then(pl.lit("red")).otherwise(pl.lit("grey")).alias("color")
            )
        else:
            df = df.with_columns(
                pl.when(
                    (pl.col("seqnames") == "chr1") & 
                    # (pl.col("start") > 207_321_000) &
                    # (pl.col("start") < 207_340_000)
                    (pl.col("start") > 207_321_678-16_000) &
                    (pl.col("end") < 207_321_678+5000)
                ).then(pl.lit("red")).otherwise(pl.lit("grey")).alias("color")
            )#.filter(pl.col("CD55") == "red")
        
        xticks = np.arange(-7, 8, 1)
        
        
        if i == "H3K9me3":
            yticks = np.arange(0, 7, 1)
            
        elif i == "H3K27ac":
            yticks = np.arange(0, 14, 2)
            
        elif i == "H3K4me3":
            if "CRISPRi" in plot_name:
                yticks = np.arange(0, 251, 50)
            else:
                yticks = np.arange(0, 12, 2)

        fig, ax = plt.subplots(1, 1, figsize=(4, 4))

        for (c,), g in df.group_by("color"):
            if c == "red":
                alpha = 1
                zorder = 10
            else:
                alpha = 0.1
                zorder = 5
            ax.scatter(g["neg_log2_Fold"], g["neg_log10_fdr"], color=c, alpha=alpha, zorder=zorder,rasterized=True,)

        ax.set_xticks(xticks)
        ax.set_yticks(yticks)

        ax.set_xlabel("-log2FoldChange")
        ax.set_ylabel("-log10(padj)")

        plt.axvline(x=-1, color='gray', linestyle='--', label='Cutoff')
        plt.axvline(x=1, color='gray', linestyle='--', label='Cutoff')
        plt.axhline(y=1.3010299956639813, color='gray', linestyle='--', label='Cutoff')
        ax.set_title(f"{plot_name} -Log2FC vs -Log10Padj", fontsize=12, color="black")
        mpl.rcParams['svg.fonttype'] = 'none'
        
        fig.savefig(OUTPUT_DIR + "png/" +f"{plot_name}.png", format="png", dpi = 500)
        fig.savefig(OUTPUT_DIR + "svg/" +f"{plot_name}.svg", format="svg", dpi = 300)

