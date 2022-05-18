from starmapper import analysis
from .starmapper import *
import starmapper.styles

sample = ["Sample",
          "create_code_button",
          "get_color_list"]

arrays = ["array_qc",
          "array_ap",
          "array_skyline",
          "array_ss",
          "array_heatmap",
          "array_mol",
          "array_circle",
          "array_linreg"]

plots = ["AP",
         "Circle",
         "DistHist",
         "HeatMaP",
         "LinReg",
         "Mol",
         "Plots",
         "QC",
         "Skyline",
         "SM",
         "SS"]

data = ["Annotation",
        "CT",
        "Data",
        "DP",
        "IJ",
        "Log",
        "PDB",
        "Profile"]

analysis = ["LowSS",
            "LogCompare"]

__all__ = sample + arrays + plots + data + analysis
