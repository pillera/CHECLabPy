import numpy as np
import matplotlib.pyplot as plt
from CHECLabPy.core.io import DL1Reader
from muon_fitter import analyze_muon_event

if __name__ == '__main__':
    
    dl1_path = "/lustrehome/pillera/CTA_CHEC/mc/Muons/run1_dl1.h5"
    reader = DL1Reader(dl1_path)
    evnr = 0
    image = reader[evnr]['photons'].values
    a = reader.read("pointing")
    alt = a['altitude_raw'][evnr]
    azi = a['azimuth_raw'][evnr]
    muon_evt = analyze_muon_event(image,reader.mapping)
