#!/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from CHECLabPy.core.io import DL1Reader
from muon_fitter import analyze_muon_event
from astropy.table import Table

if __name__ == '__main__':
    
    dl1_path = "/lustrehome/pillera/CTA_CHEC/mc/Muons/run1_dl1.h5"
    reader = DL1Reader(dl1_path)
    
    iev = reader.select_column('iev')
    # evnr = 0
    # image = reader[evnr]['photons'].values
    # a = reader.read("pointing")
    # alt = a['altitude_raw'][evnr]
    # azi = a['azimuth_raw'][evnr]
    # muon_evt = analyze_muon_event(image,reader.mapping,alt,azi)

    index = []
    par = {'Index': [],
            'Size': [],
            'RingComp': []}
    for i in np.unique(iev):
        if i%10 == 0:
            print(i)
        image = reader[i]['photons'].values
        a = reader.read("pointing")
        alt = a['altitude_raw'][i]
        azi = a['azimuth_raw'][i]
        muon_evt = analyze_muon_event(image,reader.mapping,alt,azi)
        if muon_evt[0] or muon_evt[1] is not None:
            par['Index'].append(i)
            par['Size'].append(np.sum(image))
            if muon_evt[1]['RingComp'] is not None:
                par['RingComp'].append(muon_evt[1]['RingComp'].ring_completeness)
            else:
                par['RingComp'].append(-1)

    tab = index.to_table()
    tab.write("/lustrehome/CTA_CHEC/mc/Muon_search_result.fits",format="fits")
    print("Found "+str(len(index)+" muons"))


