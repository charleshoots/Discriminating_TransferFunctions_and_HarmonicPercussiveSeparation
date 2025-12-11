# Author: Charles Hoots
# This code was developed as part of my PhD research in the
# Department of Earth Sciences, University of Hawai‘i at Mānoa.
# Unless otherwise noted, the code is my own original work.
# External libraries and standard research software packages are used as cited.

### Every single line of the following code was written directly by Charles Hoots 
### in service of his PhD dissertation research at the University of Hawaii Manoa 
### with zero external assistance or contributions from others.
### ---
### Use of any data, analysis, or codes contained or produced within 
### is open-source, covered under the MIT License:
### https://github.com/charleshoots/ATACR_HPS_Comp/blob/main/LICENSE.txt
##################################################################################

import sys;from pathlib import Path;sys.path.append(str(Path(__file__).parent.parent.parent))
import os,sys;from source.imports import *;from source.modules import *

# dirs value
dirs=io.dir_libraries()
# cat value
cat = catalog.copy()
# cat = cat[cat.StaName.isin(['7D.G17B', '7D.G25B'])]
outfold = dirs.Data/'Analysis'/'Metrics'/'noise'

# comps value
comps = ['c12', 'c1Z', 'c1P', 'c2Z', 'c2P', 'cZP']
# PSDs value
PSDs = lambda d: AttribDict({k[-1]:d.power.__dict__[k] for k in list(d.power.__dict__.keys())})
# coherences value
coherences = lambda d:{''.join(sorted(c[1:],reverse=True)):lt.math._calc_coherence(d.cross.__dict__[c],d.power.__dict__[f'c{c[1]}{c[1]}'],d.power.__dict__[f'c{c[2]}{c[2]}']) for c in comps}
# admittances value
admittances = lambda d:{''.join(sorted(c[1:],reverse=True)):lt.math._calc_admittance(d.cross.__dict__[c],d.power.__dict__[f'c{c[2]}{c[2]}']) for c in comps}
phases = lambda d:{''.join(sorted(c[1:],reverse=True)):lt.math._calc_phase(d.cross.__dict__[c]) for c in comps}
staavg = lambda d:AttribDict({i:m(d) for i,m in zip(['coh','adm','ph','psd'],[coherences,admittances,phases,PSDs])})
dayavg = lambda fold,d:AttribDict({i:[m(day) for day in [load_pickle(fold/f) for f in d.day_files[d.gooddays]]] for i,m in zip(['coh','adm','ph','psd'],[coherences,admittances,phases,PSDs])})

issues = []
for si,stanm in enumerate(cat.StaName):
    try:
        print(f'{si+1}/{len(cat)} | {stanm}')
        d = load_pickle(list((dirs.SpectraAvg/stanm).glob('*.avg_sta.pkl'))[0])
        fold = dirs.Spectra/ stanm
        sta = staavg(d)
        day = dayavg(fold,d)
        cpa={}
        for m in ['coh','adm','ph','psd']:
            if m=='psd':im=['Z','1','2','P']
            else:im=['21', 'Z1', 'P1', 'Z2', 'P2', 'ZP']
            tmp={}
            for c in im:
                tmp.update({c:np.array([e[c] for e in day[m]])})
            cpa.update({m:tmp})
        day = AttribDict(cpa)
        data = AttribDict({'sta':sta,'day':day,'f':d.f})
        write_pickle(outfold/f'{stanm}.Noise.pkl',data)
    except:
        print('Issue:' + stanm)
        issues.append(stanm)
if len(issues)>0:print('ISSUES:');issues