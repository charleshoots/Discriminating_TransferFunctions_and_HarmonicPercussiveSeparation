from modules import *
from modules import *
def unravel(lst):return list(itertools.chain.from_iterable(lst))
def unravel_report(file,AVG=False):
    r = load_pickle(file)
    nets=[n for n in list(r.keys()) if n[0]=='n']
    d=AttribDict()
    for n in nets:
        stas=list(r[n].keys())
        for s in stas:
            stanm=f'{n[1:]}.{s}'
            d[stanm]=AttribDict()
            events=np.array([e.replace('.HDH.SAC','').replace('.H2.SAC','').replace('.HZ.SAC','').replace('.H1.SAC','') 
            for e in r[n][s].events])
            d[stanm].events=events
            d[stanm].coh=r[n][s].coh
            keep=np.isnan(d[stanm].coh).sum(axis=1)==0
            d[stanm].coh=d[stanm].coh[keep,:]
            d[stanm].events=d[stanm].events[keep]
            d.f=r.f
            d[stanm].f = d.f
            if AVG:d[stanm].coh=d[stanm].coh.mean(axis=0)
    return d
def load_pickle(file):
    import pickle
    with open(file, 'rb') as handle:
        b = pickle.load(handle)
    return b


import numpy as np
from functools import reduce

def reshapedict(c):
    n = np.array(list(c.keys()));n=n[n!='f']
    sn = unravel([[f'{ni[1:]}.{s}' for s in c[ni].keys()] for ni in n])
    r = AttribDict({f'{s}':c[f'n{s.split('.')[0]}'][s.split('.')[1]] for s in sn})
    r.f = c.f
    return r
def get_reports(sta=None,
    Archive=None,
    methods=['ATaCR','NoiseCut'],evna=None):
    from local_tools.io import dir_libraries
    if Archive==None:dirs=dir_libraries();Archive = dirs.Data
    methods = np.array(methods)
    Report=AttribDict()
    AnalysisFolder = Archive/'Analysis'/'Coherence'
    lp = lambda f:reshapedict(load_pickle(f))
    atacrfile = lambda:AnalysisFolder/f'{m.lower()}'/f'{m.lower()}.{tf}.zz.coherence.pkl'
    hpsfile = lambda :AnalysisFolder/f'{m.lower()}'/f'{m.lower()}.{c}.coherence.pkl'
    unionevents = lambda a:reduce(np.intersect1d,a)
    methods = [m.replace('NoiseCut','HPS') for m in methods]
    for m in methods:
        if m.lower()=='atacr':
            tfs = ['zp_21']
            d = AttribDict()
            for tf in tfs:d[tf]= lp(atacrfile())
            Report[m] = AttribDict()
            for tf in tfs:
                Report[m][tf] = AttribDict()
                stas = list(d[tf].keys()) if sta==None else [sta]
                for s in stas:
                    Report[m][tf][s] = d[tf][s]
                    if s=='f':continue
                    evs=unionevents([Report[m][tf][s].events for i in tfs])
                    ind = np.isin(Report[m][tf][s].events,evs)
                    Report[m][tf][s].coh=Report[m][tf][s].coh[ind,:]
                    Report[m][tf][s].events=np.array(Report[m][tf][s].events)[ind]
                    Report.f = d[tf].f
        elif m.lower()=='hps':
            Report[m]=AttribDict()
            d = AttribDict()
            comps = ['zz','11','22']
            for c in comps:d[c]= lp(hpsfile())
            for c in comps:Report[m][c]=AttribDict()
            Report.f = d[c].f
            stas = list(d[c].keys()) if sta==None else [sta]
            for c in comps:
                for s in stas:Report[m][c][s] = d[c][s]
            for s in stas:
                if s=='f':continue
                evs=unionevents([Report[m][c][s].events for c in comps])
                for c in comps:
                    ind = np.isin(Report[m][c][s].events,evs)
                    Report[m][c][s].coh=Report[m][c][s].coh[ind,:]
                    Report[m][c][s].events=np.array(Report[m][c][s].events)[ind]
    for m in methods:
        c = list(Report[m].keys())[0]
        for s in stas:
            if s=='f':continue
            evs=unionevents([Report[m][list(Report[m].keys())[0]][s].events for m in methods])
            ind = np.isin(Report[m][c][s].events,evs)
            for i in Report[m].keys():
                Report[m][i][s].coh = Report[m][i][s].coh[ind,:]
                Report[m][i][s].events = np.array(Report[m][i][s].events)[ind]
    if len(stas)==1:
        for m in methods:
            for k in Report[m].keys():
                c = Report[m][k][stas[0]].copy()
                Report[m][k] = AttribDict()
                Report[m][k].coh = c.coh
                Report[m][k].events = c.events
    
    if not (evna==None):
        Report  = lt.cat.evcoh(Report,evna)
    return Report