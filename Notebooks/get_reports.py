from modules import *

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
def get_reports(comp,catalog,Archive,AVG=True):
    Report=AttribDict();Mirror=AttribDict()
    AnalysisFolder = Archive/'Analysis'/'NetworkCoherences'
    method='atacr';file = AnalysisFolder /f'{method.lower()}'/'Complete'/f'complete_{method}.{comp}_coh.report.pkl'
    Report.ATaCR=load_pickle(file)
    # method='Uncorrected';file = AnalysisFolder /f'{method.lower()}'/'Complete'/f'complete_{method}.{comp}_coh.report.pkl'
    # Report.Uncorrected=load_pickle(file)
    method='hps';file = AnalysisFolder /f'{method.lower()}'/'Complete'/f'complete_{method}.{comp}_coh.report.pkl'
    Report.Noisecut=load_pickle(file)
    for stanm in catalog.StaName:
        net='n'+stanm.split('.')[0]
        sta=stanm.split('.')[1]
        state=f'{stanm} not mirrored in update'
        if not np.isin(net,np.intersect1d(list(Report.ATaCR.keys()),list(Report.Noisecut.keys()))):print(state);continue
        elif not np.isin(sta,np.intersect1d(list(Report.ATaCR[net].keys()),list(Report.Noisecut[net].keys()))):print(state);continue

        c,iatcr,inc=np.intersect1d(
        Report.ATaCR['n'+stanm.split('.')[0]][stanm.split('.')[1]].events,
        Report.Noisecut['n'+stanm.split('.')[0]][stanm.split('.')[1]].events,
        return_indices=True)
        # e=[e.replace('.HDH.SAC','').replace('.H2.SAC','').replace('.HZ.SAC','').replace('.H1.SAC','') for e in Report.Uncorrected['n'+stanm.split('.')[0]][stanm.split('.')[1]].events]
        # Report.Uncorrected['n'+stanm.split('.')[0]][stanm.split('.')[1]].events=e
        # c,ic,uni = np.intersect1d(c,e,return_indices=True)
        # Mirror[stanm]=[[],iatcr,inc]
        n,s=('n'+stanm.split('.')[0]),stanm.split('.')[1]
        # Report.Uncorrected[n][s].coh=Report.Uncorrected[n][s].coh[uni,:]
        # Report.Uncorrected[n][s].events=np.array(Report.Uncorrected[n][s].events)[uni]
        Report.Noisecut[n][s].coh=Report.Noisecut[n][s].coh[inc,:]
        Report.Noisecut[n][s].events=np.array(Report.Noisecut[n][s].events)[inc]
        Report.ATaCR[n][s].coh=Report.ATaCR[n][s].coh[iatcr,:]
        Report.ATaCR[n][s].events=np.array(Report.ATaCR[n][s].events)[iatcr]
    Report.f=Report.ATaCR.f;f=Report.f
    Coherences = AttribDict();Coherences.f = Report.f;Coherences.Depths = catalog.StaDepth
    # Coherences.Uncorrected=AttribDict()
    Coherences.ATaCR=AttribDict();Coherences.Noisecut=AttribDict()
    for stanm in catalog.StaName:
        n,s = stanm.split('.')
        state=f'{stanm} not mirrored in update'
        if not np.isin('n'+n,np.intersect1d(list(Report.ATaCR.keys()),list(Report.Noisecut.keys()))):print(state);continue
        elif not np.isin(s,np.intersect1d(list(Report.ATaCR['n'+n].keys()),list(Report.Noisecut['n'+n].keys()))):print(state);continue
        fn=fnotch(catalog[catalog.StaName==stanm].iloc[0].StaDepth);fni=f>fn
        # Coherences.Uncorrected[stanm]=AttribDict();
        Coherences.ATaCR[stanm]=AttribDict();Coherences.Noisecut[stanm]=AttribDict()
        # Coherences.Uncorrected[stanm].coh=Report.Uncorrected[f'n{n}'][s].coh
        Coherences.ATaCR[stanm].coh=Report.ATaCR[f'n{n}'][s].coh
        if comp=='ZZ':Coherences.ATaCR[stanm].coh[:,fni]=1
        Coherences.Noisecut[stanm].coh= Report.Noisecut[f'n{n}'][s].coh
        # Coherences.Uncorrected[stanm].events=Report.Uncorrected[f'n{n}'][s].events
        Coherences.ATaCR[stanm].events=Report.ATaCR[f'n{n}'][s].events
        Coherences.Noisecut[stanm].events= Report.Noisecut[f'n{n}'][s].events
        id=np.sum([(np.isinf(Coherences.ATaCR[stanm].coh).sum(axis=1)>0),
        (np.isinf(Coherences.Noisecut[stanm].coh).sum(axis=1)>0)],axis=0)==0
        # Coherences.Uncorrected[stanm].coh=Coherences.Uncorrected[stanm].coh[id,:]
        Coherences.ATaCR[stanm].coh=Coherences.ATaCR[stanm].coh[id,:]
        Coherences.Noisecut[stanm].coh=Coherences.Noisecut[stanm].coh[id,:]
        if AVG:
            # Coherences.Uncorrected[stanm].coh=Coherences.Uncorrected[stanm].coh.mean(axis=0)
            Coherences.Noisecut[stanm].coh=Coherences.Noisecut[stanm].coh.mean(axis=0)
            Coherences.ATaCR[stanm].coh=Coherences.ATaCR[stanm].coh.mean(axis=0)
    Report=Coherences
    return Report