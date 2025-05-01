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


def get_reports(comp,catalog,Archive,dirs,AVG=True,methods=['ATaCR','NoiseCut']):
    Report=AttribDict();Mirror=AttribDict()
    AnalysisFolder = Archive/'Analysis'/'NetworkCoherences'
    if np.isin('ATaCR'.lower(),[m.lower() for m in methods]):
        method='atacr';file = AnalysisFolder /f'{method.lower()}'/'Complete'/f'complete_{method}.{comp}_coh.report.pkl'
        Report.ATaCR=load_pickle(file)
    if np.isin('NoiseCut'.lower(),[m.lower() for m in methods]):
        method='hps';file = AnalysisFolder /f'{method.lower()}'/'Complete'/f'complete_{method}.{comp}_coh.report.pkl'
        Report.NoiseCut=load_pickle(file)
    for stanm in catalog.StaName:
        n,s=('n'+stanm.split('.')[0]),stanm.split('.')[1]
        net='n'+stanm.split('.')[0]
        sta=stanm.split('.')[1]
        state=f'{stanm} not mirrored in update'
        if len(methods)==2:
            if not np.isin(net,np.intersect1d(list(Report.ATaCR.keys()),list(Report.NoiseCut.keys()))):print(state);continue
            elif not np.isin(sta,np.intersect1d(list(Report.ATaCR[net].keys()),list(Report.NoiseCut[net].keys()))):print(state);continue
            c,iatcr,inc=np.intersect1d(
            Report.ATaCR['n'+stanm.split('.')[0]][stanm.split('.')[1]].events,
            Report.NoiseCut['n'+stanm.split('.')[0]][stanm.split('.')[1]].events,
            return_indices=True)
            Report.NoiseCut[n][s].coh=Report.NoiseCut[n][s].coh[inc,:]
            Report.NoiseCut[n][s].events=np.array(Report.NoiseCut[n][s].events)[inc]
            Report.ATaCR[n][s].coh=Report.ATaCR[n][s].coh[iatcr,:]
            Report.ATaCR[n][s].events=np.array(Report.ATaCR[n][s].events)[iatcr]
        else:
            inc = np.intersect1d(Report[methods[0]][n][s].events,[e.Name for e in lt.cat.mirror(stanm,catalog,dirs)],return_indices=True)[1]
            m=methods[0]
            Report[m][n][s].coh=Report[m][n][s].coh[inc,:]
            Report[m][n][s].events=np.array(Report[m][n][s].events)[inc]
    Report.f=Report[methods[0]].f;f=Report.f
    for m in methods:Report.f=Report[m].f;f=Report.f
    temp=AttribDict();temp.f=Report.f
    for m in methods:
        temp[m]=AttribDict()
        for stanm in catalog.StaName:
            n,s = stanm.split('.')
            temp[m][stanm]=Report[m][f'n{n}'][s]
    Report = temp
    return Report

def get_reports_adm(comp,catalog,Archive,dirs,AVG=True,methods=['ATaCR','NoiseCut']):
    Report=AttribDict();Mirror=AttribDict()
    AnalysisFolder = Archive/'Analysis'/'NetworkAdmittance'
    if np.isin('ATaCR'.lower(),[m.lower() for m in methods]):
        method='atacr';file = AnalysisFolder /f'{method.lower()}'/'Complete'/f'complete_{method}.{comp}_adm.report.pkl'
        Report.ATaCR=load_pickle(file)
    if np.isin('NoiseCut'.lower(),[m.lower() for m in methods]):
        method='hps';file = AnalysisFolder /f'{method.lower()}'/'Complete'/f'complete_{method}.{comp}_adm.report.pkl'
        Report.NoiseCut=load_pickle(file)
    for stanm in catalog.StaName:
        n,s=('n'+stanm.split('.')[0]),stanm.split('.')[1]
        net='n'+stanm.split('.')[0]
        sta=stanm.split('.')[1]
        state=f'{stanm} not mirrored in update'
        if len(methods)==2:
            if not np.isin(net,np.intersect1d(list(Report.ATaCR.keys()),list(Report.NoiseCut.keys()))):print(state);continue
            elif not np.isin(sta,np.intersect1d(list(Report.ATaCR[net].keys()),list(Report.NoiseCut[net].keys()))):print(state);continue
            c,iatcr,inc=np.intersect1d(
            Report.ATaCR['n'+stanm.split('.')[0]][stanm.split('.')[1]].events,
            Report.NoiseCut['n'+stanm.split('.')[0]][stanm.split('.')[1]].events,
            return_indices=True)
            Report.NoiseCut[n][s].adm=Report.NoiseCut[n][s].adm[inc,:]
            Report.NoiseCut[n][s].events=np.array(Report.NoiseCut[n][s].events)[inc]
            Report.ATaCR[n][s].adm=Report.ATaCR[n][s].adm[iatcr,:]
            Report.ATaCR[n][s].events=np.array(Report.ATaCR[n][s].events)[iatcr]
        else:
            inc = np.intersect1d(Report[methods[0]][n][s].events,[e.Name for e in lt.cat.mirror(stanm,catalog,dirs)],return_indices=True)[1]
            m=methods[0]
            Report[m][n][s].adm=Report[m][n][s].adm[inc,:]
            Report[m][n][s].events=np.array(Report[m][n][s].events)[inc]
    Report.f=Report[methods[0]].f;f=Report.f
    for m in methods:Report.f=Report[m].f;f=Report.f
    temp=AttribDict();temp.f=Report.f
    for m in methods:
        temp[m]=AttribDict()
        for stanm in catalog.StaName:
            n,s = stanm.split('.')
            temp[m][stanm]=Report[m][f'n{n}'][s]
    Report = temp
    return Report