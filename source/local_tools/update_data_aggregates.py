import sys;from pathlib import Path
# sys.path.append(str(Path(__file__).parent.parent))
import os,sys
# from imports import * #Standard imports for doing anything in this project. Approx. ~19 seconds.
from modules import *
from scipy.stats import iqr
from local_tools.quick_class import *
from local_tools.math import spectra
from obspy.geodetics import kilometers2degrees
# cat = catalog.copy()
octavg=lt.math.octave_average
import statsmodels.api as sm
import local_tools.dataspace as ds

cat=ds.dataspace(bulk=False)

dirs = io.dir_libraries()

bulkfile='BulkLoad.SR.Coherences_50125.pkl'
df=load_pickle(dirs.Analysis/bulkfile)
df=df.iloc[np.where(np.isin(np.array([df.Name + '_' + df.StaName]),np.array([cat.sr.Name + '_' + cat.sr.StaName])))[1]]
cat.sr['Coherence']=[[] for _ in np.arange(len(cat.sr))]
coh=[]
for sr in cat.sr.iloc:
    ii=df.loc[sr.Name][df.loc[sr.Name].StaName==sr.StaName]
    assert len(ii)==1
    ii=ii.iloc[0]
    assert (ii.Name==sr.Name)&(ii.StaName==sr.StaName)
    coh.append(ii.Coherence.copy())
cat.sr['Coherence']=coh
cat.sr[df.columns].to_pickle(dirs.Analysis/bulkfile)

assert sum(np.array([d.Coherence.event==d.Name for d in cat.sr.iloc])==False)==0
assert sum(np.array([d.Coherence.StaName==d.StaName for d in cat.sr.iloc])==False)==0

# c_coh=[]
# for si,s in enumerate(cat.sr.iloc):
#     print(f'{si+1}/{len(cat.sr)}')
#     c=s.Data.Coherence()
#     coh=AttribDict({})
#     coh['TF']=c.ATaCR.zp_21.coh
#     coh['HPS_Z']=c.HPS.zz.coh
#     coh['HPS_H']=np.mean([c.HPS['11'].coh,c.HPS['22'].coh],axis=0)
#     coh['event']=s.Name
#     coh['StaName']=s.StaName
#     c_coh.append(coh)
# cat.sr['Coherence']=c_coh
# cat.sr[df.columns].to_pickle(dirs.Analysis/'BulkLoad.SR.Coherences_50125.pkl')
# if cat.bulk_coherence:


filtertype='acausul';note='V04';fold=dirs.SNR/'SNR.Models';file =f'SNR_{filtertype}.filter_{note}.pkl'
SNR=load_pickle(fold/file)
SNR=SNR.iloc[np.where(np.isin(np.array([SNR.Name + '_' + SNR.StaName]),np.array([cat.sr.Name + '_' + cat.sr.StaName])))[1]]
assert len(SNR)==len(cat.sr)
SNR=SNR.iloc[[np.where(np.array([SNR.Name + '_' + SNR.StaName])==(sr.Name+'_'+sr.StaName))[-1][0] for sr in cat.sr.iloc]]
assert np.sum(np.array([SNR.Name+'_'+SNR.StaName])==np.array([cat.sr.Name+'_'+cat.sr.StaName]))==len(cat.sr)
SNR.to_pickle(fold/file)

k=1