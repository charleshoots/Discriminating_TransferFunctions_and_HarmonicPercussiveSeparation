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

import os,sys
import matplotlib.cm as cm



# cat value
cat = catalog.r.copy()
# plotfold value
plotfold = dirs.P01.S01
# fold value
fold = dirs.Analysis/'Metrics'/'noise'
# get stanm noise value
get_stanm_noise = lambda fold,stanm:load_pickle(fold/f'{stanm}.Noise.pkl')
# load sta atacrnoise value
load_sta_atacrnoise = lambda stanm: load_pickle(list((dirs.SpectraAvg/stanm).glob('*.pkl'))[0])
# cat = cat[cat.Network=='2D']
get_stanoise = lambda stanm: load_pickle(list((dirs.TransferFunctions/stanm).glob('*-*.pkl'))[0])

# variable
_=[(dirs.P01.S01/i).mkdir(exist_ok=True) for i in ['COH_PSD','PH_ADM','TF']]

for si,sta in enumerate(cat.iloc):
    # noise value
    noise = get_stanm_noise(dirs.Analysis/'Metrics'/'noise',sta.StaName)
    fo = plotfold/'COH_PSD'/f'{sta.StaName}.Noise_COH_PSD.png'
    save_tight(fo,coh_psd(noise,sta),dpi=300)
    plt.close();print(f'{si+1}/{len(cat)} | {sta.StaName} | Saved to: {f'{fo.parent.parent.name}/{fo.parent.name}/{fo.name}'}')
    fo = plotfold/'PH_ADM'/f'{sta.StaName}.Noise_PH_ADM.png'
    save_tight(fo,ph_adm(noise,sta),dpi=300)
    plt.close();print(f'{si+1}/{len(cat)} | {sta.StaName} | Saved to: {f'{fo.parent.parent.name}/{fo.parent.name}/{fo.name}'}')
    fo = plotfold/'TF'/f'{sta.StaName}.TFs.png'
    save_tight(fo,plot_TFs(sta.StaName),dpi=300)
    plt.close();print(f'{si+1}/{len(cat)} | {sta.StaName} | Saved to: {f'{fo.parent.parent.name}/{fo.parent.name}/{fo.name}'}')

from PIL import Image
fold = plotfold
for tag in ['COH_PSD','PH_ADM']:
    print(f'Creating PDFs for {tag}')
    for n in cat.Network.unique():
        icat = cat[cat.Network==n]
        icat.sort_values(by='StaDepth',inplace=True)
        # files = list(fold.glob(f'{n}.*{tag}.png'))
        files = [fold/tag/f'{stanm}.Noise_{tag}.png' for stanm in icat.StaName]
        if tag=='COH_PSD':
            files.extend([fold/'TF'/f'{stanm}.TFs.png' for stanm in icat.StaName])
        # files = [Path(f) for f in sorted([str(g) for g in files])]
        files = [files[i] for i in np.argsort([f.name for f in files])]
        # Ensure all images are in RGB mode
        images = [Image.open(png).convert('RGB') for png in files]
        # Save the images as a single PDF
        output_pdf = fold/'_PDFs'/f'{tag}_Noise_{n}.pdf'
        images[0].save(output_pdf, save_all=True, append_images=images[1:])


lt.cat.banner('S01 COMPLETE!')