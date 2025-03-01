import os
from PIL import Image
from imports import *

# for si,sta in enumerate(cat.StaName):
files = list(dirs.P01.S07.glob('*.png'))

files=np.unique([f.name.replace('_traces__NoiseCut.png','').replace('_traces__ATaCR.png','') for f in files])

for fi,f in enumerate(files):
    print(f'{fi+1}/{len(files)} | {np.round(100*(fi+1)/len(files),3)}%')
    pdf_filename=str(dirs.P01.S07/'_PDFs'/f'{f}.pdf')
    if Path(pdf_filename).exists():continue
    print(f'{fi+1}/{len(files)}')
    png_files = list(dirs.P01.S07.glob(f'{f}*.png'))
    png_files = [str(i) for i in png_files]
    try:
        # pfold=dirs.P01.S03
        # fold=pfold/sta
        # png_files=list(fold.glob('*.png'))
        # pdf_filename=pfold/'_PDFs'/f'{sta}.pdf'
        # Load images
        images = [Image.open(img).convert('RGB') for img in png_files]
        # Save as PDF
        if images:images[0].save(pdf_filename, save_all=True, append_images=images[1:], format='PDF')
    except:
        continue

# could be decompression bomb DOS attack