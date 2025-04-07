from imports import *

w = load_pickle(dirs.Analysis/'Latest_DataAvailabilityWarnings.pkl')

aef = dirs.Events/'corrected'
hef = dirs.Events_HPS/'corrected'

not_in_catalog = [wi for wi in w if wi.split(':')[-1]=='Event not found in catalog']
not_found = [wi for wi in w if wi.split(':')[-1]=='Catalog event not found']

# for i,wi in enumerate(not_in_catalog):
#     m,s,e,reason = wi.split(':')
#     if reason=='Event not found in catalog':
#         for c in ['raw','rmresp','corrected','plots']:
#             fold = (dirs.Events if m=='atacr' else dirs.Events_HPS)
#             fold=fold/c/s
#             files=list((fold).glob(f'*{e}*'))
#             if fold.exists():
#                 deprfold=(fold/'_depreciated')
#                 deprfold.mkdir(parents=True,exist_ok=True)
#                 [shutil.move(fold/f.name,deprfold/f.name) for f in files]
#     print(f'{i+1} / {len(w)} : {s} : {np.round(100*((i+1)/len(w)),2)}%')


# for i,wi in enumerate(not_found):
#     m,s,e,reason = wi.split(':')
#     if reason=='Catalog event not found':
#         for c in ['raw','rmresp','corrected','plots']:
#             fold = (dirs.Events if m=='atacr' else dirs.Events_HPS)
#             fold=fold/c/s
#             files=list((fold).glob(f'*{e}*'))
#             if fold.exists():
#                 deprfold=(fold/'_depreciated')
#                 deprfold.mkdir(parents=True,exist_ok=True)
#                 [shutil.move(fold/f.name,deprfold/f.name) for f in files]
#     print(f'{i+1} / {len(w)} : {s} : {np.round(100*((i+1)/len(w)),2)}%')