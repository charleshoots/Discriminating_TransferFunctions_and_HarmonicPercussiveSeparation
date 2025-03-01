from imports import *

eventfolder = dirs.Events_HPS
# eventfolder = dirs.Events
[(eventfolder/'raw'/stanm).mkdir(parents=True,exist_ok=True) for si,stanm in enumerate(catalog.StaName)]
[(eventfolder/'rmresp'/stanm).mkdir(parents=True,exist_ok=True) for si,stanm in enumerate(catalog.StaName)]
[(eventfolder/'corrected'/stanm).mkdir(parents=True,exist_ok=True) for si,stanm in enumerate(catalog.StaName)]
[(eventfolder/'plots'/stanm).mkdir(parents=True,exist_ok=True) for si,stanm in enumerate(catalog.StaName)]
for si,stanm in enumerate(catalog.StaName):
    print(f'{si+1}/{len(catalog)} | {stanm}')
    # --------------------
    # source_fold = eventfolder/stanm/'CORRECTED'
    # dest_fold = eventfolder/'corrected'/stanm
    # files = [f.name for f in list(source_fold.glob('*')) if (not f.is_dir()) and not (f.name=='.DS_Store')]
    # [shutil.move(source_fold/f,dest_fold/f) for f in files]
    # --------------------
    # source_fold = eventfolder/stanm/'raw'
    # dest_fold = eventfolder/'raw'/stanm
    # files = [f.name for f in list(source_fold.glob('*')) if (not f.is_dir()) and not (f.name=='.DS_Store')]
    # [shutil.move(source_fold/f,dest_fold/f) for f in files]
    # --------------------
    source_fold = eventfolder/stanm/'rmresp'
    dest_fold = eventfolder/'rmresp'/stanm
    files = [f.name for f in list((source_fold).glob('*.SAC*'))]
    [shutil.move(source_fold/f,dest_fold/f) for f in files]
    # [load_sac(source_fold/f,rmresp=True)[0].write(str(dest_fold/f),format='SAC') for f in files if not (dest_fold/f).exists()]
    # --------------------
    # source_fold = eventfolder/stanm/'CORRECTED'/'EventRecords'
    # dest_fold = eventfolder/'plots'/stanm
    # dest_fold.mkdir(parents=True,exist_ok=True)
    # files = [f.name for f in list((source_fold).glob('*.png'))]
    # [shutil.move(source_fold/f,dest_fold/f) for f in files]
    clear_output(wait=False);os.system('cls' if os.name == 'nt' else 'clear')