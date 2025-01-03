from imports import *


# def AuditEventFolder(cat,eventsfolder,parseby='*Z.SAC',Minmag=6.0,Maxmag=8.0):
#     client = Client()
#     for ista,Station in enumerate(cat.iloc):
#         stanm=Station.StaName
#         stafolder = Path(eventsfolder) / stanm
#         files = list(stafolder.glob(parseby))
#         evna = [f.name for f in files]
#         evna = [f.split(f'{stanm}.')[-1] for f in evna]
#         evna = ['.'.join(f.split('.')[:4]) for f in evna]
#         assert len(evna)==len(np.unique(evna))
#         print(f'{ista+1}/{len(cat)} | {stanm} | {len(evna)} events found. Collecting metadata from IRIS..')
#         events=[]
#         for ev in evna:
#             timedelta = 60
#             start = UTCDateTime.strptime(str(ev),'%Y.%j.%H.%M')
#             end = start + timedelta
#             event = client.get_events(starttime=start, endtime=end,minmagnitude=Minmag, maxmagnitude=Maxmag,orderby='magnitude')[0]
#             event.Name = ev
#             events.append(event)
#         events = Catalog(events)
#         Station.Events = None
#         Station.Events = events
#     all_evnames=[[e.Name for e in Station.Events] for Station in cat.iloc]
#     for ista,Station in enumerate(cat.iloc):
#         evna = [e.Name for e in Station.Events]
#         for Event in Station.Events:
#             Event.Stations = list(np.array(cat.StaName[[i for i,sev in enumerate(all_evnames) if np.isin(Event.Name,sev)]]))
#     return cat

# eventsfolder=dirs.Events/'raw'
# cat = catalog.copy().iloc[:4]
# cat = AuditEventFolder(cat,eventsfolder)
# cat


meta_hist_plot(pd.read_pickle(dirs.Catalogs / 'sta_catalog_111524c.pkl'))