from imports import *
_=os.system('cls' if os.name == 'nt' else 'clear')
import time as pause
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
ATaCR_Parent = dirs.ATaCR
cat = catalog.copy()
# windows=[['2010/06/01','2010/11/30'],
# ['2012/08/01','2012/11/30'],
# ['2013/01/01','2013/07/30'],
# ['2014/06/01','2014/08/30'],
# ['2015/01/01','2015/05/30']]
# minmag,maxmag=6.5,6.75


windows=[['2009/12/01','2019/05/01']] #<---Entire dataset window
minmag,maxmag=6.26,6.4

nsta=0
windows=np.array([[UTCDateTime.strptime(t[0],'%Y/%m/%d'),UTCDateTime.strptime(t[1],'%Y/%m/%d')+(24*3600)] for t in windows])
minsta,maxsta=1,5
maxev=300;nev=42;norm_deviation=1.3417556034384266
client = Client("IRIS")
message = lambda:print(f'Total events collected {nev}/{maxev} \n| Current norm-dist. score {norm_deviation} \n| Search Window {wi+1}/{len(windows)} \n| Candidate Event {ei+1}/{len(events)} \n| Current station density {nsta}/{stabar} \n| Candidate Station {si+1}/{len(ecat)}')
for wi,(starttime,endtime) in enumerate(windows):
    runningcatfile=dirs.Catalogs/'TempRunningEventCatalog.pkl'
    runningcatalog=pd.read_pickle(runningcatfile)
    evcat=[e for e in runningcatalog.Events]

    binary=1
    # if nev>=40:minmag=6.35
    # events = client.get_events(starttime=starttime, endtime=endtime, minmagnitude=minmag,maxmagnitude=maxmag,orderby='magnitude')
    events = client.get_events(starttime=starttime, endtime=endtime, minmagnitude=minmag,maxmagnitude=maxmag,orderby='time')
    # Enforce a minimum 2-hour distance between subsequent events
    keep=[]
    for i,ev in enumerate(events):
        if i==0:keep.append(i);t0=ev.origins[0].time-7200
        elif ev.origins[0].time<t0:keep.append(i);t0=ev.origins[0].time-7200
        ev.Name=ev.origins[0].time.strftime('%Y.%j.%H.%M')
    events = Catalog([events[i] for i in keep])
    # Randomly re-order the events
    events = Catalog([events[i] for i in np.random.choice([a for a in range(len(events))],len(events),replace=False)])
    for ei,ev in enumerate(events):

        if np.isin(ev.Name,[e.Name for e in evcat]):continue

        stations_collected = []
        binary=-1
        if nev>=maxev:continue
        stabar=maxsta
        ecat=cat[(cat.Start<ev.origins[0].time) & (cat.End>(ev.origins[0].time+7200))].copy()
        if len(ecat)<minsta:si=0;message();print(f'Insufficient number of stations. Skipping');continue
        nsta=0
        stabar=np.random.randint(minsta,maxsta+1,1)[0]
        if binary<0:ecat=ecat.iloc[::-1]
        for si,Station in enumerate(ecat.iloc):
            print(f'{Station.StaName}')
            icatalog = Station.to_frame().T
            icatalog.Events.iloc[0] = Catalog([ev])
            # _=os.system('cls' if os.name == 'nt' else 'clear')
            if nsta>=stabar:continue
            print('=='*15);print('=='*15);print(' ');message();print(' ');print('=='*15);print('=='*15)
            ObsQA.TOOLS.io.DownloadEvents(icatalog,
            ATaCR_Parent=ATaCR_Parent,
            Minmag=minmag,Maxmag=maxmag,
            logoutput_subfolder='',
            staquery_output='./sta_query.pkl',
            chan='H',
            log_prefix=Station.StaName,
            channels='Z,P,12')
            pause.sleep(0.1)
            Result=pd.read_pickle(dirs.ATaCR/'RunTest.pkl').Result
            if not Result=='Fail':
                nsta+=1
                stations_collected.append(Station.StaName)
            print(f'Result:{Result}')
        if nsta>=0:
            # nev+=1
            ev.Stations=stations_collected
            runningcatfile=dirs.Catalogs/'TempRunningEventCatalog.pkl'
            runningcatalog=pd.read_pickle(runningcatfile)
            evcat=[e for e in runningcatalog.Events]
            if not np.isin(ev.Name,[e.Name for e in evcat]):
                evcat.append(ev)
                nev=len(evcat)
                runningcatalog.Events=Catalog(evcat)
                runningcatalog.to_pickle(runningcatfile)
                norm_deviation=np.mean([np.std([len(e.Stations) for e in evcat]),np.mean([len(e.Stations) for e in evcat])/np.mean([a for a in np.arange(minsta,maxsta+1)])])
            fig,axes=plt.subplots(nrows=2,ncols=1,figsize=(19,8))
            past=[j.magnitudes[0].mag for j in [unravel([e for e in cat.Events if e is not None])[k] for k in np.unique([e.Name for e in unravel([e for e in cat.Events if e is not None])],return_index=True)[1]]]
            update=past.copy()
            update.extend([e.magnitudes[0].mag for e in evcat])
            axes[0].hist(update,bins=[6,6.26,6.5,6.75,7.0,7.25,7.5,7.75,8.0,8.25])
            axes[0].hist(past,bins=[6,6.26,6.5,6.75,7.0,7.25,7.5,7.75,8.0,8.25])
            axes[0].set_xlabel('Mag')
            axes[0].set_ylabel('Events')
            past=[len(j.Stations) for j in [unravel([e for e in cat.Events if e is not None])[k] for k in np.unique([e.Name for e in unravel([e for e in cat.Events if e is not None])],return_index=True)[1]]]
            update=past.copy()
            update.extend([len(e.Stations) for e in evcat])

            axes[1].hist(update,bins=[i+1 for i in range(45)])
            axes[1].hist(past,bins=[i+1 for i in range(45)])
            axes[1].set_xticks([i+1 for i in range(45)])
            axes[1].set_xlabel('Station density')
            axes[1].set_ylabel('Events')
            save_tight(dirs.Plots/'hist.running.updates.png',fig)
            plt.close()