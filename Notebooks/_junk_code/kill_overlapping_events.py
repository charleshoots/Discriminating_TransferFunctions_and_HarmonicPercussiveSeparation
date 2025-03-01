from imports import *
k = 1
# The hps dead trace cull ultimately killed 7 events from the entire catalog
Events = unravel([s.Events for s in cat.iloc])

unique_events = Catalog([Events[i] for i in np.unique([e.Name for e in Events],return_index=True)[1]])
event_names = np.sort(np.unique([e.Name for e in unique_events]))
events_utc = np.array([UTCDateTime.strptime(e,dateformat) for e in event_names])


throw = []
for ei,e in enumerate(events_utc):
    suspects = abs(np.array((events_utc - (e+7200))))<7200
    inds = np.where(suspects)[0]
    suspects = [i for i in inds if i>ei]

    if len(suspects)>0:
        i_sus = 0
        sus = suspects[i_sus]
        sus_mag = unique_events[sus].magnitudes[0].mag
        ei_mag = unique_events[ei].magnitudes[0].mag
        if ei_mag>=sus_mag:
            throw.append(sus)
        else:
            throw.append(ei)
kill_events = np.unique([unique_events[i].Name for i in throw])
# for si,sta in enumerate(cat.StaName):
#     ii = np.where(cat.StaName==sta)[0][0]
#     Events = cat.iat[ii,np.where(cat.columns=='Events')[0][0]]
#     Events = Catalog([e for e in Events if e.Name not in kill_events])
#     if len(Events)<17:
#         print('hey')
#     cat.iat[ii,np.where(cat.columns=='Events')[0][0]] = Events