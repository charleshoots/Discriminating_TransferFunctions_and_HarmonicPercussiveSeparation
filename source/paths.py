from pathlib import Path;import shutil,sys,os
from obspy.core import AttribDict
def dir_libraries(project_path=None,mkdir=False):
        if project_path==None:project_path=Path(os.getcwd())
        project_path=Path(project_path)
        ATaCR_Paths = dict()
        ATaCR_Paths['Py_DataParentFolder'] = project_path /'_DataArchive'/'ATaCR_Data'/'ATaCR_Python'
        ATaCR_Paths['Py_RawDayData'] = ATaCR_Paths['Py_DataParentFolder'] / 'Data'
        ATaCR_Paths['Py_StaSpecAvg'] = ATaCR_Paths['Py_DataParentFolder'] / 'AVG_STA'
        ATaCR_Paths['Py_CorrectedTraces'] = ATaCR_Paths['Py_DataParentFolder'] / 'EVENTS'
        ATaCR_Paths['Py_b1b2_StaSpectra'] = ATaCR_Paths['Py_DataParentFolder'] / 'SPECTRA'
        ATaCR_Paths['Py_TransferFunctions'] = ATaCR_Paths['Py_DataParentFolder'] / 'TF_STA'
        ATaCR_Paths['Py_Logs'] = ATaCR_Paths['Py_DataParentFolder'] / 'Logs'
        d = AttribDict()
        d.Root = project_path
        d.Data = d.Root /'_DataArchive'
        d.ATaCR = Path(ATaCR_Paths['Py_DataParentFolder'])
        d.Catalogs = d.Data / 'Catalogs'
        d.Events = Path(ATaCR_Paths['Py_CorrectedTraces'])
        d.Spectra = Path(ATaCR_Paths['Py_b1b2_StaSpectra'])
        d.SpectraAvg = Path(ATaCR_Paths['Py_StaSpecAvg'])
        d.TransferFunctions = Path(ATaCR_Paths['Py_TransferFunctions'])
        d.Noise = Path(ATaCR_Paths['Py_RawDayData'])
        d.NoiseTrace = Path(ATaCR_Paths['Py_DataParentFolder']) / 'Noise'
        d.Logs = Path(ATaCR_Paths['Py_Logs'])

        plotfolder = Path(str(d.Root / '_FigureArchive'))
        d.Scripts=d.Root/'Notebooks'/'01_Analysis'
        d.Plots = plotfolder
        d.Events_HPS = d.Data / 'HPS_Data' / 'Data'
        d.Analysis=d.Data/'Analysis'
        d.SNR=d.Scripts/'S00_Collect.SNR'
        d.P01 = AttribDict({})
        d.P01.Parent = d.Data/'P01_Analysis_Figures'
        d.P01.S01 = d.P01.Parent/'S01_NoisePlots'
        d.P01.S02 = d.P01.Parent/'S02_StationPages'
        d.P01.S03 = d.P01.Parent/'S03_HPSPlots'
        d.P01.S04 = d.P01.Parent/'S04_CoherenceConour'
        d.P01.S05 = d.P01.Parent/'S05_CohvCohScatters'
        d.P01.S06 = d.P01.Parent/'S06_StemPlots'
        d.P01.S07 = d.P01.Parent/'S07_EventRecords_Traces'
        d.P01.S08 = d.P01.Parent/'S08_EventRecords_Metrics'
        d.P01.S09 = d.P01.Parent/'S09_Station_Noise_QC'
        d.P01.S10 = d.P01.Parent/'S10_Coherence_Histograms'
        d.P01.S11 = d.P01.Parent/'S11_SpectraPlot.Coherence_ColoredByDepth.SortedByMeta'
        d.P01.S12 = d.P01.Parent/'S12_Scatters_Noise.and.Z.vs.Coherence'
        d.P01.S99 = d.P01.Parent/'S99_Summaries'
        d.Papers=d.Root/'_FigureArchive'/'_Papers'
        d.Ch1=d.Root/'Notebooks'/'_03_Run_Paper.Figures'/'ImageOutputs'

        if mkdir:
                # Generates complete directory structure. Does nothing if it already exists
                _=[d[fo].mkdir(exist_ok=True,parents=True) for fo in list(d.keys()) if not fo=='P01']
                _=[d.P01[fo].mkdir(exist_ok=True,parents=True) for fo in list(d.P01.keys())]

        return d