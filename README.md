# Comparing noise reduction methods using transfer function or harmonic percussive separation. 

<center><img src="_docs/OverviewThumbnail.png" width=1000/></center>

## Installation: 

```
git clone https://github.com/charleshoots/Discriminating_TransferFunctions_and_HarmonicPercussiveSeparation.git
```

```
cd Discriminating_TransferFunctions_and_HarmonicPercussiveSeparation
conda env create -f environment.yml
```

```
conda activate Seismic_TF_HPS_Comparison
```

## Download models (optional):
#### - After installation, our models containing coherence and SNR can be downloaded and unzipped into the repositories directory with the following script. 

```
python source/download_models.py
```

## Overview - 

All analysis and data navigation is done with four classes:



| DataSpace |
|----------|

#### A highly detailed catalog containing all receiver and source receiver level data and metadata. Anything in this analysis (traces, entire record sections, coherence, snr, noise spectra, etc.) can be quickly accessed using just this catalog. An indexing method called .loc is built into the Pandas API when calling this class for very quick navigation without requiring knowledge of specific column names.



```
# Example
from source.imports import *
cat = catalog.r.copy() #Catalog of receiver level data/metadata
cat = catalog.sr.copy() #Catalog of source-receivers level data/metadata
trm = cat.loc['TRM'] #All source-receivers using a TRM instrument design
record = cat.loc['2015.115.06.11'] #All source receivers for a specific event name
sr = record.loc['7D.FS42D'] #A specific source-receiver from that event
noise = sr.iloc[0].Data.Noise.Averaged() #Station averaged noise spectra for each component at the associated receiver
hps_record = Stream([r.Traces().select(location='*NoiseCut*')[0] for r in record.iloc]) #All HPS corrected traces for this event
hps_record.plot() # Record section plot
```

| Signal |
|----------|

#### A simple class for calculating basic spectral measurements used in comparing signals such as coherence, phase, admittance, PSD, and CSD.

```
# Example
from source.imports import *
cat = catalog.sr.copy() #Catalog of source-receivers
sr = cat.iloc[0] #Select a source-receiver
st = sr.Traces() #Load the traces
original,corrected=st.select(location='*Original*')[0],st.select(location='*NoiseCut*')[0] #Specify traces to compare
sn = Signal(original,corrected) #Instatiate a Signal class
f,coh = sn.coherence() #Magnitude-squared Coherence between original and corrected
f,phase = sn.phase() #Spectral phase between original and corrected
f,adm = sn.admittance() #Spectral admittance between original and corrected
```

| AggregateMeasurements |
|----------|

#### A fast method of aggregating coherence and snr averages within a band. Can specify ingravity limit sensitivity.

```
# Example
from source.imports import *
cat = catalog.sr.copy() #Catalog of source-receivers
cohsnr=unpack_metrics(cat) #Outputs a AggregateMeasurements class containing all coherence, SNR, signal, and noise measurements.

f = 1/cohsnr.coh.bands #coherence frequency vector
igsensitive = True #Specify infragravity sensitivity. Arguments can be True, False, or None.

#coherence averages for each source-receiver pair within the 1-10s band subset within sensitivity to the infragravity wave.
tfzcoh = cohsnr.coh.TF_Z.Average((1,10),fn='IG' if igsensitive else None)
hpszcoh = cohsnr.coh.HPS_Z.Average((1,10),fn='IG' if igsensitive else None)
hps1coh = cohsnr.coh.HPS_1.Average((1,10),fn='IG' if igsensitive else None)
hps2coh = cohsnr.coh.HPS_2.Average((1,10),fn='IG' if igsensitive else None)

#all source-receiver Rayleigh wave (Rg) measured SNR (snr) ratios (R() argument) averaged over the 30-100s period band.
tfzcsnr = cohsnr.snr.TF_Z.R().Rg.Average((30,100),fn='IG' if igsensitive else None)
hpszsnr = cohsnr.snr.HPS_Z.R().Rg.Average((30,100),fn='IG' if igsensitive else None)
```



| dirs |
|----------|
#### A generic wrapper for quickly referencing all input and output directories used in analysis and data management.

```
# Example
from source.imports import *
dirs = io.dir_libraries
repo_root = dirs.Root #Repository root
datafolder = dirs.Data #Parent directory of all data inputs
plotfolder = dirs.Plots #Parent directory of all plot outputs
atacr_event_folder = dirs.Events #Event data from ATaCR
hps_event_folder = dirs.Events_HPS #Event data from NoiseCut
```

## All analysis in comparing noise reduction methodologies is completed sequentialy in five sections:


---
### **S00: Download events and remove noise.**
---

| S00.00_DownloadData.py |
|----------|
> ### A job script that tasks the obstools API with downloading events and noise data for use in NoiseCut and ATaCR.

| S00.01_Remove.Instrument.Response.from.Data.py |
|----------|
> ### Deconvolves the instrument response from all downloaded traces.

| S00.02_AutoRunATaCR.py |
|----------|
> ### A job script for removing noise from event data through transfer functions (TF) of cross-component coherent noise spectra using the ATaCR framework.

| S00.03_AutoRunNoiseCut.py |
|----------|
> ### A job script for removing noise from event data through Harmonic-Percussive Separation (HPS) using the NoiseCut framework.


---
### **S01: Collect signal measurements - coherence, admittance, phase, SNR, and SNR ratios.**
---
| S01.00_Collect.Event.Metrics.py |
|----------|
> ### Collects measurements of coherence, phase, and admittance between original and corrected traces using either HPS or TF.

| S01.01_Collect.SNR.py |
|----------|
> ### Collects measurements of SNR in user-defined windows around P/Pdiff, S/Sdiff, and Rayleigh surface waves before and after noise reduction using HPS or TF.

| S01.02_Record.Sections.w.SNR.py |
|----------|
> ### Detailed record sections showing the processing structure of SNR measurements. Namely phase and noise windows and resulting measurements for each event.

| S01.03_Collect.Noise.Metrics.py |
|----------|
> ### A detailed cross-component analysis of all noise data.

| S01.04.RunNoiseDataQualityReport.py |
|----------|
> ### Power-spectral densities of all noise data including quality controls from the ATaCR framework.
<details>
  <summary>(Example)</summary>
  <img src="_docs/S01.04.png" />
</details> 



---
### **S02: Data quality analysis based, automated QC, signal measurements, and record sections.**
---

| S02.00_RunStationPages.py |
|----------|
> ### Single station scatter plots of coherence spectra following corrections with HPS or TF.
<details>
  <summary>(Example)</summary>
  <img src="_docs/S02.00_7D.FS42D.both.Coherence.png"/>
</details> 

__
| S02.01_RunQuickHPSPlots.py |
|----------|

> ### Single single-receiver plots of traces after use of HPS to remove noise in all three seimic component data.
<details>
  <summary>(Example)</summary>
  <img src="_docs/S02.01_7D.FS42D.HPSEventCheck.2015.115.06.11.png"/>
</details> 

__
| S02.02_RunEventRecord_Trace_Sections.py |
|----------| 

> ### Record sections for the qualitative comparison of changes following noise reduction using HPS and TF.
<details>
  <summary>(Example)</summary>
  <img src="_docs/S02.02_2015.150.11.23.Mw7.9.685km_traces__NoiseCut.png"/>
</details> 


---
### **S03: Figures of final results used for publication**
---


| S03.Figure01_RunMapsPlot.py |
|----------|
> ### Global and regional maps of all deployments in the catalog.
<details>
  <summary>(Example)</summary>
  <img src="_docs/S03.Figure01.png" />
</details> 

__
| S03.Figure02_Methods.Example.py |
|----------|
> ### Directly comparing the two methods in a single example event using traces, spectrograms, and coherence.
<details>
  <summary>(Example)</summary>
  <img src="_docs/S03.Figure02.png" />
</details> 

__
| S03.Figure03_DirectMetricComparisons_byDepth_within.or.regardless.of.IG.py |
|----------|
> ### Scatter plots showing a direct comparison of average coherence and SNR ratios clustered by regular depth intervals both within and regardless of infragravity sensitivity.
<details>
  <summary>(Example)</summary>
  <img src="_docs/S03.Figure03.png" />
</details> 

__
| S03.Figure04_ConsolidatedPlot.py |
|----------|
> ### Consolidated overview plot of results.
<details>
  <summary>(Example)</summary>
  <img src="_docs/S03.Figure04.png" />
</details> 

__
| S03.Figure05_CoherenceContourPlots.py |
|----------|
> ### Contour plot of coherence measurement averages with depth.
<details>
  <summary>(Example)</summary>
  <img src="_docs/S03.Figure05.png" />
</details> 

__
| S03.Figure06_DirecComparisonAverages_by_DeploymentParams.py |
|----------|
> ### Scatter plots directly comparing coherences for each deployment parameter.
<details>
  <summary>(Example)</summary>
  <img src="_docs/S03.Figure06.png" />
</details> 

__
| S03.Figure07_Coherence.Spectra.by.DeploymentParam.py |
|----------|
> ### Multiple coherence spectra for each deployment and event parameter.
<details>
  <summary>(Example)</summary>
  <img src="_docs/S03.Figure07.png" />
</details> 

__
| S03.Figure08_SNR.SpectraBands.by.DeploymentParam.py |
|----------|
> ### Multiple snr ratio spectral averages for each deployment and event parameter.
<details>
  <summary>(Example)</summary>
  <img src="_docs/S03.Figure08.png" />
</details> 

__
| S03.Figure09_HighTilt.Plots.py |
|----------|
> ### Coherence and SNR ratio distributions grouped by severity of vertical tilt as approximated from ZH coherence averages.
<details>
  <summary>(Example)</summary>
  <img src="_docs/S03.Figure09.png" />
</details> 

__
| S03.Figure10_Coherence.with.HPS.Traces.Magnitude.py |
|----------|
> ### Spectra coherence between traces corrected with HPS and those corrected with TF.
<details>
  <summary>(Example)</summary>
  <img src="_docs/S03.Figure10.png" />
</details> 

__
| S03.Figure11_MetricComparisonScatterPlots.py |
|----------|
> ### Various scatter plots comparing different signal measurements beyond just coherence and snr (ie phase signal amplitude ratios).
<details>
  <summary>(Example)</summary>
  <img src="_docs/S03.Figure11.png" />
</details> 

__
| S03.Figure11.part2_ExampleEvents.with.LostStructure.py |
|----------|
> ### Example plot of selected events that demonstrate structure lost by noise reduction methods.
<details>
  <summary>(Example)</summary>
  <img src="_docs/S03.Figure11part2.png" />
</details> 

__
| S03.Table02_Correlograms.py |
|----------|
> ### Correlation tables of measurement and deployment parameters.
<details>
  <summary>(Example)</summary>
  <img src="_docs/S03.Table02.png" />
</details> 

__

---
### **S04: Supplemental figures**
---

| S04.FigureS01_MetaHistPlots.py |
|----------|
> ### Basic histograms of deployment and event parameters.
<details>
  <summary>(Example)</summary>
  <img src="_docs/S04.TableS01.png" />
</details> 

__
| S04.FigureS05_and_S08_CoherenceSpectraAverages.py |
|----------|
> ### Coherence and SNR ratio spectral averages with depth.
<details>
  <summary>(Example)</summary>
  <img src="_docs/S04.FigureS05S08.png" />
</details> 

__
| S04.FigureS06_ComparativeHistogram.py |
|----------|
> ### Histograms of measurement comparisons (ie which method produces lower coherence or snr on average).
<details>
  <summary>(Example)</summary>
  <img src="_docs/S04.FigureS06.png" />
</details> 

__
| S04.FigureS09_Histograms.py |
|----------|
> ### Histograms of coherence and snr ratio distrubutions for each deployment and event parameter.
<details>
  <summary>(Example)</summary>
  <img src="_docs/S04.FigureS09.png" />
</details> 

__
| S04.FigureS10_Coh.vs.Coh.Scatter.py |
|----------|
> ### Scatter plots directly comparing each method's coherence when grouped by different deployment and event parameters.
<details>
  <summary>(Example)</summary>
  <img src="_docs/S04.FigureS10.png" />
</details> 

__
| S04.FigureS13_CoherenceScatter_by_Meta.py |
|----------|
> ### Scatter plots of coherence and SNR ratio spectral averages with different deployment parameters
<details>
  <summary>(Example)</summary>
  <img src="_docs/S04.FigureS13.png" />
</details> 

__
| S04.FigureS14_NarrowSymmetryAnalysisPlots.py |
|----------|
> ### Scatter plots for simultaneous cross-measurement comparison of TF and HPS using coherence and snr ratios.
<details>
  <summary>(Example)</summary>
  <img src="_docs/S04.FigureS14.png" />
</details> 

__

## License 
##### With the exception of the ATaCR and NoiseCut scripts, this code was developed as part of the PhD research of Charles Hoots in the Department of Earth Sciences, University of Hawai‘i at Mānoa. Use of any analysis, data, or codes contained within is open-source, covered under the MIT License: https://github.com/charleshoots/ATACR_HPS_Comp/blob/main/LICENSE.txt

## References
Bell, S. W., D. W. Forsyth, and Y. Ruan (2014), Removing noise from the vertical component records of ocean-bottom seismometers: Results from year one of the Cascadia Initiative, Bull. Seismol. Soc. Am., 105, 300-313, https://doi.org/10.1785/0120140054

Crawford, W.C., Webb, S.C., (2000). Identifying and removing tilt noise from low-frequency (0.1 Hz) seafloor vertical seismic data, Bull. seism. Soc. Am., 90, 952-963, https://doi.org/10.1785/0119990121

Janiszewski, H A, J B Gaherty, G A Abers, H Gao, Z C Eilon, Amphibious surface-wave phase-velocity measurements of the Cascadia subduction zone, Geophysical Journal International, Volume 217, Issue 3, June 2019, Pages 1929-1948, https://doi.org/10.1093/gji/ggz051

Zali, Zahra, Theresa Rein, Frank Krüger, Matthias Ohrnberger, and Frank Scherbaum. “Ocean Bottom Seismometer (OBS) Noise Reduction from Horizontal and Vertical Components Using Harmonic–Percussive Separation Algorithms.” Solid Earth 14, no. 2 (2023): 181–95. https://doi.org/10.5194/se-14-181-2023.

NumPy – Harris, C. R. et al. (2020). Array programming with NumPy. Nature, 585, 357–362. https://doi.org/10.1038/s41586-020-2649-2

SciPy – Virtanen, P. et al. (2020). SciPy 1.0: Fundamental algorithms for scientific computing in Python. Nature Methods, 17, 261–272. https://doi.org/10.1038/s41592-019-0686-2

Pandas – McKinney, W. (2010). Data structures for statistical computing in Python. In Proceedings of the 9th Python in Science Conference (SciPy 2010), 51–56. https://doi.org/10.25080/Majora-92bf1922-00a

Librosa – McFee, B. et al. (2015). librosa: Audio and music signal analysis in Python. In Proceedings of the 14th Python in Science Conference (SciPy 2015), 18–25. https://doi.org/10.25080/Majora-7b98e3ed-003

ObsPy – Beyreuther, M., Barsch, R., Krischer, L., Megies, T., Behr, Y., & Wassermann, J. (2010). ObsPy: A Python toolbox for seismology/seismological observatories. Seismological Research Letters, 81(3), 530–533. https://doi.org/10.1785/gssrl.81.3.530

