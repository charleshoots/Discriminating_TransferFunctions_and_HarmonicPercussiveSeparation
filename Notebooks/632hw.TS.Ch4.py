# __________________________________________________________________
# QUESTION - 3
# __________________________________________________________________
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
from matplotlib import cm
# === PARAMETERS ===
Delta_T_max = 200  # °C
h_l = 80           # km (top of plume anomaly)
h_p = 80          # km (plume thickness)
kappa = 1e-6 * (60 * 60 * 24 * 365.25 * 1e6) / 1e6  # mm²/s → km²/Myr
k = 3.1            # W/m/K
# Depth vector
z = np.linspace(0, 500, 1000)
dz_km = z[1] - z[0]
dz_m = dz_km * 1e3  # meters
# === FUNCTION: Temperature profile ===
def T_profile(z, t):
    sqrt_kt = 2 * np.sqrt(kappa * t)
    term1 = erf((z - h_l) / sqrt_kt)
    term2 = erf((z + h_l) / sqrt_kt)
    term3 = erf((z - (h_l + h_p)) / sqrt_kt)
    term4 = erf((z + (h_l + h_p)) / sqrt_kt)
    return (Delta_T_max / 2) * ((term1 + term2) - (term3 + term4))
# === TEMPERATURE VS DEPTH: Times to highlight ===
tmax = 50
labelsize = 20
highlight_times = [0, 5, 15, 30, tmax]
colors = {0: "blue",5: "green",15: "red",30: "deepskyblue",tmax: "magenta"}
cmap=cm.get_cmap('turbo')
colors={k:cmap(k/55) for k in colors.keys()}
# Full background contour times
time_contours = np.linspace(0.01, tmax, 50)
T_matrix = np.array([T_profile(z, t) for t in time_contours])
# === HEAT FLUX VS DISTANCE ===
times_flux = np.linspace(0.01, 60, 300)   # Myr
times_flux = times_flux.tolist()
times_flux.extend([0,5,15,30,tmax])
times_flux = np.sort(np.unique(np.array(times_flux)))
distances_km = times_flux * 80           # 80 km/Myr
flux_mW_per_m2 = []
for t in times_flux:
    T = T_profile(z, t)
    dT_dz = np.gradient(T, dz_m)
    q_surface = -k * dT_dz[0] * 1e3  # W/m² → mW/m²
    flux_mW_per_m2.append(q_surface)
# === PLOTTING ===
titlesize = 20
fig, axs = plt.subplots(1, 2, figsize=(20, 8), sharey=False)
# -- LEFT: Temperature vs. Depth --
for T in T_matrix:
    axs[0].plot(T, -z, color='black', alpha=0.2, linewidth=0.5)

for t in highlight_times:
    T = T_profile(z, t if t > 0 else 1e-6)
    axs[0].plot(T, -z, label=f"{t} Ma, {str("%.0e" % (t * 80)).replace('+','') if t>0 else 0} km from Kilauea", color=colors[t], linewidth=2)
axs[0].set_xlim(0, 300)
axs[0].set_ylim(-500, 0)
axs[0].set_xlabel("ΔTemperature (°C)",fontsize=labelsize)
axs[0].set_ylabel("Depth (km)",fontsize=labelsize)
axs[0].set_title("$\Delta T_{max}$:" + str(Delta_T_max) + "°C"  "  |  Thermal Diffusion from Mantle Plume",fontsize=titlesize)
axs[0].legend(title="Time",fontsize=15)
axs[0].grid(True)

# -- RIGHT: Heat Flux vs. Distance --
axs[1].plot(distances_km, flux_mW_per_m2, label="Surface Heat Flux Anomaly")
highlight_flux = [flux_mW_per_m2[np.where(times_flux==k)[0][0]] for k in highlight_times]
_=[axs[1].scatter(xt*80,hf,color=colors[xt],s=150,zorder=100,edgecolor='k') for xt,hf in zip(highlight_times,highlight_flux)]
axs[1].axhline(-2, color="red", linestyle="--", label="Detection Threshold (-2 mW/m²)")

if np.array(flux_mW_per_m2).min()<= -2:
    intercept = round(distances_km[np.where(np.array(flux_mW_per_m2)<=-2)[0].min()])
    axs[1].axvline(intercept, color="k", linestyle="--", label=f"Detection Distance: {intercept} km")

axs[1].set_xlabel("Distance from Plume Center (km)",fontsize=labelsize)
axs[1].set_ylabel("Surface Heat Flux Anomaly (mW/m²)",fontsize=labelsize)
axs[1].set_title( "$\Delta T_{max}$:" + str(Delta_T_max) + "°C"  "  |  Heat Flux Anomaly vs Distance",fontsize=titlesize)
axs[1].legend(fontsize=13)
axs[1].grid(True)
axs[1].set_ylim(-3.8,.1)

ticksize=14
for i in [0,1]:axs[i].set_xticklabels(axs[i].get_xticklabels(),fontsize=ticksize)
for i in [0,1]:axs[i].set_yticklabels(axs[i].get_yticklabels(),fontsize=ticksize)

plt.tight_layout()
plt.show()


# __________________________________________________________________
# QUESTION - 4
# __________________________________________________________________

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
import pandas as pd

#Read data
df = pd.read_csv('/Users/charlesh/Downloads/Hw11_data.txt',delim_whitespace=True,names=['Age','Elevation'])

# Constants
alpha = 3.28e-5           # 1/K
Delta_T = 1300            # °C
kappa_m2_s = 1e-6         # m²/s
d0 = 2.6                  # km

# Time arrays
# t_myr = np.linspace(0.1, 160, 500)
t_myr = df.Age.values
t_sec = t_myr * 1e6 * 365.25 * 24 * 3600

# Plate cooling model function
def plate_model_depth(t_sec, yL0_m=95e3,Delta_T=1350, N_terms=200):
    sum_series = np.zeros_like(t_sec)
    for n in range(N_terms):
        term = (1 / (2*n + 1)) * np.exp(-((2*n + 1)**2 * np.pi**2 * kappa_m2_s * t_sec) / (4 * yL0_m**2))
        sum_series += term
    multiplier = (2 * Delta_T * alpha / np.pi)
    depth_km = d0 + multiplier * (yL0_m / 1e3) * (1 - (4/np.pi) * sum_series)
    return depth_km

plt.figure(figsize=(10, 7))
Delta_T = 1350
ylo = 95e3
# HSCM model
def hscm_depth(t_sec,Delta_T=1350):
    return d0 + 2 * alpha * Delta_T * np.sqrt(kappa_m2_s * t_sec / np.pi)
depth_HSCM = hscm_depth(t_sec,Delta_T)
Delta_T_list = [1350,1450]
ylo_list = [95e3, 130e3]


# These models seem to under-predict depths by around 200m. 
# See what happens with this adjustment?
adjust = 200


#There is steps in the data at ~90 and 105 Myr. Try t_reset's at those values.
t_reset = 105 #Myr
# t_reset = 91 #Myr

# ------Plot outputs------
t_reset_seconds = (t_reset* 1e6 * 365.25 * 24 * 3600) #convert reset to seconds
for Delta_T in Delta_T_list:
    for ylo in ylo_list:
        # Evaluate models
        depth_PM = plate_model_depth(t_sec + t_reset_seconds, yL0_m=ylo,Delta_T=Delta_T)
        # Plotting
        plt.plot(t_myr, (-depth_PM*1e3) - adjust, label='$\Delta T$'+f"{Delta_T}°C, yL0 = {int(ylo/1000)} km", linestyle='--', linewidth=2)
#Plot HSCM for reference
plt.plot(t_myr, -depth_HSCM, label="HSCM", color='black', linestyle='-', linewidth=2)
# Plot field data
plt.scatter(df.Age,(df.Elevation)*1000,c='k')
# Vertical line denotes t_reset age.
plt.axvline(t_reset, color='k', linestyle='--', label=f"Reset Age ({int(t_reset)}Myr)")
# Label stuff
plt.xlabel("Age, t (Myr)")
plt.ylabel("Depth of Seafloor (km)")
plt.title("Seafloor Depth Predictions: HSCM vs Plate Models")
plt.grid(True)
plt.legend()
plt.tight_layout()
