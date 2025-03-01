from imports import *
files = ['/Volumes/ECHO/Research/_DataArchive/HPS_Data/Data/YO.X09/2014.192.19.22.HZ.SAC',
'/Volumes/ECHO/Research/_DataArchive/HPS_Data/Data/YO.X09/2014.200.12.27.HZ.SAC',
'/Volumes/ECHO/Research/_DataArchive/HPS_Data/Data/YO.X09/2014.202.14.54.HZ.SAC',
'/Volumes/ECHO/Research/_DataArchive/HPS_Data/Data/7D.G34D/2015.144.04.53.HZ.SAC',
'/Volumes/ECHO/Research/_DataArchive/HPS_Data/Data/7D.G34D/2015.058.13.45.HZ.SAC',
'/Volumes/ECHO/Research/_DataArchive/HPS_Data/Data/7D.J28C/2013.289.10.31.HZ.SAC',
'/Volumes/ECHO/Research/_DataArchive/HPS_Data/Data/7D.J28C/2013.268.16.42.HZ.SAC',
'/Volumes/ECHO/Research/_DataArchive/HPS_Data/Data/7D.J28C/2013.304.12.02.HZ.SAC',
'/Volumes/ECHO/Research/_DataArchive/HPS_Data/Data/7D.J28C/2014.109.13.28.HZ.SAC',
'/Volumes/ECHO/Research/_DataArchive/HPS_Data/Data/7D.J42C/2014.109.13.28.HZ.SAC']
# st=Stream([load_sac(f)[0][0] for f in files])
st=Stream([read(f)[0] for f in files])
# [plt.plot(tr.times(),tr.data) for tr in st]
for tr in st:tr.plot()