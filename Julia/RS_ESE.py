# Some constants and the planck function
import numpy as np
import netCDF4
import numpy as np

pi = np.pi
h = 6.626e-34
c = 299792458
k = 1.38e-23

def planck(wav, T):
    c1 = 2.0*h*c**2
    c2 = h*c/(wav*k*T)
    intensity = c1/ ( (wav**5)*(np.exp(c2) - 1.0) )

    return intensity*1.e-6

def getAtmosphere(myLat,myLon, url='files/MERRA300.prod.assim.inst6_3d_ana_Nv.20150613.hdf.nc4'):
	#url = 'data/MERRA300.prod.assim.inst6_3d_ana_Nv.20150613.hdf.nc4'
	file = netCDF4.Dataset(url)
	# What latitude do we want? Take Caltech as example
	#myLat = 34.1377
	#myLon = 118.1253

	lat  = file.variables['YDim'][:]
	lon  = file.variables['XDim'][:]
	# Temperature profile
	T =    file.variables['T'][:]
	# specific humidity profile
	q =    file.variables['QV'][:]
	# mean pressure profile:
	p = file.variables['Height'][:]
	# Surafce pressure
	psurf = file.variables['PS'][:]
	# Time in UTC
	time = file.variables['TIME'][:]

	# AK and BK global attributes (important to calculate pressure half-levels)
	ak = file.getncattr('HDF_GLOBAL.ak')[:]
	bk = file.getncattr('HDF_GLOBAL.bk')[:]
	
	iLat = np.argmin(np.abs(lat-myLat))
	iLon = np.argmin(np.abs(lon-myLon))

	# Let's choose one index in the time-domain
	index = 1

	# Surface pressure at Caltech and time slice 1 (scalar):
	ps_local = psurf[index,iLat,iLon]
	# q and T profiles at Caltech and time slice 1 (vector):
	q_local = q[index,:,iLat,iLon]
	T_local = T[index,:,iLat,iLon]

	# Half-level (at the boundaries, one index more), 
	# factor 100 is for conversion between Pa and hPa
	p_half = (ak + bk*ps_local)/100.
	# Full-levels, layer centered values:
	p_full = (p_half[0:-1]+p_half[1:])/2.
	NLEV = len(p_full)

	# Let us ignore gravity changes for simplicity
	go = 9.8196 #(m/s**2) 
	Rd = 287.04 # specific gas constant for dry air
	R_universal = 8.314472;
	Na = 6.0221415e23;
	    
	# array for dz
	dz = np.zeros((NLEV,))
	# array for molar density of dry air for each layer:
	rho_N = np.zeros((NLEV,))
	rho_N_h2o = np.zeros((NLEV,))
	# also get a VMR vector of H2O (volumetric!)
	vmr_h2o = np.zeros((NLEV,))
	# Now actually compute the layer thickness in meters
	for i in range(NLEV):
	    # Compute $\Delta$-z:
	    dz[i]    =  np.log(p_half[i+1]/p_half[i])*Rd*T_local[i]*(1+0.608*q_local[i])/go
	    # provide dry air density in molec/cm2/m (by diving molec/m3 by (100*100))
	    # For dry air, we actually have to subtract water vapor from the volumetric mixing ratio, i.e. (1-VMR_{volume}(H2O))
	    rho_N[i] =  p_full[i]*(1-q_local[i]*1.6068)*100./(R_universal*T_local[i])*Na/10000.0
	    # Do the equivalent for the water vapor contribution
	    rho_N_h2o[i] =  p_full[i]*(q_local[i]*1.6068)*100./(R_universal*T_local[i])*Na/10000.0
	    vmr_h2o[i] = q_local[i]*1.6068
	print('Total column density of dry air: ' +str(np.sum(dz*rho_N))+' molec/cm^2')
	print('Total column density of water vapor: ' + str(np.sum(dz*rho_N_h2o))+' molec/cm^2')
	VCD_dry = dz*rho_N
	return p_full, T_local, VCD_dry, NLEV, vmr_h2o,dz

def getAtmosphere400(myLat,myLon,index=1, url='../Week2/notebooks/Data/MERRA2_400.inst6_3d_ana_Nv.20150422.nc4'):
	#url = 'data/MERRA300.prod.assim.inst6_3d_ana_Nv.20150613.hdf.nc4'
	file = netCDF4.Dataset(url)

	lat  = file.variables['lat'][:]
	lon  = file.variables['lon'][:]
	# Temperature profile
	T =    file.variables['T'][:]
	# specific humidity profile
	q =    file.variables['QV'][:]
	# mean pressure profile:
	dp = file.variables['DELP'][:]
	#p = file.variables['Height'][:]
	# Surface pressure
	psurf = file.variables['PS'][:]
	# Time in UTC
	time = file.variables['time'][:]


	iLat = np.argmin(np.abs(lat-myLat))
	iLon = np.argmin(np.abs(lon-myLon))

	# Let's choose one index in the time-domain
	#index = 1

	# Surface pressure at Caltech and time slice 1 (scalar):
	ps_local = psurf[index,iLat,iLon]
	# q and T profiles at Caltech and time slice 1 (vector):
	q_local = q[index,:,iLat,iLon]
	T_local = T[index,:,iLat,iLon]
	dp_local = dp[index,:,iLat,iLon]
	print(dp_local)
	p_half = np.zeros((len(T_local)+1,))
	p_half[0]=0.01
	for i in range(len(T_local)):
	    p_half[i+1]=p_half[i]+dp_local[i]/100   
	# Half-level (at the boundaries, one index more), 
	# factor 100 is for conversion between Pa and hPa
	#p_half = (ak + bk*ps_local)/100.
	# Full-levels, layer centered values:
	p_full = (p_half[0:-1]+p_half[1:])/2.
	NLEV = len(p_full)

	# Let us ignore gravity changes for simplicity
	go = 9.8196 #(m/s**2) 
	Rd = 287.04 # specific gas constant for dry air
	R_universal = 8.314472;
	Na = 6.0221415e23;
	    
	# array for dz
	dz = np.zeros((NLEV,))
	# array for molar density of dry air for each layer:
	rho_N = np.zeros((NLEV,))
	rho_N_h2o = np.zeros((NLEV,))
	# also get a VMR vector of H2O (volumetric!)
	vmr_h2o = np.zeros((NLEV,))
	# Now actually compute the layer thickness in meters
	for i in range(NLEV):
	    # Compute $\Delta$-z:
	    dz[i]    =  np.log(p_half[i+1]/p_half[i])*Rd*T_local[i]*(1+0.608*q_local[i])/go
	    # provide dry air density in molec/cm2/m (by diving molec/m3 by (100*100))
	    # For dry air, we actually have to subtract water vapor from the volumetric mixing ratio, i.e. (1-VMR_{volume}(H2O))
	    rho_N[i] =  p_full[i]*(1-q_local[i]*1.6068)*100./(R_universal*T_local[i])*Na/10000.0
	    # Do the equivalent for the water vapor contribution
	    rho_N_h2o[i] =  p_full[i]*(q_local[i]*1.6068)*100./(R_universal*T_local[i])*Na/10000.0
	    vmr_h2o[i] = q_local[i]*1.6068
	print('Total column density of dry air: ' +str(np.sum(dz*rho_N))+' molec/cm^2')
	print('Total column density of water vapor: ' + str(np.sum(dz*rho_N_h2o))+' molec/cm^2')
	VCD_dry = dz*rho_N
	return p_full, T_local, VCD_dry, NLEV, vmr_h2o,dz