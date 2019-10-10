from math import sqrt

def calcMoments(params,vS0m,vP0m,vP1m,vD0m,vD1m,vP1p,vD1p):
	nBins, numMoms, numAmps = params
	# This step trades readibility for generality but it is OK if we keep the order of amps the same
	allMoms_Mass = [[0 for col in range(nBins)] for row in range(numMoms)]
	for massNum in range(nBins):
		# ********* We will neglect which waveset the M=0 waves belongs to since it seems from the eq that the moments dont care *************
		S0 = vS0m[massNum]
		c_S0 = S0.conjugate()
		P0 = vP0m[massNum]
		c_P0 = P0.conjugate()
		P1m = vP1m[massNum]
		c_P1m = P1m.conjugate()
		D0 = vD0m[massNum]
		c_D0 = D0.conjugate()
		D1m = vD1m[massNum]
		c_D1m = D1m.conjugate()
		P1p = vP1p[massNum]
		c_P1p = P1p.conjugate()
		D1p = vD1p[massNum]
		c_D1p = D1p.conjugate()
		# calculating the moments 
		# "H00"
		allMoms_Mass[0][massNum] = (S0*c_S0+P0*c_P0+P1m*c_P1m+D0*c_D0+D1m*c_D1m+P1p*c_P1p+D1p*c_D1p).real
		#"H10"
		allMoms_Mass[1][massNum] = 2./sqrt(3.)*(S0*c_P0).real+4./sqrt(15.)*(P0*c_D0).real+2./sqrt(5.)*(P1m*c_D1m).real+2./sqrt(5.)*(P1p*c_D1p).real
		#"H11"
		allMoms_Mass[2][massNum] = 2./sqrt(6.)*(S0*c_P1m).real+2./sqrt(10.)*(P0*c_D1m).real-2./sqrt(30.)*(P1m*c_D0).real
		#"H20"
		allMoms_Mass[3][massNum] = 2./sqrt(5.)*(S0*c_D0).real+2./5.*(P0*c_P0).real-1./5.*(P1m*c_P1m).real-1./5.*(P1p*c_P1p).real + 2./7.*(D0*c_D0).real+1./7.*(D1m*c_D1m).real+1./7.*(D1p*c_D1p).real
		#"H21"
		allMoms_Mass[4][massNum] = 2./sqrt(10.)*(S0*c_D1m).real+2./5.*sqrt(3./2.)*(P0*c_P1m).real+2./7./sqrt(2.)*(D0*c_D1m).real
		#"H22"
		allMoms_Mass[5][massNum] = 1./5.*sqrt(3./2.)*(P1m*c_P1m).real-1./5.*sqrt(3./2.)*(P1p*c_P1p).real+1./7.*sqrt(3./2.)*(D1m*c_D1m).real-1./7.*sqrt(3./2.)*(D1p*c_D1p).real
		#"H30"
		allMoms_Mass[6][massNum] = 6./7.*sqrt(3./5.)*(P0*c_D0).real-6./7./sqrt(5.)*(P1m*c_D1m).real-6./7./sqrt(5.)*(P1p*c_D1p).real
		#"H31"
		allMoms_Mass[7][massNum] = 4./7.*sqrt(3./5.)*(P0*c_D1m).real+6./7./sqrt(5.)*(P1m*c_D0).real
		#"H32"
		allMoms_Mass[8][massNum] = 2./7.*sqrt(3./2.)*(P1m*c_D1m).real-2./7.*sqrt(3./2.)*(P1p*c_D1p).real
		#"H40" even though a term is like |D0|**2 we still have to take the real part since python will interpret it as a complex with 0j 
		allMoms_Mass[9][massNum] = 2./7.*(D0*c_D0).real-4./21.*(D1m*c_D1m).real-4./21.*(D1p*c_D1p).real
		#"H41"
		allMoms_Mass[10][massNum] = 2./7.*sqrt(5./3.)*(D0*c_D1m).real
		#"H42"
		#allMoms_Mass[11][massNum] = sqrt(10./21.)*(D1m*c_D1m).real-sqrt(10./21.)*(D1p*c_D1p).real
		allMoms_Mass[11][massNum] = 1./6.44*(D1m*c_D1m).real-1./6.44*(D1p*c_D1p).real
	return allMoms_Mass

# we will define the moments for the photon beam that vincent recalculated.
def calcMoments_ph(params,vS0,vP0,vP1,vD0,vD1,vD2):
	nBins, numMoms_ph, numAmps_ph = params 
	# This step trades readibility for generality but it is OK if we keep the order of amps the same
	allMom_Mass_ph = [0 for col in range(nBins)]
	allMoms_Mass_ph = [[0 for col in range(nBins)] for row in range(numMoms_ph)]
	for massNum in range(nBins):
		S0 = vS0[massNum]
		c_S0 = S0.conjugate()
		P0 = vP0[massNum]
		c_P0 = P0.conjugate()
		P1 = vP1[massNum]
		c_P1 = P1.conjugate()
		D0 = vD0[massNum]
		c_D0 = D0.conjugate()
		D1 = vD1[massNum]
		c_D1 = D1.conjugate()
		D2 = vD2[massNum]
		c_D2 = D2.conjugate()
		# We have to define the relevant quantity before using it obviously. That is why the order of these equations are done this way
		allMoms_Mass_ph[1][massNum]=-2.0*((S0*c_S0).real+(P0*c_P0).real+(D0*c_D0).real)
		allMoms_Mass_ph[0][massNum]=-allMoms_Mass_ph[1][massNum]+2*((P1*c_P1).real+(D1*c_D1).real+(D2*c_D2).real)
		allMoms_Mass_ph[3][massNum]=-4./5./sqrt(3.)*(2.*sqrt(5.)*(P0*c_D0).real+5*(S0*c_P0).real)
		allMoms_Mass_ph[2][massNum]=-allMoms_Mass_ph[3][massNum]+4./sqrt(5.)*(P1*c_D1).real
		allMoms_Mass_ph[5][massNum]=-2./15.*(3.*sqrt(5.)*(P0*c_D1).real-sqrt(15.)*(P1*c_D0).real+5.*sqrt(3.)*(S0*c_P1).real)
		allMoms_Mass_ph[4][massNum]=-allMoms_Mass_ph[5][massNum]+2.*sqrt(2./5.)*(P1*c_D2).real
		allMoms_Mass_ph[7][massNum]=-4./35.*(7.*(P0*c_P0).real+5.*(D0*c_D0).real+7.*sqrt(5.)*(S0*c_D0).real)
		allMoms_Mass_ph[6][massNum]=-allMoms_Mass_ph[7][massNum]-2./35.*(7.*(P1*c_P1).real-5.*(D1*c_D1).real+10.*(D2*c_D2).real)
		allMoms_Mass_ph[9][massNum]=-2./35.*(7.*sqrt(5.)*(S0*c_D1).real+7.*sqrt(3.)*(P0*c_P1).real+5*(D0*c_D1).real)
		allMoms_Mass_ph[8][massNum]=-allMoms_Mass_ph[9][massNum]+2./7.*sqrt(6.)*(D1*c_D2).real
		allMoms_Mass_ph[10][massNum]=2./35.*(7.*sqrt(5.)*(S0*c_D2).real-10.*(D0*c_D2).real)
		allMoms_Mass_ph[11][massNum]=-allMoms_Mass_ph[10][massNum]-sqrt(6.)/35.*(5.*(D1*c_D1).real+7.*(P1*c_P1).real)
		allMoms_Mass_ph[13][massNum]=-12./7.*sqrt(3./5.)*(P0*c_D0).real
		allMoms_Mass_ph[12][massNum]=-allMoms_Mass_ph[13][massNum]-12./7./sqrt(5.)*(P1*c_D1).real
		allMoms_Mass_ph[15][massNum]=-2./7.*sqrt(2./5.)*(2.*sqrt(3)*(P0*c_D1).real+3*(P1*c_D0).real)
		allMoms_Mass_ph[14][massNum]=-allMoms_Mass_ph[15][massNum]-2./7.*sqrt(3./5.)*(P1*c_D2).real
		allMoms_Mass_ph[17][massNum]=-2./7.*sqrt(3.)*((P0*c_D2).real+sqrt(2.)*(P1*c_D1).real)
		allMoms_Mass_ph[16][massNum]=-allMoms_Mass_ph[17][massNum]-2./7.*sqrt(6.)*(P1*c_D1).real
		allMoms_Mass_ph[19][massNum]=-6./7.*(P1*c_D2).real
		allMoms_Mass_ph[18][massNum]=0
		allMoms_Mass_ph[21][massNum]=-4./7.*(D0*c_D0).real
		allMoms_Mass_ph[20][massNum]=-allMoms_Mass_ph[21][massNum]-2./21.*(4.*(D1*c_D1).real-(D2*c_D2).real)
		allMoms_Mass_ph[23][massNum]=-2./7.*sqrt(10./3.)*(D0*c_D1).real
		allMoms_Mass_ph[22][massNum]=-allMoms_Mass_ph[23][massNum]-2./21.*sqrt(5.)*(D1*c_D2).real
		allMoms_Mass_ph[25][massNum]=-2./21.*sqrt(5.)*(sqrt(3.)*(D0*c_D2).real+sqrt(2)*(D1*c_D1).real)
		allMoms_Mass_ph[24][massNum]=-2./7.*sqrt(5./3.)*(D0*c_D2).real
		allMoms_Mass_ph[27][massNum]=-2./3.*sqrt(5./7.)*(D1*c_D2).real
		allMoms_Mass_ph[26][massNum]=0
		allMoms_Mass_ph[29][massNum]=-1./3.*sqrt(10./7.)*(D2*c_D2).real
		allMoms_Mass_ph[28][massNum]=0
	return allMoms_Mass_ph
