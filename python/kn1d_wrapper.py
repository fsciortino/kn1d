import pidly
import numpy as np
import pexpect
import matplotlib.pyplot as plt
import math

# Set_ps and set_x paths are here to be used in main so 
# that the kn1d output isn't just lost
# optionally could just remove idl.close() in main instead
# to just view the final output.

#CHANGE ME IF MOVED
setPsPath = '~/idl/set_ps.pro'
setXPath =  '~/idl/set_x.pro'



# We will set this if xExtend is called at some point
numXElem = 0

# function to extend the input arrays x, nE, tE and tI
# from 0 to the minimum value of x for passing to kn1d.
# Returns new x, Modifies the other input arrays.
# Asuumes elements of x are ascending
# Currently tE,tI and , nE are just constant extensions
# of their value at min(x). Maybe will make more physical
# eventually - July 2017
def xExtend(x):

	dx = 0.01   #seems small enough
	xmin = min(x)

	# We are filling in from 0 to min(x), so find number of 
	# elements required
	global numXElem 
	numXElem =  int(xmin / dx)

	x_extend = np.arange(numXElem)#map(lambda y: (y + 1) * dx, range(numXElem))
	vfunc = np.vectorize(lambda y: (y+1)*dx )
	x_extend=vfunc(x_extend)
	x = np.append(x_extend,x)
	return x


# adds an adElem number of elements to array of value array[0]
# returns the altered array
def extendArray(adElem,array):

	extend = np.full((1,adElem),array[0])
	array = np.append(extend,array)

	return array


def main():
	#initiate idl envi
	idl = pidly.IDL()

	idl('.r ~/idl/get_edgets.pro')
	shot = 1160718012 #raw_input('Shot Number:')
	tA = 1.0#raw_input('Start Time:')
	tB = 1.1 #raw_input('End Time:')

	filename = '"shot'+ str(shot) +'_'+str(tA).replace(".","")+'_' +str(tB).replace(".","") +'.sav"'
	print(filename)

	idl('get_edgets,'+str(shot)+','+str(tA)+','+str(tB)+
		',ts_ne=dens,ts_te=te,ts_rmid=rmid,r_lcfs=rsep')

	#load the cmod parameters in idl

	idl('restore, "../savfiles/kn1d_out_basecase_iter1.sav"')
	#now combine the required cmod parameters and those from the shot
	#time frame of interest



	idl.x = idl.rsep - idl.rmid +idl.x_sep
	idl.t_e = idl.te
	idl.t_i = idl.te
	idl.n_e = idl.dens
	idl('vx = findgen(n_elements(t_e))*0.0')
	idl('d_pipe = findgen(n_elements(t_e))*0.0')
	idl('lc = findgen(n_elements(t_e))*0.0')
	
	#alright now we have to reverse all these elements which we will use
	
	idl('x = rotate(x,2)')
	idl('t_e = rotate(t_e,2)')
	idl('t_i = rotate(t_i,2)')
	idl('n_e = rotate(n_e,2)')

	

	
	
	#get the filename and load it
	#filename = raw_input('.sav filename: ')
	idl('save,x,x_lim,x_sep,p_wall,mu,t_i,t_e,n_e,vx,lc,d_pipe,filename = '+filename)
	#idl('restore,"'+ filename +'.sav"')

	#extract the arrays, assuming convential naming scheme

	extend = raw_input('Extend arrays to 0:')
	if (extend == True) or ('y' == extend.lower()) or ('yes' == extend.lower()):
		#unpack the arrays we will be manipulating
		#these are numpy narrays so ...
		nE = idl.n_e
		print(nE)
		x = idl.x
		print(x)
		tI = idl.t_I
		print(tI)
		tE = idl.t_e
		vx = idl.vx
		lc = idl.lc
		dPipe = idl.d_pipe
	

		#manipulate away


		# now maniuplate them and reset the IDL array in IDL to avoid 
		# passing large arrays to idl
		idl.x = xExtend(x)
		idl.d_pipe = extendArray(numXElem,dPipe)
		idl.n_e = extendArray(numXElem,nE)
		print(idl.n_e)
		idl.t_e = extendArray(numXElem,tE)
		idl.t_i = extendArray(numXElem,tI)
		idl.vx = extendArray(numXElem,vx)
		idl.lc = extendArray(numXElem,lc)

			#now we are gonna add in a connection length outside the limiter

		for i in range(len(idl.t_e)):
			if (idl.x[i] >=  idl.x_lim):
				break
			else:
				idl('lc['+str(i)+'] = 0.8')

		
	idl('.r '+setPsPath)
	
	#we remove the .sav from the filename for saving the file
	idl('set_ps, '+filename)

	idl('kn1d,x,x_lim,x_sep,p_wall,mu,t_i,t_e,n_e,vx,lc,d_pipe,xH,nH,lyman, balmer,/plot')
	idl('lyman')
	idl('')
	#saving the lyman and blamer outputs for plotting later with x for plotting
	filename = '"lb'+ str(shot) +'_'+str(tA).replace(".","")+'_' +str(tB).replace(".","") +'.sav"'
	idl('save,xH,nH,lyman,balmer,filename = '+ filename)

	# finally we close idl

	idl.close()
	


	return


def mProf(teArr,neArr, rArr):

	rMax = np.max(rArr)

	print('rMAx: '+str(rMax))

	rMan = 2.34

	"""
	rExt = 2.2
	rArray = np.linspace(rMax,rExt,10)

	teAdd = np.zeros(len(rArray))

	"""
	iMan = 0
	for i in range(len(rArr)):

		cR = rArr[i]
		if cR > rMan:
			if iMan == 0:
				iMan = i

			if teArr[i] < 10.0:
				teArr[i] = 10.0


			#teArr[i] = teArr[iMan]*math.exp(-2*(cR-rMan))
			neArr[i] = neArr[iMan]*math.exp(-100*(cR-rMan))




	

	return teArr,neArr



def opacTest():

	idl = pidly.IDL()
	stop = 0
	end =  130

	filename = '"../savfiles/OpacTestOutputsIncGauge"'


	nE = np.load('nEnew.npy')

	rMid = nE[:,1]
	# tHis is a workaround because broadcasting causes some error
	for i in range(len(nE)):
		nE[i,0] = nE[i,0]*10**20

	nE = nE[:,0]
	



	tE = np.load('tEnew.npy')

	rMidTE = tE[:,1]
	tE = tE[:,0]*1000




	#remove some of the core
	nE = nE[stop:end]
	tE = tE[stop:end]
	rMid = rMid[stop:end]

	tE,nE = mProf(tE,nE,rMid)


	#flip the arrays
	nE = nE[::-1]
	tE = tE[::-1]
	rMid = rMid[::-1]







	tI = tE

	vx = np.zeros(len(tE))
	#dPipe = np.full(len(tE),0.02)
	#lc = np.full(len(tE),6)
	dPipe =np.zeros(len(tE))
	lc = np.zeros(len(tE))
	




	"""

	np.savetxt('kn1dFiles/lc.txt',lc)
	np.savetxt('kn1dFiles/vx.txt',vx)
	np.savetxt('kn1dFiles/dPipe.txt',dPipe)

	np.savetxt('kn1dFiles/nE.txt',nE)
	np.savetxt('kn1dFiles/tE.txt',tE)
	np.savetxt('kn1dFiles/tI.txt',tI)

	np.savetxt('kn1dFiles/nRho.txt',nRho)

	

	"""
	print('loaded')

	#idl.x = idl.rsep - idl.rmid + idl.x_sep
	
	print('opened idl')

	#manipulate away

	rsep = 2.298
	xlim = 2.360

	wall = max(rMid)

	rKN1D = abs(rMid - wall)
	print(rKN1D)
	print(nE)
	print(tE)
	"""
	fig,ax = plt.subplots(2,1)
	ax[0].plot(rKN1D,nE)
	ax[1].plot(rKN1D,tE)
	#plt.show()
	"""


	# now maniuplate them and reset the IDL array in IDL to avoid 
	# passing large arrays to idl
	idl.x = rKN1D 
	idl.d_pipe = dPipe
	print(nE.shape)
	idl.n_e = nE[:]
	idl.xsep = abs(rsep-wall)

	idl.t_e = tE
	idl.t_i = tI
	idl.vx = vx
	idl.lc = lc

	idl.x_lim = abs(xlim - wall)
	idl.gaugeH2 = 1 #mTorr
	idl.mu = 2.00 #deuterium 1 for hydrogen


	#manipulating outside the limiter
	connL = 0.5 #meters
	pipeDia = 0.02 #meters

	#idl.d_pipe[0] = pipeDia


	#now we are gonna add in a connection length outside the limiter

	for i in range(len(idl.x)):
		"""
		if (idl.x[i] >= abs(xlim-wall)):
			continue
		"""
		if (idl.x[i] >= abs(xlim-wall))and (idl.x[i]<abs(rsep-wall)):
			idl('lc['+str(i)+'] = '+str(connL*2*math.pi*2.3))
		if (idl.x[i] >= abs(xlim-wall)):
			continue

		else:
			idl('lc['+str(i)+'] = '+str(connL))
			idl('lc['+str(i+1)+'] = '+str(connL))
			#idl('d_pipe['+str(i)+'] = '+str(pipeDia))
			#idl('d_pipe['+str(i+1)+'] = '+str(pipeDia))

	print(idl.lc)
	print(idl.d_pipe)



	print('all prints done')

		
	idl('.r '+setPsPath)
	
	#we remove the .sav from the filename for saving the file
	idl('set_ps, '+filename)

	idl('save,x,x_lim,xsep,gaugeH2,mu,t_i,t_e,n_e,vx,lc,d_pipe,filename = "../savfiles/OpacTestOutputsIncGauge.sav"')

	#idl('kn1d,x,x_lim,xsep,gaugeH2,mu,t_i,t_e,n_e,vx,lc,d_pipe,xH2,nH2, nHP,xH,nH,lyman, balmer,/plot')
	idl('kn1d,x,x_lim,xsep,gaugeH2,mu,t_i,t_e,n_e,vx,lc,d_pipe, xh2,nh2,gammaxh2,th2,qxh2_total,nhp,thp,sh,sp, xh,nh,gammaxh,th,qxh_total,nethsource,sion,qh_total,sidewallh,lyman,balmer,gammahlim,/plot')
	idl('lyman')
	idl('')

	print('finished running')

	#saving the lyman and blamer outputs for plotting later with x for plotting

	idl('save,xH,nH,nHP,lyman,balmer,filename = "../savfiles/OpacTestOutputsIncGauge.sav"')
	"""
	print(idl.xH)
	print(idl.nH)
	print(idl.lyman)
	print(idl.balmer)
	print(idl.nHP)
	print(idl.xH2)
	"""
	# finally we close idl


	x = wall - idl.xH


	resDict = {}

	resDict['nH'] = idl.nH
	resDict['midX'] = x
	resDict['lyman'] = idl.lyman
	resDict['balmer'] = idl.balmer
	resDict['tH'] = idl.th

	subDict = {}
	subDict['rmid'] = rMid
	subDict['te'] = tE
	subDict['ne'] = nE
	subDict['ti'] = tI

	resDict['inputs'] = subDict

	fig,ax = plt.subplots(2,1)
	ax[0].scatter(idl.xH,idl.nH)
	print(idl.xH)
	print(idl.nH)
	print(idl.xH2)
	print(idl.nH2)

	ax[1].scatter(idl.xH2,idl.nH2)
	plt.show()

	import pickle

	with open("../savfiles/NeutProfGuage.dat", "wb") as f:
		pickle.dump(resDict, f)

	idl.close()
	
	

	return 0 

if __name__ == "__main__":
	#main()


	opacTest()
	"""
	import pickle
	with open("../savfiles/NeutProf.dat", "rb") as f:
		ddict = pickle.load(f)

	print(ddict.keys())
	print(ddict['midX'])
	"""



