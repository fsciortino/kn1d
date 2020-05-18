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

def mProf(teArr,neArr, rArr):

    rMax = np.max(rArr)
    rMan = 2.34
    iMan = 0
    for i in range(len(rArr)):
        cR = rArr[i]
        if cR > rMan:
            if iMan == 0:
                iMan = i

            if teArr[i] < 10.0:
		teArr[i] = 10.0

	    neArr[i] = neArr[iMan]*math.exp(-100*(cR-rMan))
    return teArr,neArr



def test_kn1d():

    idl = pidly.IDL()
    stop = 0
    end =  130

    nE = np.load('python/nEnew.npy')
    
    rMid = nE[:,1]
    # tHis is a workaround because broadcasting causes some error
    for i in range(len(nE)):
        nE[i,0] = nE[i,0]*10**20
            
    nE = nE[:,0]
            
    # load pre-saved profiles
    tE = np.load('python/tEnew.npy')
    
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
    
    print('loaded')
    
    # Radial grid:
    rsep = 2.298
    xlim = 2.360
    wall = max(rMid)
    rKN1D = abs(rMid - wall)
    print(rKN1D)
    print(nE)
    print(tE)
    
    # debug plots
    fig,ax = plt.subplots(2,1)
    ax[0].plot(rKN1D,nE)
    ax[1].plot(rKN1D,tE)
    #plt.show()
    
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
    idl.gaugeH2 = 0.1 #mTorr
    idl.mu = 2 #deuterium 1 for hydrogen
    
    #manipulating outside the limiter
    connL = 0.5 #meters
    pipeDia = 0.02 #meters
    
    # set connection length in/out of limiter
    
    for i in range(len(idl.x)):
        if (idl.x[i] >= abs(xlim-wall))and (idl.x[i]<abs(rsep-wall)):
            idl('lc['+str(i)+'] = '+str(connL*2*math.pi*2.3))
        if (idl.x[i] >= abs(xlim-wall)):
            continue

        else:
            idl('lc['+str(i)+'] = '+str(connL))
            idl('lc['+str(i+1)+'] = '+str(connL))


    print(idl.lc)
    print(idl.d_pipe)
    
    print('Completed setup')
    
    idl('.r '+setPsPath)
    
    #we remove the .sav from the filename for saving the file
    filename = '"./savfiles/OpacTestOutputsIncGauge"'
    idl('set_ps, '+filename)
    
    #idl('save,x,x_lim,xsep,gaugeH2,mu,t_i,t_e,n_e,vx,lc,d_pipe,'+\
    #     'filename = "./savfiles/OpacTestOutputsIncGauge.sav"')
    
    #from IPython import embed
    #embed()
    
    # run KN1D:
    idl('kn1d,x,x_lim,xsep,gaugeH2,mu,t_i,t_e,n_e,vx,lc,d_pipe,'+\
        'xh2,nh2,gammaxh2,th2,qxh2_total,nhp,thp,sh,sp, xh,nh,'+\
        'gammaxh,th,qxh_total,nethsource,sion,qh_total,sidewallh,'+\
        'lyman,balmer,gammahlim,'+\
        'refine=1, File="pytest", NewFile=1, compute_errors=1,'+\
        'plot=1, debrief=1, pause=0, Hdebrief=1, H2debrief=1')
        #'/plot')
    #idl('lyman')
    idl('')
    
    print('finished running')
    
    #saving the lyman and blamer outputs for plotting later with x for plotting
    
    #idl('save,xH,nH,nHP,lyman,balmer,'+\
    #     'filename = "./savfiles/OpacTestOutputsIncGauge.sav"')
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
    #print(idl.xH)
    #print(idl.nH)
    #print(idl.xH2)
    #print(idl.nH2)
    
    ax[1].scatter(idl.xH2,idl.nH2)
    plt.show()
    
    import pickle
    
    with open("./savfiles/NeutProfGuage.dat", "wb") as f:
        pickle.dump(resDict, f)
            
    idl.close()

    return resDict


if __name__ == "__main__":
    res = test_kn1d()

    """
    import pickle
    with open("./savfiles/NeutProf.dat", "rb") as f:
        ddict = pickle.load(f)
    
    print(ddict.keys())
    print(ddict['midX'])
    """



