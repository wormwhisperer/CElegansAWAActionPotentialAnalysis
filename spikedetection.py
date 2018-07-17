from scipy import *
from pylab import *
import os
import cPickle
from scipy.ndimage import filters
from scipy.signal import gaussian
ion()


##### first some random utility functions

def diffthresh(ts,xs): # this function actually identifies individual spikes

    dt = ts[1]-ts[0] ## time bin
    nframes = int(round(0.04/dt)) #0.04 for wt, 0.06 for shk, size of averaging window for calculating derivative
    winspace = int(round(0.005/dt)) # spacing to check for derivatives above threshhold
    initpoint = int(round(0.9/dt)) # start time in seconds
    endpoint = int(round(6.0/dt)) # end time
    
    ns = (endpoint-initpoint)
    
    dthresh = 0.006 #0.006 for wt, 0.0085 for shk, size for rising derivative threshhold
    dthresh2 = 0.004  # 0.004 for wt, 0.006 for shk, size for falling derivative threshhold
    vthresh = -0.03  # overall voltage threshhold
    upcheck = list(ones([winspace,])-1)  # index vector for upswings
    downcheck = list(ones([winspace,])-1) # index vector for downswings
    indicator = 0
    
    point = initpoint
    
    while point<endpoint-nframes-2:
        diff = xs[point+nframes]-xs[point]
        if (diff>dthresh2)&(max(xs[point:point+nframes])>vthresh):
            upcheck.extend(list(ones([winspace,])))
            indicator = 1
        else:
            upcheck.extend(list(ones([winspace,])-1))
        if (diff<-dthresh):
            downcheck.extend(list(ones([winspace,])))
        else:
            downcheck.extend(list(ones([winspace,])-1))
        if indicator&(downcheck[-1]-downcheck[-winspace-1]==-1):
            indicator = 0
        point += winspace
    
            
    upcheck = upcheck[winspace:min(ns+winspace,len(upcheck)-1)]
    downcheck = downcheck[winspace:min(ns+winspace,len(downcheck)-1)]
    

    
    if ns>len(upcheck):
        diff = ns-len(upcheck)
        upcheck.extend(list(ones([diff,])-1))
        downcheck.extend(list(ones([diff,])-1))
        

    return upcheck,downcheck,ts[initpoint:endpoint],xs[initpoint:endpoint],initpoint

def tabdata(filename,n,char='\t'):
    
    strn = open(filename,'r')
    outd = {}
    for i in range(n):
        outd[i] = []
        
    
    for ln in strn.readlines():
        ln2 = ln.strip('\n')
        vec = ln2.split(char)
        for i in range(n):
            try:
                outd[i].append(float(vec[i]))
            except:
                pass
        
    return outd

def jumps(vec):  ## function for finding the indices of up and downswings
    
    dupc = array(vec[1:])-array(vec[0:-1])
    upinds = nonzero(abs(dupc))[0]
    
    return upinds

### now the core spikedetection functions


def analysePlot(dic):  ## this function plots the first order identification of upswings and downswings
    
    ### input dic is the traces dictionary returned by readData
    
    num = 0
    for j in dic.keys():
        num+=len(dic[j][0,:])-1
    
    fig = figure()
    axd = {}
    n1 = int(sqrt(num))
    n2 = int(ceil((num+0.0)/n1))
    c1 = 1
    c2 = 0
    outd = {}
    n = 0
    
    chk = 0
    for i in range(n2):
        for j in range(n1):
            if chk<num:
                axd[n1*i+j+1] = fig.add_subplot(n2,n1,n1*i+j+1)
                chk+=1
    for k in range(1,num+1):
        if c1>len(dic[c2][0,:])-1:
            c1=1
            c2+=1
        upc,downc,ts,xs,initpoint = diffthresh(dic[c2][:,0],dic[c2][:,c1])
        outd[n] = dic[c2],upc,downc,ts,xs
        n+=1
        xs = 1000*array(xs)
        axd[k].plot(ts,xs,'-')
        ups = nonzero(upc)[0]
        downs = nonzero(downc)[0]
        tsu = [ts[i] for i in ups]
        xsu = [xs[i] for i in ups]
        tsd = [ts[i] for i in downs]
        xsd = [xs[i] for i in downs]
        axd[k].plot(tsu,xsu,'.')
        axd[k].plot(tsd,xsd,'.',markersize=1)
        c1+=1
        print c1,c2
    
    return fig,axd,outd

def cleanTimes2(upc,downc,ts,xs,initpoint): #### this function cleans up the initial guess of spikes
    
    ### input dic is the traces dictionary returned by readData
    
    upinds = jumps(upc)
    downinds = jumps(downc)
    upi = []
    di = []
    spiketimes = []
    dt = ts[1]-ts[0]
    spikewidth = int(round(0.02/dt))  ## 0.5 for shk, 0.1 for wt, approx width of a spike 
    spikedic = {}
    
    
    
    check=1
    
    while check:  ## check for consecutive upswing/downswing alternation within spec'd time limits
        if len(upinds)*len(downinds)==0:
            print 'breaking'
            break
        check1 = 1
        while check1:
            if downinds[0]<upinds[0]:
                print 'chopping initial down'
                downinds = downinds[2:]
            else:
                check1 = 0
        if len(upinds)<4:
            print 'adding to uplist'
            uplist = list(upinds)
            uplist.append(len(xs)-2)
            uplist.append(len(xs)-1)
            upinds = array(uplist)
            check=0
        if len(downinds)<4:
            print 'adding to downlist'
            downlist = list(downinds)
            downlist.append(len(xs)-2)
            downlist.append(len(xs)-1)
            downinds = array(downlist)
            #check=0
        check1 = 1
        while check1:
            if (upinds[3]<downinds[0])*(upinds[2]-upinds[1]<spikewidth):
                upinds = delete(upinds,[1,2])
                print 'joining ups'
                if len(upinds)<4:
                    print 'adding to uplist'
                    uplist = list(upinds)
                    uplist.append(len(xs)-2)
                    uplist.append(len(xs)-1)
                    upinds = array(uplist)
                    check=0
            else:
                check1 = 0
        check1 = 1
        while check1:
            if (downinds[3]<upinds[2])*(downinds[2]-downinds[1]<spikewidth):
                downinds = delete(downinds,[1,2])
                print 'joining downs'
                if len(downinds)<4:
                    print 'adding to downlist'
                    downlist = list(downinds)
                    downlist.append(len(xs)-2)
                    downlist.append(len(xs)-1)
                    downinds = array(downlist)
            else:
                check1 = 0
        if (upinds[2]>downinds[1])*((downinds[0]-upinds[1])<spikewidth):
            if max(xs[upinds[0]:downinds[1]])>-0.026:
                print 'real spike found!',upinds[2]-downinds[1]
                upi.append(upinds[0])
                upinds = upinds[2:]
                di.append(downinds[1])
                downinds = downinds[2:]
            else:
                print 'below threshold'
                upinds = upinds[2:]
                downinds = downinds[2:]
        else:
            print 'chopping fake initial upswing'
            upinds = upinds[2:]
            
    print 'leaving loop',len(di),len(upi)
      
    if len(di):
        if max(di)>len(ts):
            print 'wtf',max(di)
    if len(upi):
        if max(upi)>len(ts):
            print 'wtf',max(upi)
            
    ts2 = []
    xs2 = []
    j = 0
    
    for i in range(len(upi)):
        ts2.extend(ts[upi[i]:di[i]])
        xs2.extend(xs[upi[i]:di[i]])
        ### these lines define what range is included in a spike for averaging        
        timebefore = 0.2 ## time before point of max slope
        timeafter = 0.15 ## time after point of max slope
        ### this code finds the point of max slope
        wind = int(round(2e-2/dt))
        c = gaussian(wind,wind/2)
        smoothed = filters.convolve1d(xs[upi[i]:di[i]],c/c.sum())
        ds = smoothed[1:]-smoothed[0:-1]
        ind = argmax(ds)+upi[i]
        
        spikedic[j] = ts[ind-int(round(timebefore/dt)):ind+int(round(timeafter/dt))],xs[ind-int(round(timebefore/dt)):ind+int(round(timeafter/dt))]
        j+=1
    
    spiketimes = 0.5*dt*(array(upi)+array(di))+dt*initpoint
    
    ### outputs are indices of upswings, downswings, and the time and potential vectors for each spike
    ### the list of spike times, and a dictionary with the traces of each spike
    
    return upi,di,ts2,xs2,spiketimes,spikedic

def cleanPlot2(dic):  ### this function plots the final identified spikes
    
    ### input dic is the traces dictionary returned by readData
    
    num = 0
    for j in dic.keys():
        num+=len(dic[j][0,:])-1
    
    fig = figure()
    axd = {}
    n1 = int(sqrt(num))
    n2 = int(ceil((num+0.0)/n1))
    c1 = 1
    c2 = 0
    outd = {}
    n = 0
    spikesdic = {}
    firstspikedic = {}
    firstspiketimes = []
    freqs = []
    jj = 0
    
    for i in range(n2):
        for j in range(n1):
            axd[n1*i+j+1] = fig.add_subplot(n2,n1,n1*i+j+1)
    for k in range(1,num+1):
        if c1>len(dic[c2][0,:])-1:
            c1=1
            c2+=1
        upc,downc,ts,xs,initpoint = diffthresh(dic[c2][:,0],dic[c2][:,c1])
        upi,di,ts2,xs2,spiketimes,spikedic = cleanTimes2(upc,downc,ts,xs,initpoint)
        testt = 1e6
        for kk in spikedic.keys():
            tsss,xsss = spikedic[kk]
            if len(tsss)>0:
                if min(tsss)<testt:
                    indo = kk
                    testt = min(tsss)
        if indo in spikedic.keys():
            tsss,xsss = spikedic[indo]
            if len (tsss)>0:
                firstspikedic[jj] = spikedic[indo]
                jj+=1
                tsss,xsss = spikedic[indo]
                del spikedic[indo]
        if len(spiketimes)>0:
            firstspiketimes.append(min(spiketimes))
        if len(spiketimes)>1:
            freqs.append(freqCalc(spiketimes))    
        outd[n] = dic[c2],upi,di,ts2,xs2,spiketimes
        if len(spikesdic.keys())>0:
            last = spikesdic.keys()[-1]
        else:
            last = 0
        for q in spikedic.keys():
            spikesdic[last+1+q] = spikedic[q]
        n+=1
        axd[k].plot(ts,xs,'-')
        axd[k].plot(ts2,xs2,'.',markersize=1)
        axd[k].plot(tsss,xsss,'m.',markersize=1)
        #axd[k].plot(ts2,array(xs)*array(downc),'.')
        c1+=1
        print k
    
    allspiketimes = []
    for j in outd.keys():
        allspiketimes.extend(outd[j][5])
    
        
    ### output is the figure object, a dictionary of axes, the dictionary of individual spikes
    ### and a list of spike times
    
    return fig,axd,outd,spikesdic,allspiketimes,firstspiketimes,freqs,firstspikedic

#### now some plotting functions


def plotAllSpikes(spikesdic):
    for j in spikesdic.keys():
        ts,xs=spikesdic[j]
        if len(ts)<10:
            del spikesdic[j]
    fig = figure()
    ax = fig.add_subplot(111)
    for j in spikesdic.keys():
        ts,xs = spikesdic[j]
        ax.plot(array(ts)-ts[0],xs)
        #ax.plot(ts,xs) ## use this one to plot spikes at their original times
    
    return

### now functions for reading in data

def readData(folder):
    
    ### input is the full path to a folder of data, eg '/Users/pkidd/Documents/data/slo-2'
    
    lst = os.listdir(folder)
    datatrace = {}
    k = 0
    for i in lst:
        if (i.find('.txt')>0):
            dic = tabdata(folder+'/'+i,50)
            fig = figure()
            axd = {}
            j = 1
            check = 1
            while check:
                try:
                    axd[j] = fig.add_subplot(4,5,j)
                    axd[j].plot(dic[0],dic[j])
                    j+=1
                except:
                    check = 0
            pause(0.05)
            axlst = input("Enter list of traces with spikes, eg [5,6,7]")
            close("all")
            dat = ones([len(dic[0]),len(axlst)+1])-1
            dat[:,0] = dic[0]
            for q in range(len(axlst)):
                dat[:,q+1] = dic[axlst[q]]
            if dat.shape[1]:
                datatrace[k] = dat
                k+=1
    
    f = file(folder+'datatrace','wb')
    cPickle.dump(datatrace,f,protocol=cPickle.HIGHEST_PROTOCOL)
    f.close()
    
    ### output is the dictionary of traces, goes into cleanPlot
    
    return datatrace

def importList(filename):
    
    spikeslist = {}
    f = file(filename,'rb')
    lst = f.readlines()
    first = lst[0]
    first.strip('\n')
    first.replace('\r','',1)
    cols = first.split('\t')
    num = len(cols)/2
    f.close()
    
    for j in range(num-1):
        ts = []
        xs = []
        for ln in lst[0:-3]:
            ln3 = ln.strip('\n')
            ln4 = ln3.replace('\r','')
            ln2 = ln4.split('\t')
            if ln2[2*j]!='' and ln2[2*j+1]!='' and ln2[2*j]!=' ' and ln2[2*j+1]!=' ':
                ts.append(float(ln2[2*j]))
                xs.append(float(ln2[2*j+1]))
        spikeslist[j] = array(ts),array(xs)
    
    return spikeslist


def fixFiles(filename):
    
    f = file(filename,'rb')
    txt = f.read()
    f.close()
    check  = 0
    j = 0
    total = len(txt)
    
    while j<total:
        if txt[j]=='.' and check==1:
            txt = txt[0:j-1]+'\t'+txt[j-1:]
            total = total+1
            check = 0
        if txt[j]=='.':
            check = 1
        if txt[j]=='\t' or txt[j]=='\n':
            check = 0
        j = j+1
        
    f = file(filename,'wb')
    f.write(txt)
    f.close()
    
    return

def spikeread(filename):
    
    strn = open(filename,'r')
    times=  []
    spike = []
    errs = []
    
    for ln in strn.readlines():
        ln2 = ln.replace('\n','')
        ln3 = ln2.replace('\r','')
        vec = ln3.split('\t')
        if len(vec)>1:
            times.append(float(vec[0]))
            spike.append(float(vec[1]))
            errs.append(float(vec[2]))
        
    
    return times,spike,errs

#### now functions for writing data

def writestuff(urspike,ts,errs,allspiketimes,firstspiketimes,prefix): ### this function writes data
    
    ### inputs are, the average spike, the time points, the standard error, the spike times list, and a prefix label
    
    f = file(prefix+'_allspiketimes.txt','wb')
    for i in range(len(allspiketimes)):
        f.write(str(allspiketimes[i])+'\n\r')
    f.close()
    
    f = file(prefix+'_firstspiketimes.txt','wb')
    for i in range(len(firstspiketimes)):
        f.write(str(firstspiketimes[i])+'\n\r')
    f.close()
    
    f = file(prefix+'_average.txt','wb')
    for i in range(len(urspike)):
        f.write(str(ts[i])+'\t'+str(urspike[i])+'\t'+str(errs[i])+'\n\r')
    f.close
    
    return



def writespikes(spikedic,filename,dec):
    
    f = file(filename,'wb')
    
    lst = spikedic.keys()
    lens = []
    for j in lst:
        lens.append(len(spikedic[j][0]))
    mx = max(lens)
    
    for i in range(mx/int(dec)):
        for j in lst:
            ts,xs = spikedic[j]
            if len(ts)>0:
                ts = array(ts)-ts[0]
                if int(dec)*i<len(ts)-1:
                    f.write(str(ts[int(dec)*i])+'\t'+str(xs[int(dec)*i]))
                else:
                    f.write(' '+'\t'+' ')
                if j<len(lst)+1:
                    f.write('\t')
        f.write('\n\r')
    
    f.close()
    
    return


#### now analysis functions, spike averaging, etc

def urspikeSlope(spkdic):  ### this function calculates the average spike
    
    ### input is the spikesdic from cleanPlot
    
    b = gaussian(100,20)
    c = gaussian(200,40)
    
    endtime = 1 # 0.26 for WT, 0.6 for SHK. this excludes spikes longer than a certain time, helps to eliminate some odd mistakes
    dts = []
    ets = []
    bts = []
    for j in spkdic.keys():        
        ts,xs = spkdic[j]
        if len(ts)>10 and (ts[-1]-ts[0])<endtime:
            dts.append(ts[1]-ts[0])
        else:
            del spkdic[j]
    dt = max(dts)
    for j in spkdic.keys():
        ts,xs = spkdic[j]
        ts = array(ts)-ts[0]
        decfac = int(round(dt/ts[1]))
        newlen = int(round(len(xs)*ts[1]/dt))
        decspike = []
        dects = []
        for k in range(newlen):
            decspike.append(xs[decfac*k])
            dects.append(k*dt)
        spkdic[j] = dects,decspike
        ga = filters.convolve1d(decspike,b/b.sum())
        ga2 = filters.convolve1d(decspike,c/c.sum())
        dga = list(ga2[1:]-ga2[0:-1])
        dga.append(0)
        dga = array(dga)
        ind = argmax(ga)
        ind2 = argmax(dga[0:ind])
        ets.append(dects[-1]-dects[ind2])
        bts.append(dects[ind2]-dects[0])
    ent = max(ets)
    bnt = max(bts)
    ln = int(round((ent+bnt)/dt))+2
    spike = ones([ln,])-1
    num = 0
    newdic = {}
    tims = linspace(0,dt*ln,ln)
    for j in spkdic.keys():
        ts,xs = spkdic[j]
        ga = filters.convolve1d(xs,b/b.sum())
        ga2 = filters.convolve1d(xs,c/c.sum())
        dga = list(ga2[1:]-ga2[0:-1])
        dga.append(0)
        dga = array(dga)
        ind = argmax(ga)
        maxpt = argmax(dga[0:ind])
        xs2 = list(xs[0]*ones([int(round(bnt/dt))-maxpt,]))
        xs2.extend(list(xs))
        xs2.extend(list(xs[-1]*ones([ln-len(xs2),])))
        newdic[num] = tims,xs2
        num+=1
    
    for j in newdic.keys():
        ts,xs = newdic[j]
        spike = spike+array(xs)/num
    
    print num
    
    stanD = ones([ln,])-1
    
    for j in newdic.keys():
        ts,xs = newdic[j]
        stanD = stanD + sqrt((spike-xs)**2/num**2)
        
    plotAllSpikes(newdic)
            
    
    ### output is the average spike, the time points, and the standard error
    
    return spike,tims,stanD





def spikeProps(times,spike,errs):
    
    dt = times[1]-times[0]
    wind = int(round(2e-3/dt))
    c = gaussian(wind,wind/2)
    smoothed = filters.convolve1d(spike,c/c.sum())
    errsS = filters.convolve1d(errs,c/c.sum())
    ind = argmax(smoothed)
    errsL = []
    
    deriv = smoothed[1:]-smoothed[0:-1]
    thresh = 0.03*max(deriv)   ### the prefactor is the proportion of the maximum derivative at which the threshold sits
    
    ind2 = argmax(deriv)
    segment = deriv[0:ind2]
    flip = segment[::-1]-thresh
    flip = flip/abs(flip)
    thpoint = ind2-argmin(flip)
    errsL.append(errsS[thpoint])
    
    
    amplitude = max(spike)-smoothed[thpoint]
    errsL.append(sqrt(errsS[thpoint]**2+errsS[argmax(spike)]**2))
    
    halfmax = amplitude/2+spike[thpoint]
    
    comp = abs(smoothed-halfmax)
    
    hfpoint1 = argmin(comp[0:ind])
    hfpoint2 = ind+argmin(comp[ind:])
    
    width = dt*(hfpoint2-hfpoint1)
    err1 = (times[hfpoint1+wind/10+1]-times[hfpoint1-wind/10-1])*errsS[hfpoint1]/(smoothed[hfpoint1+wind/10+1]-smoothed[hfpoint1-wind/10-1]) 
    err2 = (times[hfpoint2+wind/10+1]-times[hfpoint2-wind/10-1])*errsS[hfpoint2]/(smoothed[hfpoint2+wind/10+1]-smoothed[hfpoint2-wind/10-1]) 
    errsL.append(sqrt(err1**2+err2**2))
    
    ahppoint = argmin(smoothed[ind:])+ind
    
    ahpsize = smoothed[thpoint]-smoothed[ahppoint]
    errsL.append(sqrt(errsS[thpoint]**2+errsS[ahppoint]**2))
    
    fig = figure()
    ax = fig.add_subplot(111)
    
    ax.plot(times,spike,'b-')
    ax.plot(times,smoothed,'g-',linewidth=2)
    
    ax.plot([times[thpoint]],smoothed[thpoint],'ro')
    ax.plot(times,smoothed[thpoint]*ones([len(times,)]),'y-')
    
    ax.plot(times[hfpoint1:hfpoint2],halfmax*ones([hfpoint2-hfpoint1,]),'k-')
    ax.plot(times[argmax(spike)]*ones([500,]),linspace(smoothed[thpoint],max(spike),500),'m-')
    ax.plot(times[ahppoint]*ones([100,]),linspace(smoothed[ahppoint],smoothed[thpoint],100),'c-')
    
    return smoothed[thpoint],amplitude,width,ahpsize,errsL



def simpleAvg(spkdic):  ### this function calculates the average spike
    
    ### input is the spikesdic from cleanPlot
    
    
    
    endtime = 100 # 0.26 for WT, 0.6 for SHK. this excludes spikes longer than a certain time, helps to eliminate some odd mistakes
    dts = []
    
    for j in spkdic.keys():        
        ts,xs = spkdic[j]
        if len(ts)>10 and (ts[-1]-ts[0])<endtime:
            dts.append(ts[1]-ts[0])
        else:
            del spkdic[j]
    dt = max(dts)
    print dt
    num = len(dts)
    for j in spkdic.keys():
        ts,xs = spkdic[j]
        ts = array(ts)-ts[0]
        decfac = int(round(dt/ts[1]))
        newlen = int(round(len(xs)*ts[1]/dt))
        print decfac,newlen
        decspike = []
        dects = []
        for k in range(newlen):
            decspike.append(xs[decfac*k])
            dects.append(k*dt)
        print dects[-1]
        spkdic[j] = dects,decspike
    
    lns = []
        
    for j in spkdic.keys():
        ts,xs = spkdic[j]
        lns.append(len(ts))
    ln = max(lns)
    for j in spkdic.keys():
        ts,xs = spkdic[j]
        if len(ts)<ln:
            diff = ln-len(ts)
            ts.extend(list(linspace(ts[-1]+dt,ts[-1]+dt*diff,diff)))
            xs.extend(xs[-1]*ones([diff,]))
            
        
    spike = ones([ln,])-1
    
    for j in spkdic.keys():
        ts,xs = spkdic[j]
        spike = spike+array(xs)/num
    
    print num
    
    stanD = ones([ln,])-1
    
    for j in spkdic.keys():
        ts,xs = spkdic[j]
        stanD = stanD + sqrt((spike-xs)**2/num**2)
        
    plotAllSpikes(spkdic)
    
    ### output is the average spike, the time points, and the standard error
    
    return spike,linspace(0,dt*ln,ln),stanD
   
def freqCalc(spiketimes):

    dts = array(spiketimes[1:])-array(spiketimes[0:-1])
    end = len(spiketimes)-1
    for i in range(len(dts)):
        if dts[i]>2.0:
            end = i
            break
    frq = (spiketimes[end]-spiketimes[0])/end
    
    return frq
            
            
            
def readSpikeTraces(filename):
    
    otd = tabdata(filename,500)
    notd = {}
    for k in otd.keys():
        if k and len(otd[k])>1:
            notd[k-1] = otd[0],otd[k]
    
    return notd        
    