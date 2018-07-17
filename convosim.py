
from pylab import *
from scipy import *
from scipy.integrate import *
from plotting import *
from scipy.ndimage import filters
from scipy.signal import gaussian


### How to use: read in folder of data using function folderread (check that correct column is being read)
### then push dictionary through filtering function (filterAll). check results with plotData. run filtered data through cleanPlot2
### be sure to adjust window sizes, detection thresholds, and spike width settings


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

def folderread(directory):
    
    outd = {}
    
    lst = os.listdir(directory)
    j = 0
    for i in lst:
        if (i.find('.txt')>0) & (len(i)<100):
            dic = tabdata(directory+'/'+i,5)
            a = i.strip('.txt')
            outd[j] = dic[2]
            j += 1
            
    return outd

def runningInt(dat,tim):
    
    DelT = 0.4
    dt = tim[1]-tim[0]
    init = int(1.0/dt)
    fin = int(7.0/dt)
    intv = int(DelT/dt)
    diffs = [0]
    
    tms = tim[init:fin]
    
    for i in range(init+1,fin):
        if i-init<intv:
            sl = (dat[i]-dat[init])/i
            corr = dat[init:i] - sl*linspace(init,i,i-init)-dat[init]
            diffs.append(sum(corr))
            #diffs.append(sum(dat[0:i]))
        else:
            sl = (dat[i]-dat[i-intv])
            corr = dat[i-intv:i] - sl*linspace(0,intv,intv)-dat[i-intv]
            diffs.append(sum(corr))
            #diffs.append(sum(dat[i-intv:i]))
    
    return diffs,tms

def runningDeriv(dat,tim):
    
    DelT = 1.0
    dt = tim[1]-tim[0]
    init = int(1.0/dt)
    fin = int(7.0/dt)
    intv = int(DelT/dt)
    diffs = []
    
    for i in range(init+intv,fin):
        a,b,c = polyfit(linspace(0,intv,intv),dat[i-intv:i],2)
        diffs.append(a)
        
    tms = tim[init+intv/2:fin-intv/2]
    
    return diffs,tms

def diffthresh(ts,xs):

    dt = ts[1]-ts[0] ## time bin
    nframes = int(0.18/dt) #NOW: 0.35 for wt, 0.48 SHK #0.18 for wt, 0.25 for shk, 0.04 for wt sim, 0.08 for shk sim
    winspace = 1 #int(0.06/dt) # spacing to check for derivatives above threshhold
    initpoint = int(1.3/dt) # start time in seconds
    endpoint = int(5.9/dt) # end time
    
    ns = (endpoint-initpoint)
    
    dthresh = 0.009 #NOW: 0.009 for WT, 0.012 SHK #0.013 for wt, 0.017 for shk, 0.006 for wt sim, 0.008 for shk sim
    dthresh2 = 0.007  #NOW: 0.008 for WT, 0.01 SHK # 0.013 for wt, 0.017 for shk, 0.0035 for wt sim, 0.005 for shk sim
    vthresh = -100  # overall voltage threshhold, don't use
    upcheck = list(ones([winspace,])-1)  # index vector for upswings
    downcheck = list(ones([winspace,])-1) # index vector for downswings
    indicator = 0
    
    point = initpoint
    
    while point<endpoint:
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
        

    return upcheck,downcheck,ts[initpoint:endpoint],xs[initpoint:endpoint]

def filterAll(dic,ts):
    
    flt = exp(-ts/1.0)
    
    otd = {}
    for i in dic.keys():
        vec = decon(smooth(norml(dic[i]),4),flt)
        #otd[i] = smooth(vec,4)
        otd[i] = smooth(vec,8)
    
    return otd

def analysePlot(dic,ts2):
    
    num = len(dic.keys())
    
    fig = figure()
    axd = {}
    n1 = int(sqrt(num))
    n2 = int(ceil((num+0.0)/n1))
    outd = {}
    n = 0
    

    for i in range(n2):
        for j in range(n1):
            axd[n1*i+j+1] = fig.add_subplot(n2,n1,n1*i+j+1)
            ts3,vs3 = dic[n1*i+j]
            #upc,downc,ts,xs = diffthresh(ts2,dic[n1*i+j])
            upc,downc,ts,xs = diffthresh(ts3,vs3)
            outd[n] = dic[n1*i+j],upc,downc,ts,xs
            n+=1
            axd[n1*i+j+1].plot(ts,xs)
            ups = nonzero(upc)[0]
            downs = nonzero(downc)[0]
            tsu = [ts[k] for k in ups]
            xsu = [xs[k] for k in ups]
            tsd = [ts[k] for k in downs]
            xsd = [xs[k] for k in downs]
            axd[n1*i+j+1].plot(tsu,xsu,'g.')
            axd[n1*i+j+1].plot(tsd,xsd,'r.')
    
    return fig,axd,outd

def jumps(vec):  ## function for finding the indices of up and downswings
    
    dupc = array(vec[1:])-array(vec[0:-1])
    if vec[0]==1:
        upinds = [0]
        upinds.extend(nonzero(abs(dupc))[0])
    else:
        upinds = nonzero(abs(dupc))[0]
    
    return upinds
    


def cleanTimes2(upc,downc,ts,xs): #### this function cleans up the initial guess of spikes
    
    ### input dic is the traces dictionary returned by readData
    
    upinds = jumps(upc)
    downinds = jumps(downc)
    upi = []
    di = []
    spiketimes = []
    dt = ts[1]-ts[0]
    spikewidth = int(round(0.5/dt))  ## 0.5 for shk, 0.1 for wt, approx width of a spike 
    spikedic = {}
    
    upDel = []
    dwnDel = []
    for i in range(len(upinds)/2):
        if (upinds[2*i+1]-upinds[2*i])<3:
            upDel.append(2*i)
            upDel.append(2*i+1)
    for i in range(len(downinds)/2):
        if (downinds[2*i+1]-downinds[2*i])<2:
            dwnDel.append(2*i)
            dwnDel.append(2*i+1)
    upinds = delete(upinds,upDel)
    downinds = delete(downinds,dwnDel)
    
    
    
    
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
                if len(upinds)*len(downinds)==0:
                    print 'breaking'
                    break
            else:
                check1 = 0
        if len(upinds)*len(downinds)==0:
                    print 'breaking'
                    break
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
        print len(upinds),len(downinds)
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
            if max(xs[upinds[0]:downinds[1]])>-0.025:
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
        timeafter = 0.5 ## time after point of max slope
        ### this code finds the point of max slope
        wind = int(round(2e-2/dt))
        c = real(gaussian(wind,wind/2))
        smoothed = filters.convolve1d(real(xs[upi[i]:di[i]]),c/c.sum())
        ds = smoothed[1:]-smoothed[0:-1]
        ind = argmax(ds)+upi[i]
        
        spikedic[j] = ts[ind-int(round(timebefore/dt)):ind+int(round(timeafter/dt))],xs[ind-int(round(timebefore/dt)):ind+int(round(timeafter/dt))]
        j+=1
    
    spiketimes = 0.5*dt*(array(upi)+array(di))
    
    ### outputs are indices of upswings, downswings, and the time and potential vectors for each spike
    ### the list of spike times, and a dictionary with the traces of each spike
    
    return upi,di,ts2,xs2,spiketimes,spikedic

def cleanPlot2(dic,ts3):  ### this function plots the final identified spikes
    
    ### input dic is the traces dictionary returned by readData
    
    num = len(dic.keys())
    
    fig = figure()
    axd = {}
    n1 = int(sqrt(num))
    n2 = int(ceil((num+0.0)/n1))
    
    outd = {}

    
    
    for i in range(n2):
        for j in range(n1):
            axd[n1*i+j+1] = fig.add_subplot(n2,n1,n1*i+j+1)
    for k in range(num):
        ts4,vs4 = dic[dic.keys()[k]]
        #upc,downc,ts,xs = diffthresh(ts3,dic[dic.keys()[k]])
        upc,downc,ts,xs = diffthresh(array(ts4)/1000,vs4[:,0]/1000)
        upi,di,ts2,xs2,spiketimes,spikedic = cleanTimes2(upc,downc,ts,xs)
        outd[k] = dic[dic.keys()[k]],upi,di,ts2,xs2,spiketimes,spikedic
        axd[k+1].plot(ts,xs,'-')
        axd[k+1].plot(ts2,real(xs2),'r.')
        axd[k+1].set_xticklabels([])
        axd[k+1].set_yticklabels([])
        print k
    

        
    ### output is the figure object, a dictionary of axes, the dictionary of individual spikes
    ### and a list of spike times
    
    return fig,axd,outd

def plotData(dic):
    
    num = len(dic.keys())
    #ts = linspace(0,len(dic[0])/100,len(dic[0]))
    
    fig = figure()
    axd = {}
    n1 = int(sqrt(num))
    n2 = int(ceil((num+0.0)/n1))
    outd = {}
    n = 0
    

    for i in range(n2):
        for j in range(n1):
            axd[n1*i+j+1] = fig.add_subplot(n2,n1,n1*i+j+1)
            #axd[n1*i+j+1].plot(norml(dic[dic.keys()[n1*i+j]]))
            ts,vs = dic[n1*i+j]
            axd[n1*i+j+1].plot(ts,vs)
    
    return fig,axd

def plotresults(ts,o1,o2,o3):
    
    upc,downc,ts1,xs = diffthresh(ts,smooth(o1,2))
    upc2,downc2,ts2,xs2 = diffthresh(ts,smooth(o2+0.2,2))
    upc3,downc3,ts3,xs3 = diffthresh(ts,smooth(o3+0.4,2))

    ups = nonzero(upc)[0]
    downs = nonzero(downc)[0]
    ups2 = nonzero(upc2)[0]
    downs2 = nonzero(downc2)[0]
    ups3 = nonzero(upc3)[0]
    downs3 = nonzero(downc3)[0]
    
    tsu = [ts1[i] for i in ups]
    xsu = [xs[i] for i in ups]
    tsd = [ts1[i] for i in downs]
    xsd = [xs[i] for i in downs]
    tsu2 = [ts2[i] for i in ups2]
    xsu2 = [xs2[i] for i in ups2]
    tsd2 = [ts2[i] for i in downs2]
    xsd2 = [xs2[i] for i in downs2]
    tsu3 = [ts3[i] for i in ups3]
    xsu3 = [xs3[i] for i in ups3]
    tsd3 = [ts3[i] for i in downs3]
    xsd3 = [xs3[i] for i in downs3]
    
    fig2 = figure()
    ax2 = fig2.add_subplot(111)
    
    ax2.plot(ts1,xs,'b-',linewidth=2)
    ax2.plot(ts2,xs2,'b-',linewidth=2)
    ax2.plot(ts3,xs3,'b-',linewidth=2)
    
    ax2.plot(tsu,xsu,'go')
    ax2.plot(tsu2,xsu2,'go')
    ax2.plot(tsu3,xsu3,'go')
    
    ax2.plot(tsd,xsd,'ro')
    ax2.plot(tsd2,xsd2,'ro')
    ax2.plot(tsd3,xsd3,'ro')
    
    return fig2,ax2

