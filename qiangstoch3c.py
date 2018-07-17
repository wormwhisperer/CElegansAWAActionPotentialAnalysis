from scipy import *
from scipy.integrate import *
from scipy.optimize import *
from scipy import signal
import time
from numpy.random import seed,poisson

from pylab import *
ion()

def minf(V,Va,Vb):

    return 0.5*(1 + tanh((V-Va)/Vb))

def tau(V,Va,Vb):
    
    return 1.0/cosh((V-Va)/(2*Vb))

def CaV(V,x):
    
    ### calcium channel fit to voltage clamp with potassium blockers and 0K+ buffer
    ### current has fast activation (Vm1/Vm2) and slower partial inactivation (Vm3/Vm4)
    ### Inactivation time scale also fit from data (T2/23/24)
    
    Vm1 = -21.6
    Vm2 = 9.17
    
    Vm3 = 16.2
    Vm4 = -16.1
    

    T1 = 1.0
    #T2 = 1.0/50.0
    T2 = 80*tau(V,23,24)

    
    mx = 3.612/4
    
    m1 = minf(V,Vm1,Vm2)**1
    m2 = minf(V,Vm3,Vm4)**1
    
    
    
    
    y = ones([3,])-1
    y[0] = (m1*m2/mx-x[0])/T1-m1*m2*x[1]/(mx*T1)-x[0]/(2*T2)+x[1]/(2*T2)
    y[1] = (x[0]-x[1])/(2*T2)
    
    return y

def CaV2(V,x):
    
    ### copy of calcium channel function without quick self-inactivaton. for testing.
    
    Vm1 = -19.6
    Vm2 = 9.17
    
    Vm3 = 19.2
    Vm4 = -16.1
    

    T1 = 1.0
    #T2 = 1.0/50.0
    T2 = 80*tau(V,23,24)
    
    mx = 3.612
    
    m1 = minf(V,Vm1,Vm2)**1
    m2 = minf(V,Vm3,Vm4)**1

    y = (m1*m2/mx-x[0])/T1
    
    return y

def SHK1(V,x):
    
    ### SHK-1 function, fit to subtractino of shk-1 current from wt at 0 Ca
    
    T = 30.0
    m = minf(V,2,10)
    
    return (m-x)/T

def BK1(V,x):
    
    ### inactivation of slow potassium current. fit to get delay to bursting and baseline
    
    T = 1200.0*1.5 # 1200
    m = minf(V,-42,5) #m = minf(V,-40,5)
    
    return (m-x)/T

def SLO12(V,x):
    
    #T = 2350.0 ## rough estimate from data
    TKL = 18000*1.5
    TKH = 2000*1.5
    vtk1 = -52. #-50
    vtk2 = 20.
    #T = 2500.0 ## rough estimate from data
    T = TKL+(TKH-TKL)*0.5*(1+tanh(V-vtk1)/vtk2)
    m = minf(V,-32.,2.) #minf(V,-30,4)
    
    return (m-x)/T

def SLO3(V,x):
    
    m = minf(V,13,20)
    T = 1000 # 1200
    
    return (m-x)/T

def SLO4(V,x):
    
    m = minf(V,-25,5)
    T = 1000.0 # 600
    
    return (m-x)/T
    
def full_model(x,t):
    
    cap = 1.5e-3 ### from qiang's measurement
    gCa = 0.1 #0.3 #0.3 ### fit to calcium channel voltage clamp
    gK = 1.5 #1.5 #0.5 ### fit to SHK-1 voltage clamp
    #gK = 0.0
    gK2 = 0.8 #0.5 #0.6 #16.5 #33.0 ### tuned to give resting potential around -75 mV
    #gK2 = 0.0
    gK3 = 0.3*1.3#1.1 #0.95 #0.7 #0.1 #10.0 ### tuned to make voltage jump after current step the right height
    gK4 = 1.0 #1.0 #0.6 #10.0 #6.4 ### fit from combined slo-1 and slo-2 currents
    gK5 = 0.7*1.3 #1.4 #1.2 #0.3 #0.3?
    gK6 = 0.0 #0.9 #0.3?
    gK7 = 0.1
    gL = 0.25 #0.4 ### fit to leak slope from calcium channel voltage clamp
    #gL = 1.0
    vL = -65.0 ### fit to leak slope from calcium channel voltage clamp
    vCa = 120.0 ### canonical
    vK = -84.0 ### canonical
    gKI = 5.0 ### fit to get spiking baseline
    fac = 0.4  ## inactivation of calcium
    

    V = x[0]
    Ca1 = x[1]
    Ca2 = x[2]
    Ca3 = x[3]
    SHK = x[4]
    BK = x[5]
    SLO = x[6]
    KB = x[7] #KB
    Ca = x[8]
    SLO2 = x[9]
    
    

    
    #gCa = It
    #KB = 0
    
    
    #ICa = gCa*Ca1*(1-Ca2)*(V-vCa) #### Calcium current with rapid self-inactivation
    ICa = gCa*(Ca1+fac*Ca2)*(V-vCa)
    IL = gL*(V-vL) #### Linear leak current
    #Kir = min((V-vK),0.05*(V-vK)+9.5) #### Inward-rectifying voltage-insensitive potassium current
    
    Kir = -log(1+exp(-0.2*(V-vK-gKI)))/0.2+gKI
    
    z = Ca/(Ca+10.0)  ### calicum-sensitive potassium current
    
    
    I = x[-1]*sin(2*pi*t/x[10])
    #I = x[-1]*(1+signal.square(2*pi*(t-1000.)/10000.,0.5))/2.0
    #I = x[-1]
    #I = 0.5*(sign(t-25000+x[10])+1)*x[-1]*(t-25000+x[10])/x[10]
    
    y = ones([12,])-1.0
    y[0] = (I-ICa-(gK*SHK+gK3*minf(V,-42,5)*(1-BK)+gK4*SLO+gK5*KB+gK6+gK7*SLO2)*(V-vK)-IL-gK2*Kir)/cap
    yCa = CaV(V,[Ca1,Ca2,Ca3])
    y[1] = yCa[0]
    y[2] = yCa[1]
    y[3] = yCa[2]
    y[4] = SHK1(V,SHK)
    y[5] = BK1(V,BK)
    y[6] = SLO3(V,SLO)
    y[7] = SLO12(V,KB)
    y[8] = (-1.25*ICa-194.75*x[8]/5)/10000.0  ### calcium ion concentration in nM, setting Caeq to 5nm
    y[9] = SLO4(V,SLO2)
           
    
    return y

def run_model(x):
    
    seed()
    
    start=time.time()
    
    I = 12.0+2.0*x[0]
    #gK = 50+x[1]
    inits= [-75.03,1.13e-5,1.13e-5,1.0,2.04e-7,0.0,3.48e-5,0.0,0.312,1.24e-9,x[2],I]
    nCa = 1000 #x[1]
    nBK = x[1]
    nKB = x[1]
    #chs = ones(n)
    tf = 100000
    ts = [0]
    t = 0
    tau = 1.0
    vals = [inits]
    
    
    
    while t<tf:
        rates = full_model(inits,0)
        dt = min([0.05/max(abs(rates[1:])),0.25])
        Vs = odeint(full_model,inits,linspace(t,t+dt,2))
        t = t+dt
        #It = I*(1+signal.square(2*pi*(t-1000)/10000,0.6))/2.0
        Ca1 = inits[1]
        Ca2 = inits[2]
        BK = inits[5]
        KB = inits[7]
        V = inits[0]
        #rns = rand(n)
        #rts = abs(tau*(chs-mK))*dt
        #chk = rts>rns
        #chs = abs(chs-chk)
        BKrate = dt*nBK*tau*BK1(V,BK)
        Carate = dt*nCa*tau*CaV(V,[Ca1,Ca2,0.5])
        Ca1flip = poisson(abs(Carate[0]))*sign(Carate[0])
        Ca2flip = poisson(abs(Carate[1]))*sign(Carate[1])
        KBrate = dt*nKB*tau*SLO12(V,KB)
        BKflip = poisson(abs(BKrate))*sign(BKrate)
        KBflip = poisson(abs(KBrate))*sign(KBrate)
        inits = Vs[-1,:]
        inits[1] = Ca1flip/nCa+Ca1
        inits[2] = Ca2flip/nCa+Ca2
        #inits[3] = sum(chs)/n
        inits[5] = BKflip/nBK+BK
        inits[7] = KBflip/nKB+KB
        #inits[-1] = It
        ts.append(t)
        vals.append(inits)
        
    vals = array(vals)
    ts2,xs2 = downsample(ts,vals[:,0],0.25)
    vals2 = ones([len(ts2),len(inits)])-1
    vals2[:,0] = xs2
    
    for i in range(1,len(vals[0,:])):
        ts2,xs2 = downsample(ts,vals[:,i],0.25)
        vals2[:,i] = array(xs2)
        
    end = time.time()
    print end-start
    
    return ts2,vals2

def outsplot(outd,ind,ind2,ex):
    
    fig = figure()
    axd = {}
    ln = len(outd[outd.keys()[0]])+0.0
    n1 = int(sqrt(ln))
    n2 = int(ceil(ln/n1))
    
    for i in range(n2):
        for j in range(n1):
            if n1*i+j<ln:
                axd[n1*i+j+1] = fig.add_subplot(n2,n1,n1*i+j+1)
                for k in range(5):
                    ts,vls = outd[k+ex,ind][n1*i+j]
                    axd[n1*i+j+1].plot(ts,vls[:,0])
            
    return fig,axd

def downsample(ts,xs,step):
    
    
    
    #start = time.time()
    dts = array(ts[1:])-array(ts[0:-1])
    maxL = int(step/min(dts))+1
    first = 0
    num = int(ts[-1]/step)
    ts2 = linspace(0,ts[-1],num)
    ts3 = []
    xs2 = []
    
    for t in ts2:
        ind = findInd(ts[first:first+maxL],t)
        first = first+ind
        ts3.append(ts[first])
        xs2.append(xs[first])
        # if ts3[-1]>7000:
        #     return ts3,xs2
    #end = time.time()
    #print end-start
    
    return ts3,xs2
    
    
    
def findInd(vec,val):
    #start = time.time()
    c1 = array(vec)-val
    #end = time.time()
    #print end-start
    return argmin(abs(c1))

def writespikesDic(dic,filename):
    
    ln = len(dic.keys())
    lns = []
    for k in dic.keys():
        ts,vs = dic[k]
        lns.append(len(ts))
    ln2 = max(lns)
    f = file(filename,'wb')
    lst = sorted(dic.keys())
    for j in range(ln2):
        if 1000*int(j/1000)==j:
            print j
        for k in range(ln-1):
            ts,vs = dic[lst[k]]
            if j<len(ts):
                f.write(str(ts[j])+'\t'+str(vs[j])+'\t')
            else:
                f.write(' '+'\t'+' '+'\t')
        ts,vs = dic[lst[-1]]
        if j<len(ts):
            f.write(str(ts[j])+'\t'+str(vs[j])+'\n\r')
        else:
            f.write(' '+'\t'+' '+'\n\r')
    f.close()
    
    return

from multiprocessing import Pool
# 
# def testfunc(x):
#     print x[0],x[1]
#     return 1
# 
#outd = {}

#start = time.time()
#fcs = [1.1,1.2,1.3,1.4]
#nss = [100,200,500,2000,5000]
#fnss = [[1.1,500],[1.1,2000],[1.2,200],[1.2,500],[1.2,2000],[1.3,200],[1.3,5000]]

#workL = []

# for i in range(5):
#     for j in range(4):
#         for k in range(5):
#             print '####### working on i = ',i,'k=',k,' #########'
#             b = list(fcs[j]*ones([14,]))
#             c = list(nss[k]*ones([14,]))
#             aL = [-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6]
#             workL.extend(zip(aL,b,c))
        
#p = Pool(20)
#outL = p.map(run_model,workL)
#p.close()
            
#end = time.time()

#print '###### full runtime',(end-start)/3600,'hours ######'

# p = Pool(20)
# outs = p.map(run_model,zip(list(6.0*ones([60,])),list(15.0*ones([60,]))))
# p.close()
