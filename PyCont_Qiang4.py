from scipy import *
from pylab import *
from PyDSTool import *
import time
import os
import shutil
from PyDSTool.Toolbox import phaseplane as pp
ion()

fc = 1.1
fc2 = 1.5

pars = {'Iapp': 15.0,
        'C': 1.5e-3,
        'vK': -84.,
        'gK': 1.5, #1.5, 
        'gK2': 0.8, #0.8 (wt)
        'gK3': 0.3*fc, #0.5 (shk1) 0.3 (wt) or 0.75
        'gK4': 1.0, #1.0 (wt), 2.0 (shk)
        'gK5': 0.7*fc, #2.2 (shk1) 0.7 (wt) or 2.9
        'gK6': 0.0,
        'gK7': 0.1,
        'vCa': 120.,
        'gCa': 0.1, #0.13 (shk), 0.1 (wt) or 0.16
        'vL': -65.,
        'gL':0.25, 
        'vm1': -21.6, #-21.6,
        'vm2': 9.17,
        'vm3': 16.2, #16.2,
        'vm4': -16.1,
        'TK': 30.0, #30
        'TC1': 1.0,
        'TC2': 80.0,
        'Tbk': 1200*fc2, #1200,
        'Tkb': 6000, #2500,
        'vk1': 2.0, #2.0,
        'vk2': 10.0,
        'vs1': 13.0, #13.0,
        'vs2': 20.0,
        'vq1': -25.0, #-25.0,
        'vq2': 5.0,
        'TS': 1000.0, # 600
        'TS2': 1000.0, # 600
        'vt1': 20.0, #20.0,
        'vt2': 24.0,
        'vb1': -42, #-40.,
        'vb2': 5.0,
        'vp1': -32, # -30,
        'vp2': 2.0, # 2
        'gKI': 5.0, #5.0,
        'mx': 3.612/4,
        'fac': 0.4,
        'TI': 2000,
        'TKL': 18000.*fc2, #18000
        'TKH': 2000.*fc2, #2000
        'vtk1': -52, #-50
        'vtk2': 20, #20
        'T1': 25000} 

#icdict = {'v': -22, 'w': 0.00501, 'c1': 0.441, 'c2': 0.441, 'bk': 0.9999, 'slo': 0.1}
icdict = {'v': -74.99, 'w': 2.05e-7, 'c1': 9.71e-6, 'c2': 9.71e-6, 'slo2': 2.069e-9, 'slo': 0.000151, 'bk': 8.35e-7, 'kb': 1.7e-10}

# Set up model
auxfndict = {'minf': (['v'], '0.5*(1 + tanh((v-vm1)/vm2))'), \
             'qinf': (['v'], '0.5*(1 + tanh((v-vq1)/vq2))'), \
			 'winf': (['v'], '0.5*(1 + tanh((v-vm3)/vm4))'), \
             'xinf': (['v'], '0.5*(1 + tanh((v-vk1)/vk2))'), \
             'yinf': (['v'], '0.5*(1 + tanh((v-vb1)/vb2))'), \
             'zinf': (['v'], '0.5*(1 + tanh((v-vs1)/vs2))'), \
			 'tau': (['v'], '1.0/cosh((v-vt1)/(2*vt2))'), \
             'pinf': (['v'], '0.5*(1+tanh((v-vp1)/vp2))'), \
             'ic': (['t'], 'heav(6000-t)*heav(t-1000)'), \
             'ic2': (['t'], 'heav(t-25000+T1)*(t-25000+T1)/T1'), \
             'gkt': (['v'], 'TKL+(TKH-TKL)*0.5*(1+tanh(v-vtk1)/vtk2)'), \
             'kir': (['v'], '-log(1+exp(-0.2*(v-vK-gKI)))/0.2+gKI') \
			}

#vstr = '(Iapp - gCa*c1*(1-c2)*(v-vCa) - (gK*w+gK3*yinf(v)*(1-bk)+gK4*slo+gK5*kb)*(v-vK) - gK2*kir(v) - gL*(v-vL))/C'
vstr = '(Iapp*ic(t) - gCa*(c1+fac*c2)*(v-vCa) - (gK*w+gK7*slo2+gK4*slo+gK6+gK3*yinf(v)*(1-bk)+gK5*kb)*(v-vK) - gK2*kir(v) - gL*(v-vL))/C'
wstr = '(xinf(v)-w)/TK'
c1str = '(minf(v)*winf(v)/mx-c1)/TC1'
c1str = '(minf(v)*winf(v)/mx-c1)/TC1-minf(v)*winf(v)*c2/(mx*TC1)-c1/(2*TC2*tau(v))+c2/(2*TC2*tau(v))'
c2str = '(c1-c2)/(2*TC2*tau(v))'
bkstr = '(yinf(v)-bk)/Tbk'
slostr = '(zinf(v)-slo)/TS'
slo2str = '(qinf(v)-slo2)/TS2'
kbstr = '(pinf(v)-kb)/gkt(v)'

DSargs = args(name='Phil6')
DSargs.pars = pars
#DSargs.varspecs = {'v': vstr, 'w': wstr, 'c1': c1str, 'c2': c2str, 'bk': bkstr, 'slo': slostr}
DSargs.varspecs = {'v': vstr, 'w': wstr, 'c1': c1str, 'c2': c2str, 'slo2': slo2str, 'slo': slostr, 'bk': bkstr, 'kb': kbstr}
DSargs.fnspecs = auxfndict
DSargs.ics = icdict
DSargs.pdomain = {'gK6': [0.0,100.0],'gK': [0.0,100.0], 'gK2': [0.0,1000.0]}

#DSargs.xdomain = {'v': [-100, 50], 'w': [0, 1], 'c1': [0, 1], 'c2': [0, 1], 'slo2': [0, 1], 'slo': [0, 1]}
DSargs.tdata = [0,7000]
DSargs.algparams = {'init_step': 0.1,'max_pts': 10000000}

testDS = Generator.Vode_ODEsystem(DSargs)

#shutil.move('dop853_temp/_dop853_Qiang_vf.so',os.getcwd())

def runOnce(I,T1,ax):

    start = time.time()
    testDS.pars['Iapp'] = I
    testDS.pars['T1'] = T1
    traj_on = testDS.compute('force')
    pts_on = traj_on.sample()
    end = time.time()
    print end-start
    ax.plot(pts_on['t'],pts_on['v'])
    #ax2.plot(pts_on['t'],pars['gK3']*0.5*(1+tanh((pts_on['v']-pars['vb1'])/pars['vb2']))*(1-pts_on['bk'])+pars['gK5']*pts_on['kb'])
    
    
    
    
    return pts_on,ax



#ln = len(pts_on['w'])
#plot(pts_on['t'][ln/2:], pts_on['v'][ln/2:])

#pp.plot_PP_vf(testDS,'v','w',scale_exp=-1)

#fp_coord = pp.find_fixedpoints(testDS, n=4, eps=1e-8)[0]
#fp = pp.fixedpoint_2D(testDS, Point(fp_coord), eps=1e-8)

# n=3 uses three starting points in the domain to find nullcline parts, to an
# accuracy of eps=1e-8, and a maximum step for the solver of 0.1 units.
# The fixed point found is also provided to help locate the nullclines.
#nulls_x, nulls_y = pp.find_nullclines(testDS, 'v', 'w', n=3, eps=1e-8, max_step=0.1,
                             #fps=[fp_coord])

# plot the fixed point
#pp.plot_PP_fps(fp)

# plot the nullclines
#plt.plot(nulls_x[:,0], nulls_x[:,1], 'b')
#plt.plot(nulls_y[:,0], nulls_y[:,1], 'g')
# 

# PCargs.name = 'LC1'
# PCargs.type = 'LC-C'
# PCargs.initpoint = 'EQ1:H1'
# PCargs.MinStepSize = 0.005
# PCargs.MaxStepSize = 1.0
# PCargs.StepSize = 0.01
# PCargs.MaxNumPoints = 5000
# PCargs.NumSPOut = 40;
# PCargs.LocBifPoints = 'LPC'
# PCargs.SolutionMeasures = ['avg','min','max','nm2']
# PCargs.SaveEigen = True
# PyCont.newCurve(PCargs)
# 
# print('Computing curve...')
# start = clock()
# PyCont['LC1'].forward()
# print('done in %.3f seconds!' % (clock()-start))


#PyCont['LC1'].display(('Iapp','v_min'),stability=True)
#PyCont['LC1'].display(('Iapp','v_max'),stability=True)


#plt.xlim([0, 300])
#plt.ylim([-75, 75])
#PyCont.plot.fig1.axes1.axes.set_title('Bifurcation Diagram')

#PyCont['LC1'].plot_cycles(figure='fig2', method='stack', exclude='P2', tlim='5T')

#PyCont['EQ1'].display(('v','w'),stability=True, figure='new')

#PyCont['LC1'].plot_cycles(coords=('v','w'), figure='fig3', exclude='P2')

#os.remove('_dop853_Qiang_vf.so')

def runtest(ics,gcs,ax):
    
    start = time.time()
    amps = []
    #prs = []
    for j in range(len(ics)):
        testDS.pars['Iapp'] = ics[j]
        testDS.pars['gK6'] = gcs[j]
        
        #testDS.update()
    
        traj_on = testDS.compute('force')
        pts_on = traj_on.sample()
        amp = findAmp(pts_on)
        #perd = findPeriod(pts_on['t'],pts_on['v'])
        amps.append(amp)
        #prs.append(perd)
    mxa = max(amps)
    #mxp = max(prs)
    for j in range(len(ics)):
        if not ((gcs[j]>0.13*(ics[j]-14)/6.5) and (amps[j]<5.0)):
            ax.plot([ics[j]],[gcs[j]],'.',color=(amps[j]/mxa,0,max(1-amps[j]/mxa,0)))
            #ax2.plot([ics[j]],[gcs[j]],'.',color=(min(prs[j]/mxp,1),0,max(1-prs[j]/mxp,0)))
    end = time.time()
    print end-start
        
    return amps,ax

def findAmp(pts_on):
    
    ln = len(pts_on['v'])
    dts1 = (pts_on['t']-2000)/abs(pts_on['t']-2000)
    dts2 = (pts_on['t']-8950)/abs(pts_on['t']-8950)
    #print min(dts2),max(dts2)
    ind1 = list(dts1).index(1.0)
    ind2 = list(dts2).index(1.0)
    mxpt = argmax(pts_on['v'][ind1:ind2])+ind1
    maxv = pts_on['v'][mxpt]
    maxt = pts_on['t'][mxpt]
    #dts3 = (pts_on['t']-maxt-1000)/abs(pts_on['t']-maxt-1000)
    #ind3 = list(dts3).index(1.0)
    #if ind3>ind2:
        #ind3 = ind2
    minv = min(pts_on['v'][mxpt:ind2+1])
    amp = maxv-minv
    
    return amp

def gKtest(gK3,gK5):
    
    gKs = []
    ics = linspace(6,20,150)
    vb1 = -49.
    vb2 = 5.
    testDS.pars['gK3'] = gK3
    testDS.pars['gK5'] = gK5
    testDS.pars['gK6'] = 0.0
    
    
    for ic in ics:
        testDS.pars['Iapp'] = ic
        traj_on = testDS.compute('force')
        pts_on = traj_on.sample()
        gKs.append(min(gK3*0.5*(1+tanh((pts_on['v']-pars['vb1'])/pars['vb2']))*(1-pts_on['bk'])+gK5*pts_on['kb']))
    
    return gKs

def modeltest(axl,axl2):
    
    start = time.time()
    
    ics = linspace(0,30,16)
    mns = []
    
    for ic in ics:
        testDS.pars['Iapp'] = ic
        traj_on = testDS.compute('force')
        pts_on = traj_on.sample()
        axl.plot(pts_on['t'],pts_on['v'])
        axl2.plot(pts_on['t'],pars['gK3']*0.5*(1+tanh((pts_on['v']-pars['vb1'])/pars['vb2']))*(1-pts_on['bk'])+pars['gK5']*pts_on['kb'])
        #axl2.plot(pts_on['t'],pts_on['kb'])
        dts1 = (pts_on['t']-1100)/abs(pts_on['t']-1100)
        dts2 = (pts_on['t']-5950)/abs(pts_on['t']-5950)
        #print min(dts2),max(dts2)
        ind1 = list(dts1).index(1.0)
        ind2 = list(dts2).index(1.0)
        mns.append(min(pars['gK3']*0.5*(1+tanh((pts_on['v'][ind1:ind2]-pars['vb1'])/pars['vb2']))*(1-pts_on['bk'][ind1:ind2])+pars['gK5']*pts_on['kb'][ind1:ind2]))
    df = max(mns[5:])-min(mns[5:])
    print df
    
    end = time.time()
    print end-start,'seconds'
        
    return axl,axl2,mns

def doCont2():
    
    PyCont = ContClass(testDS)

    PCargs = args(name='EQ1', type='EP-C')
    PCargs.freepars = ['Iapp']
    PCargs.StepSize = 0.01
    PCargs.MaxNumPoints = 2000
    PCargs.MaxStepSize = 0.1
    PCargs.LocBifPoints = ['H','LP']
    PCargs.StopAtPoints = 'B'
    PCargs.verbosity = 2
    PCargs.SaveEigen = True
    PyCont.newCurve(PCargs)
    
    print('Computing curve...')
    start = clock()
    PyCont['EQ1'].forward()
    PyCont['EQ1'].backward()
    print('done in %.3f seconds!' % (clock()-start))
    
    #PyCont.display(('v','Iapp'),stability=True)

    
    PCargs = args(name='HO2', type='H-C1')
    PCargs.initpoint = 'EQ1:H1'
    PCargs.freepars = ['Iapp', 'gK6']
    PCargs.StopAtPoints = 'B' #['B', 'DH']
    PCargs.MaxNumPoints = 1200
    PCargs.MaxStepSize = 0.1
    PCargs.MinStepSize = 0.00000001
    PCargs.LocBifPoints = 'B' #'all'
    PCargs.verbosity = 2
    PyCont.newCurve(PCargs)
    
    print('Computing Hopf curve...')
    start = clock()
    PyCont['HO2'].backward()
    PyCont['HO2'].forward()
    print('done in %.3f seconds!' % (clock()-start))
    # 
    # 
    # 
    # 
    # PyCont['HO2'].display(('Iapp','gK6'), stability=True)
    
    icurve = PyCont['HO2'].sol['Iapp']
    gkcurve = PyCont['HO2'].sol['gK6']
    
    #sol = PyCont['EQ1'].sol
    
    return icurve,gkcurve

def localMax(vec):
    
    mn = mean(vec)
    dt = array(vec[1:])-array(vec[0:-1])
    dt2 = dt/abs(dt)
    ddt = list(-dt2[1:]+dt2[0:-1])
    indices = [i for i, x in enumerate(ddt) if (x == 2 and vec[i+2]>mn)]
    
    return array(indices)+2

def findPeriod(ts,vs):
    
    tm = ts[-1]
    tind = array(ts)-tm/3
    tchk = list(tind/abs(tind))
    start = tchk.index(1)
    ids = localMax(vs[start:])
    tms = []
    for ind in ids:
        tms.append(ts[ind+start])
    if len(tms)>3:
        a,b = polyfit(range(len(ids)),tms,1)
    else:
        a = 0
        
    return a
    
    