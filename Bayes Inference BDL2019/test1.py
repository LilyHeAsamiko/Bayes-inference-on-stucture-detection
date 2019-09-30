# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 10:58:49 2019

@author: LilyHeAsamiko
"""
#from obspy import read, read_inventory
#from obspy.io.xseed import Parser
#from obspy.signal import PPSD
#sd = read("https://examples.obspy.org/BW.KW1..EHZ.D.2011.038")
#datasd = sd.select(id="BW.KW1..EHZ")[0]
#
#inv = read_inventory("https://examples.obspy.org/BW_KW1.xml")
#print(datasd.stats) 
#ppsd = PPSD(datasd.stats, metadata=inv)
#ppsd.add(sd)
#ppsd.abs
import unicodedata
import math
test = math.inf
from obspy.signal.tf_misfit import plot_tfr
import seaborn as sns


def SPECPS(RL, xx, yy, rho, winL, T, prop):
#corner:5: fcp = 0.32*vs/rho  3 fcp = 1.53*vp/2pi*rho   
#    n = np.linspace(0, T, len(ref))
if prop = 'R':
    RL1 = RL
    if T % winL > 0:
        T = int(T/winL)*winL
#        ref = ref[0:T-1]
        RL = RL[0:int(T)]
        Re = 6371
#    WW = []
#    for dw in dW:
#    winL = dw/360*2*np.pi*Re
       
    freq = np.fft.fftfreq(T)[0:int(T/2)]
    freq = freq[freq >= 0]
    t0 = np.linspace(0, T-winL, int(len(RL)/winL))
    fmin = min(freq)
    fmax = max(freq)
    RLT = []
    v = []
    k = []
    Vg = []
    x = []
    a = []
    TRL = []
    Phii = []
    vavg = 0
    fi = 1/winL
    
    for i in range(len(t0)): 
        print(i)
        t = t0[i]
        pti = np.zeros((winL, 1))
        if i+1 < len(t0):
            ti = np.array(range(int(t), int(t0[i+1])))
        else:
            break
        p0 = RL[ti]
#        py0 = yy[ti]
        dt = (t0[i+1]-t)/winL
        Amax = max(p0)
        Amin = min(p0)
        A0 = p0[0]
        At = p0[-1]
        if abs(A0 - At)/abs(Amax - Amin):
#            A = max(abs(p0))
            h0 = 8
            b = 25            
            pf = np.fft.fft(p0)
            w0  = 80.897/180*np.pi
            phi0 = 15.754/180*np.pi
            Re = 6371
            dx = xx[winL]*360*2*np.pi*Re
            dy = yy[winL]*360*2*np.pi*Re
            delta = np.sqrt(dx**2+dy**2)
            wi= np.arctan(pf.real/pf.imag)
            Hi = (np.sign(ti-t)+1)/2            
#            pti = A*np.exp(-2*(ti-t))
            pti = h0*np.exp(-(delta/b)**4)   
#            plt.figure(),
#            plt.plot(pti)
            pti *= np.fft.fft(np.cos(2.*np.pi*fi*(ti-t) + wi * np.pi) * Hi)  
#            plt.figure(),
#            plt.plot(pti)
#        plot_tfr(pti, dt=dt, fmin=0, fmax=fi+winL*dt)
            pti += p0
        else:
            if len(t0) > i+2:
                ti2 = np.array(range(int(t), int(t0[i+2])))
            else:
                ti2 = np.array(range(int(t), int(2*t0[-1]-t0[-2])))
            p0 = xx[ti2]
            p0 = yy[ti2]
 #           A = max(abs(p0))            
            pf = np.fft.fft(p0)
            wi = np.arctan(pf.real/pf.imag)
            Hi = (np.sign(ti2-t)+1)/2
            pti = h0*np.exp(-(delta/b)**4) 
            pti *= np.fft.fft(np.cos(2.*np.pi*fi*(ti-t) + wi * np.pi) * Hi)   
            pti += p0
            pti = pti[0:winL]
            wi = wi[0:winL]
            #            x = pti
#            x = list(x)
        vi = pti/(ti-t+0.000001)
#            v = list(v)
#        phii = list(phii)
        ki = wi/ vi
        if abs(vi).all()< 1:
            x.append(0)
            v.append(0)
        else:
            print(ki)
            print(max(ki)-min(ki))
            if max(ki)-min(ki) < 2:
                Phii.append(wi)
                phii = wi
                v.append(vi)
                k.append(phii/vi)
                x.append(pti)
                PT.append(pti)
        #            k = np.array(k)
                if i > 0:
                    vg=(Phii[i]-Phii[i-1])/(k[i]-k[i-1])
                    ai = [vi[0]-v[i-1][-1],vi[1]-vi[0],vi[2]-vi[1],vi[3]-vi[2],vi[4]-vi[3]]
                else:
                    vg=Phii[i]/k[i]  
                    ai = []
                    ai.append(vi[0])
                    if len(vi) > 1:
                        for i in range(len(vi)-1):
                            ai.append(vi[i+1]-vi[i])
#                    D.append(x)
                Vg.append(vg)
                TP.append(ti)
                a.append(ai)
            else:
                PT1 = PT
                v1 = v
                k1 = k
                Vg1 = Vg
                x1 = x
                a1 = a
                TP1 = TP
                Phii1 = Phii
                vavg1 = vavg
                print([PT1,v1,k1,Vg1,x1,a1,TP1,Phii1,vavg1])
                pti, vi, ki, vgi, xi, ai, ti, phii, vavg = P(p0, rho, int(len(p0)/2), int(len(p0)))
                print([pti, vi, ki, vgi, xi, ai, ti, phii, vavg])
                
                Phii.append(phii)
                v.append(vi)
                k.append(ki)
                x.append(pti)
                PT.append(pti)
                Vg.append(vgi)
                TP.append(ti)
                a.append(ai)
        if i == len(t0) - 1:
            break
        
        r = np.sqrt(xy**2+xx**2)
        alpha = -np.log(y/r**0.5)/xz
        AL = r**(-0.5)*exp(-alpha*xz)
        X,Y = np.meshgrid(xx,xy)
       
        plt.figure(),
        plt.subplot(211)
        plt.pcolor(X,Y,AL)
        
        plt.figure(),
        plt.subplot(211)
        plt.pcolor(X,Y,AL)        
        plt.title('Love wave grad')
        #inversion
        accept, AR, Residule, dK, dTheta = MCMC_2D(x, y, Tx, alpha, iters)
        std_Residule = []
        std_dK = []
        std_dTheta =[]
        for i in range(len(AL)):
            std_Residule.append(3*np.std(Residuel[0:i]))
            std_dK.append(3*np.std(dK[0:i])) 
            std_dTheta.append(3*np.std(dTheta[0:i]))
        

        plt.subplot(212),
        plt.pcolor(X, Y, AR)
        plt.title('Rayleigh wave grad')
        
        SD1, SD2, SR1, SR2, p1, p2 = evaluate(AR,AL)
        df = pd.DataFrame([SD1, SD2, SR1, SR2, p1, p2])
        plt.figure(),
        plt.matshow(df.corr())
        plt.colorbar()
        
        ignored, p = stats.ttest(df)
        
        RMSE1 = np.sqrt((Residule/np.mean(AR))**2)/len(AR)


def PS(ref, rho, winL, T, prop):
#corner:5: fcp = 0.32*vs/rho  3 fcp = 1.53*vp/2pi*rho   
#    n = np.linspace(0, T, len(ref))
#    ref = y, T = Ty, prop = ''
    ref1 = ref
    freq = []
    if T % winL > 0:
        T = int(T/winL)*winL
#        ref = ref[0:T-1]
        ref = ref[0:int(T)]
    freq = np.fft.fftfreq(int(T))[0:int(T/2)]
    freq = freq[freq >= 0]
    t0 = np.linspace(0, T-winL, int(len(ref)/winL))
    fmin = min(freq)
    fmax = max(freq)
    PT = []
    v = []
    vi = []
    k = []
    Vg = []
    x = []
    a = []
    TP = []
    Phii = []
    vavg = 0
    fi = 1/winL
    for i in range(len(t0)): 
#    for i in range(197, len(t0)):
        print(i)
        t = t0[i]
        pti = np.zeros((winL, 1))
        if i+1 < len(t0):
            ti = np.array(range(int(t), int(t0[i+1])))
#        else: break
        p0 = ref[ti]
#        dt = (t0[i+1]-t)/winL
        Amax = max(p0)
        Amin = min(p0)
        A0 = p0[0]
        At = p0[-1]
        if abs(A0 - At)/abs(Amax - Amin) < 0.5:
            A = max(abs(p0))            
            pf = np.fft.fft(p0)
            wi= np.arctan(pf.real/pf.imag)
            Hi = (np.sign(ti-t)+1)/2
            pti = A*(ti-t)*np.exp(-2*(ti-t))
#            plt.figure(),
#            plt.plot(pti)
            pti *= np.cos(2.*np.pi*fi*(ti-t) + wi * np.pi) * Hi  
#            plt.figure(),
#            plt.plot(pti)
#        plot_tfr(pti, dt=dt, fmin=0, fmax=fi+winL*dt)
            pti += p0
        else:
            if len(t0) > i+2:
                ti2 = np.array(range(int(t), int(t0[i+2])))
            else:
                ti2 = np.array(range(int(t), int(2*t0[-1]-t0[-2])))
            p0 = ref[ti2]
            Amax = max(p0)
            Amin = min(p0)
            A0 = p0[0]
            At = p0[-1]
            A = max(abs(p0))            
            pf = np.fft.fft(p0)
            wi= np.arctan(pf.real/pf.imag)
            Hi = (np.sign(ti2-t)+1)/2
            pti = A*(ti2-t)*np.exp(-2*(ti2-t))
            pti *= np.cos(2.*np.pi*fi*(ti2-t) + wi * np.pi) * Hi  
            pti += p0
            pti = pti[0:winL]
            wi = wi[0:winL]
            ti = ti2
            #            x = pti
#            x = list(x)
        pti = np.nan_to_num(pti)
        wi = np.nan_to_num(wi)
        ti = np.nan_to_num(ti)
        if len(pti) < len(ti):
            ti = ti[0:len(pti)]
        vi = pti/(ti-t+0.000001)
#            v = list(v)
#        phii = list(phii)
        ki = wi[0:len(vi)]/ (vi+0.000001)
#        if abs(vi).all()< 0.00001:
#            x.append(0)
#            v.append(0)
#            Phii.append(0)
#            k.append(0)
#            PT.append(0)
#            Vg.append(0)
#            TP.append(0)
#            a.append(0)
#            print(vi)
#        else:
        print(ki)
        print(max(ki)-min(ki))
        vg = []
        ai = []
        if max(ki)-min(ki) < 2:
        #            k = np.array(k)
            if i >= len(Phii)-1:
                Phii.append(wi)
                phii = wi
                v.append(vi)
                k.append(phii/vi)
                x.append(pti)
                PT.append(pti)
#                if i > 1 & (len(Phii[i]) == len(Phii[i-1])):
            if i > 1 and (len(Phii[i]) == len(Phii[i-1])):  
                Phii1 = np.array(Phii)
                k1 = np.array(k)
                vg=(Phii1[i]-Phii1[i-1])/(k1[i]-k1[i-1])
                ai = [vi[0]-v[i-1][-1],vi[1]-vi[0],vi[2]-vi[1],vi[3]-vi[2],vi[4]-vi[3]]
            else:
                vg=np.array(Phii[i])/np.array(k[i])  
                ai = []
                ai.append(vi[0])
            if len(vi) > 1:
                for i in range(len(vi)-1):
                    ai.append(vi[i+1]-vi[i])
#                    D.append(x)
                    Vg.append(vg)
                    TP.append(ti)
                    a.append(ai)
        else:
            PT1 = PT
            v1 = v
            k1 = k
            Vg1 = Vg
            x1 = x
            a1 = a
            TP1 = TP
            Phii1 = Phii
            vavg1 = vavg
            print([PT1,v1,k1,Vg1,x1,a1,TP1,Phii1,vavg1])
            pti, vi, ki, vgi, xi, ai, ti, phii, vavg = PS(p0, rho, int(len(p0)/2), int(len(p0)), prop)
            print([pti, vi, ki, vgi, xi, ai, ti, phii, vavg])
                
            Phii.append(phii)
            v.append(vi)
            k.append(ki)
            x.append(pti)
            PT.append(pti)
            Vg.append(vgi)
            TP.append(ti)
            a.append(ai)
        if i == len(t0) - 1: break
#    PT.merge(method = 1)
    if prop == 'P':
        PT = np.array(PT[:]).ravel()
        plt.figure(),
        plt.subplot(311)
        plt.plot(PT)
        plt.ylabel('simulated Pwave')
        plt.title(['residule:', np.mean(np.array(PT) - ref1[0:len(np.array(np.array(PT[:]).ravel()))])]) 
        plt.subplot(312)
        plt.plot(ref1[0:len(np.array(np.array(PT[:]).ravel()))])
        plt.ylabel('original Pwave')
        plt.subplot(313)
        plt.plot(np.array(PT) - ref1[0:len(np.array(np.array(PT[:]).ravel()))])  
        plt.ylabel('residule')
    #    plot_tf_misfits(np.array(np.array(PT[:]).ravel()), ref1[0:len(np.array(np.array(PT[:]).ravel()))], dt=1, fmin=fmin, fmax=fmax)
    #    vavg = freq*2*np.pi*rho[sum(np.array(np.array(x).ravel()))]/1.53
        vavg = freq*2*np.pi*sum(np.array(np.array(x).ravel()))/1.53
    elif prop == 'S':
        PT = np.array(PT[:]).ravel()
        plt.figure(),
        plt.subplot(311)
        plt.plot(PT)
        plt.ylabel('simulated Swave')
        plt.title(['residule:', np.mean(np.array(PT) - ref1[0:len(np.array(np.array(PT[:]).ravel()))])]) 
        plt.subplot(312)
        plt.plot(ref1[0:len(np.array(np.array(PT[:]).ravel()))])
        plt.ylabel('original Swave')
        plt.subplot(313)
        plt.plot(np.array(PT) - ref1[0:len(np.array(np.array(PT[:]).ravel()))])  
        plt.ylabel('residule')
    #    plot_tf_misfits(np.array(np.array(PT[:]).ravel()), ref1[0:len(np.array(np.array(PT[:]).ravel()))], dt=1, fmin=fmin, fmax=fmax)
    #    vavg = freq*2*np.pi*rho[sum(np.array(np.array(x).ravel()))]/1.53
    #corner:6: fcs = 0.21*vs/rho
        vavg = freq/0.21
    else:
        PT = np.array(PT[:]).ravel()
        plt.figure(),
        plt.subplot(311)
        plt.plot(PT)
        plt.ylabel('simulated wave')
        plt.title(['residule:', np.mean(np.array(PT) - ref1[0:len(np.array(np.array(PT[:]).ravel()))])]) 
        plt.subplot(312)
        plt.plot(ref1[0:len(np.array(np.array(PT[:]).ravel()))])
        plt.ylabel('original wave')
        plt.subplot(313)
        plt.plot(np.array(PT) - ref1[0:len(np.array(np.array(PT[:]).ravel()))])  
        plt.ylabel('residule')
    return PT, v, k, Vg, x, a, TP, Phii, vavg 
#for t = int(s*fs+3*(s-p)*fs+1):int(s*fs+(3*(s-p)+6)*fs+1)
def Sf(x, y, z, rho, winL, Tx, Ty, Tz, iters):
#Sf(x = datasete[int(s*fs+3*(s-p)*fs+1):int(s*fs+(3*(s-p)+6)*fs+1)], y = datasetn[int(s*fs+3*(s-p)*fs+1):int(s*fs+(3*(s-p)+6)*fs+1)],  z= datasetz[int(s*fs+3*(s-p)*fs+1):int(s*fs+(3*(s-p)+6)*fs+1)], rho, 5, Tx = 3*fs, Ty = 3*fs, Tz = 3*fs, iters = 25)
    w0  = 80.897/180*np.pi
    phi0 = 15.754/180*np.pi
    Re = 6371
    dW = np.linspace(3, 12, int(9/1.5+1))
#    WW = []
    for dw in dW:
#        winL = dw/360*2*np.pi*Re
        ref = x
        ref1 = y
        ref2 = z
        #forward
        LT, vl, kl, Vgl, xy, al, TL, PhiLi, vavgl = PS(y,rho, winL, Ty, '')
        ZT, vz, kz, Vgz, xz, az, TZ, PhiZi, vavgz = PS(z,rho, winL, Tz, '')
        RT, vr, kr, Vgr, xx, ar, TR, PhiRi, vavgr = PS(x,rho, winL, Tx, '')
        r = np.sqrt(xy**2+xx**2)
        alpha = -np.log(y/r**0.5)/xz
        AL = r**(-0.5)*exp(-alpha*xz)
        AL,Y = np.meshgrid(xx,xy)
        plt.figure(),
        plt.subplot(211)
        plt.pcolor(X,Y,AL)
        plt.title('Love wave grad')
              
        #inversion
        accept, AR, Residule, dK, dTheta = MCMC_2D(x, y, Tx, alpha, iters)
        std_Residule = []
        std_dK = []
        std_dTheta =[]
        for i in range(len(AL)):
            std_Residule.append(3*np.std(Residuel[0:i]))
            std_dK.append(3*np.std(dK[0:i])) 
            std_dTheta.append(3*np.std(dTheta[0:i]))
        
        plt.subplot(212),
        plt.pcolor(X, Y, AR)
        plt.title('Rayleigh wave grad')
        
        #3D
        RT_3D, vl, kl, Vgl, xy, al, TL, PhiLi, vavgl = SPECPS_3D(xx, xy, rho, winL, Tx, 'L')
        
        
        SD1, SD2, SR1, SR2, p1, p2 = evaluate(AR,AL)
        df = pd.DataFrame([SD1, SD2, SR1, SR2, p1, p2])
        plt.figure(),
        plt.matshow(df.corr())
        plt.colorbar()
        
        ignored, p = stats.ttest(df)
        
        RMSE1 = np.sqrt((Residule/np.mean(AR))**2)/len(AR)
#        RMSE2 = np.sqrt((Residule/np.mean(AL))**2)/len(AL) 
              
        plt.figure(),
        plt.subplots(311)
        plt.plot(range(accept), Residule, 'r')
        plt.fillbetween(range(accept), Residule - std_Residule, Residule + std_Residule, color = 'gray', alpha = 0.05)
        plt.subplots(312)
        plt.plot(range(accept), dK, 'r')
        plt.fillbetween(range(accept), dK - std_dK, dK + std_dK, color = 'gray', alpha = 0.05)
        plt.subplots(313)
        plt.plot(range(accept), dTheta, 'r')
        plt.fillbetween(range(accept), dTheta - std_dTheta, dTheta + std_dTheta, color = 'gray', alpha = 0.05)
   
        ResultMatrix = [LT, ZT, RT, kl, kz, kr, Vgl, Vgz, Vgr, xy,xz,xx,al,az, ar, TL, TZ, TR, PhiLi, PhiZi, PhiRi, AL, AR] 
        corr = ResultMatrix.corr()
        
        ax = sns.heatmap(
            corr,
            vmin = -1, vmax=1,center =0,
            cmap = sns.diverging_pallette(20, 220, n = 200),
            square = True,
        )
        ax.set_xticklabels(
            ax.get_xticklabels(),
            rotation = 45,
            horizontalalignment = 'right'
        )
#        vr = 0.9*vsx
#        vl = (vs - vr) 
#       n = np.linspace(0, T, len(refx))
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""        
        Wi = np.linspace(w0, w02*np.pi, int(2*np.pi/dw+1))    
        Phi = np.linspace(phi0-np.pi, phi0+np.pi, int(2*np.pi/dw+1))
        dC = []
        C = []
        GCG = []
        CG = []
        for i in range(0, len(Wi)):
            Ci = 2*np.pi*Re/(i+0.5)/winL
            dC.append(-Ci/VgR[i]*winL/T)            
#            dcL = 0.88 + 2.3*10**(-15)*(np.cos(2*Phi) - np.sin(2*Phi)) + 8.6*10**(-4)*(np.cos(4*Phi) - np.sin(4*Phi))               
            dt = 13
            coe = (dt/winL)**2
            delta = 5
            stdm = np.std(wx)
            varm = stdm**2
            c = varm*np.exp(-dw**2/2/delta**2)
            C.append(c)
            gcg = varm*np.sqrt(2*np.pi)*delta/winL
            GCG.append(gcg)
            cg = varm*np.exp(-(i+0.5)**2*delta**2/2)
            CG.append(cg)
        WW.append(wx + C*GCG/(Cov(wx)+GCG)*(wx-Wi))  
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""      
         
        r = np.sqrt(xy**2+xx**2)
        alpha = -np.log(y/r**0.5)/xz
        AL = r**(-0.5)*exp(-alpha*xz)
        X,Y = np.meshgrid(xx,xy)
       
        plt.figure(),
        plt.subplot(211)
        plt.pcolor(X,Y,AL)
        
        plt.figure(),
        plt.subplot(211)
        plt.pcolor(X,Y,AL)        
        plt.title('Love wave grad')
        #inversion
        accept, AR, Residule, dK, dTheta = MCMC_2D(x, y, Tx, alpha, iters)
        std_Residule = []
        std_dK = []
        std_dTheta =[]
        for i in range(len(AL)):
            std_Residule.append(3*np.std(Residuel[0:i]))
            std_dK.append(3*np.std(dK[0:i])) 
            std_dTheta.append(3*np.std(dTheta[0:i]))
        

        plt.subplot(212),
        plt.pcolor(X, Y, AR)
        plt.title('Rayleigh wave grad')
        
        SD1, SD2, SR1, SR2, p1, p2 = evaluate(AR,AL)
        df = pd.DataFrame([SD1, SD2, SR1, SR2, p1, p2])
        plt.figure(),
        plt.matshow(df.corr())
        plt.colorbar()
        
        ignored, p = stats.ttest(df)
        
        RMSE1 = np.sqrt((Residule/np.mean(AR))**2)/len(AR)
#        RMSE2 = np.sqrt((Residule/np.mean(AL))**2)/len(AL) 
        
        
        plt.figure(),
        plt.subplots(311)
        plt.plot(range(accept), Residule, 'r')
        plt.fillbetween(range(accept), Residule - std_Residule, Residule + std_Residule, color = 'gray', alpha = 0.05)
        plt.subplots(312)
        plt.plot(range(accept), dK, 'r')
        plt.fillbetween(range(accept), dK - std_dK, dK + std_dK, color = 'gray', alpha = 0.05)
        plt.subplots(313)
        plt.plot(range(accept), dTheta, 'r')
        plt.fillbetween(range(accept), dTheta - std_dTheta, dTheta + std_dTheta, color = 'gray', alpha = 0.05)
   
        ResultMatrix = [LT, ZT, RT, kl, kz, kr, Vgl, Vgz, Vgr, xy,xz,xx,al,az, ar, TL, TZ, TR, PhiLi, PhiZi, PhiRi, AL, AR] 
        corr = ResultMatrix.corr()
        
        ax = sns.heatmap(
            corr,
            vmin = -1, vmax=1,center =0,
            cmap = sns.diverging_pallette(20, 220, n = 200),
            square = True,
        )
        ax.set_xticklabels(
            ax.get_xticklabels(),
            rotation = 45,
            horizontalalignment = 'right'
        )    


def evaluate(P1,P2,Pref,title1, title2, title3):
    P1 = np.array(P1)
    P2 = np.array(P2) 
    Pref = np.array(Pref) 
    SD1 = abs(P1 - Pref).ravel()/P1
    SD2 = ((P1-Pref)**2).ravel()/P1**2
    SR1 =abs(P1 - P2).ravel()/P2
    SR2 = ((P1-P2)**2).ravel()/P2**2
    
    ignored1, p1 = stats.ttest(P1, P2)
    ignored2, p2 = stats.ttest(P1, Pref)
    
    plt.figure,
    plt.subplot(221)
    plt.bar(range(len(SD1[0,])), SD1[0,])
    plt.plot(np.mean(SD1[0,]),0,'rx')
#    plt.title([title1, SD1[1,], SD1[2,],SD1[3,]])
    plt.title([title1, SD1[1,], SD1[2,]])
    plt.ylabel(['SD1',title2])
    plt.subplot(222)
    plt.bar(range(len(SD2[0,])), SD2[0,])
    plt.plot(np.mean(SD2[0,]),0,'rx')
#    plt.title([SD2[1,], SD2[2,],SD2[3,]])
    plt.title([SD2[1,], SD2[2,]])    
    plt.ylabel(['SD2',title2])
    plt.subplot(223)
    plt.bar(range(len(SR1[0,])), SR1[0,])
    plt.plot(np.mean(SR1[0,]),0,'rx')
#    plt.title([ SR1[1,], SR1[2,],SR1[3,]])
    plt.title([SR1[1,], SR2[2,]])
    plt.ylabel(['SR1',title3])
    plt.subplot(224)
    plt.bar(range(len(SR2[0,])), SR2[0,])
    plt.plot(np.mean(SR2[0,]),0,'rx')
#    plt.title([SR2[1,], SR2[2,],SR2[3,]])
    plt.title([SR2[1,], SR2[2,]])
    plt.ylabel(['SR2',title3])
    
    return SD1,SD2,SR1,SR2,p1,p2

    
def Gaussian(mu, sigma, G0):
    p0 = 1
    mu0 = 0
    sigma0 = 0
    for mui in np.linspace(mu-4,mu+5,10):
        for sigmai in np.linspace(sigma-1,sigma+1,10):
            G = np.random.normal(mui, sigmai,50)
            count0, bins0, ignored0 = plt.hist(G+ mean0, 50, density=True)
            y = np.exp(-(bins0-mu)/(2*sigma**2))/np.sqrt(2*pi*sigma**2)
            d,p = stats.ttest_ind(G[0:49,], G0[0:49,])
            if p< p0:
                p0 = p
                mu0 = mui
                sigma0 = sigmai
    return p0, mu0, sigma0, bins0, y    
    
def Gamma(k, theta, G0):
    p0 = 1
    theta0 = 0
    k0 = 0
    for ki in np.linspace(k-4,k+5,10):
        for thetai in np.linspace(theta-0.1,theta+0.1,10):
            G = np.random.gamma(ki, thetai,50)
            mean0 = np.mean(G) + ki*thetai
            count0, bins0, ignored0 = plt.hist(G + mean0, 50, density=True)
            y = (bins0-mean0)**(ki-1)*(np.exp(-(bins0-mean0)/thetai)/(sps.gamma(ki)*thetai**ki))
            d,p = stats.ttest_ind(G[0:49,], G0[0:49,])
            if p< p0:
                p0 = p
                k0 = ki
                theta0 = thetai
    return p0, k0, theta0, bins0, y    
#data = datf
    
def MCMC_2D(x, y, T, alpha, iters):
    Ay = abs(y)
    if sum(Ay > 1000) > 0:
        Ay[y > 0] = np.log2(y[y > 0])
        Ay[y <= 0] = np.log2(-y[y <= 0])
    k =(2/skew(Ay))**2
    theta = np.sqrt(np.var(Ay)/k)
    mean = np.mean(Ay)+k*theta
    count, bins, ignored = plt.hist(Ay+mean,50,density = 'TRUE')   

    accept = 0
    dK = []
    dTheta = []
    X = []
    Residule = []
    residule = 10000
    for i in range(iters):
#        pi,ki,thetai,binsi,yi = Gamma(k, theta, y)
        pi,ki,thetai,binsi,yi = Gamma(k, theta, bins)
        q = np.exp(np.log2(binsi) - np.log2(bins))
        A = min(1, np.mean(q))
        if (A > np.random.rand()): 
            accept += 1.0 
            rih = (x**2+yi**2)**0.25
            expZ_al = yi/rih
            xi = rih*expZ_al
            if sum(xi - x) < residule:
                residule = xi -x
                Residule.append(residule)
                X.append(xi)
                dK.append(ki - k)
                dTheta.append(thetai - theta)
        else:
            break 
    return accept, X, Residule, dK, dTheta

def MCMC_Amp(data, iters):
    # all permutation without order is L C N
#    row = factor(L)/(factor(N)*factor(N-L))
# Gamma prior
    N= len(data)
    data = data.real
    accept = 0
    dt_df = pd.DataFrame(data)
    #replace blank cells with mean
    dt_df.fillna(np.mean(dt_df.iloc[:,0]), inplace=True)
    dt_df = np.array(dt_df, dtype = float)  
#amplitude only positive    
    dt_df = abs(dt_df)
    
    k =(2/skew(dt_df))**2
    theta = np.sqrt(np.var(dt_df)/k)
    mean = np.mean(dt_df)+k*theta
    count, bins, ignored = plt.hist(dt_df+mean,50,density = 'TRUE')
#    bins = np.log2(bins)
#    mean = np.log2(mean)
#    k = np.log2(k)
#    y = (abs(bins - 2.1*mean)**k)**2/abs(bins - mean)**k * np.exp(-(bins - mean) / theta**2)* np.exp(theta) / (sps.gamma(2**(k-10))**22* (theta**k)**2*4*(10**62))

#    y = abs(bins - mean)**(k - 1)*np.exp(-(bins -mean) / theta)/(sps.gamma(k) * theta**k)   
#    plt.plot(y,linewidth=2, color='r') 
     
    K = []
    Theta = []
    Bins = [[]]
    Y = [[]]
    for i in range(iters):
#        pi,ki,thetai,binsi,yi = Gamma(k, theta, y)
        pi,ki,thetai,binsi,yi = Gamma(k, theta, bins)
        q = np.exp(np.log2(binsi) - np.log2(bins))
        A = min(1, np.mean(q))
        if (A > np.random.rand()): 
            accept += 1.0 
            K .append(ki)
            Theta.append(thetai)
            Bins.append(binsi)
            Y.append(yi)
        else:
            break
    K1 = np.mean(K)
    Theta1 = np.mean(Theta)
#    Bins1 = np.mean(Bins[Bins != []],axis = 0)
    Bins[np.shape(Bins) == 0] = Bins[-1]
    Bins1 = np.mean(Bins, axis = 0)
    Mean1 = np.mean(Bins1) + K1 * Theta1
    y = abs(Bins1 - Mean1)**(K1 - 1)*np.exp(-(Bins1 - Mean1) / Theta1)/(sps.gamma(K1) * Theta1**K1)
    Gs = np.random.gamma(K1, Theta1, 50)
    count_, bins_, ignored_ = plt.hist(Gs+ Mean1, 50, density = 'TRUE')
#    plt.plot(y, linewidth=2, color='r') 
    k_ = K1
    theta_ = Theta1
    mean_ = Mean1    
#
#use log scale 
#    if sum(abs(dt_df) > 1000) > 0:
#        dt_df[dt_df > 0] = np.log2(dt_df[dt_df > 0])
#        dt_df[dt_df <= 0] = -np.log2(-dt_df[dt_df <= 0])
        
    if sum(abs(dt_df) > 1000) > 0:
        dt_df[dt_df > 0] = np.log2(dt_df[dt_df > 0])
        dt_df[dt_df <= 0] = np.log2(-dt_df[dt_df <= 0])

        
    k1 =(2/skew(dt_df))**2
    theta1 = np.sqrt(np.var(dt_df)/k1)
    mean1 = np.mean(dt_df) + k1*theta1
    count1, bins1, ignored = plt.hist(dt_df + mean1,50,density = 'TRUE')
#    y1 = abs(bins1 - mean1)**(k1 - 1)*np.exp(-(bins1 - mean1) / theta1)/(sps.gamma(k1) * theta1**k1)
#    plt.plot(y1,linewidth=2, color='r')  
    K2 = 0
    Theta2 = 0
    Bins2 = 0
    for i in range(iters):
        pi1,ki1,thetai1,binsi1,yi1 = Gamma(k1, theta1, bins1)
        q = yi1/y1
        A = min(1, np.mean(q))
        if (A > np.random.rand()): 
            accept1 += 1.0 
            K2 += ki1
            Theta2 += thetai1
            Bins2 += binsi1
        else:
            break
    K2 /= accept
    Theta2 /= accept
    Bins2 /= accept
    Mean2 = np.mean(yi1) + K1 * Theta1
    K2 = 2**(K1)
    Theta2 = 2**(Theta2)  
    Bins2 = 2**(Bins2)  
    Mean2 = 2**(Mean2)  
    y2 = abs(Bins2 - Mean2)**(K2 - 1) * np.exp(-(Bins2 - Mean2) / Theta2) / (sps.gamma(K2) * Theta2**K2)
    Gs2 = np.random.gamma(K2, Theta2, 50)
    count2_, bins2_, ignored2_ = plt.hist(Gs2 + Mean2, 50, density = 'TRUE')
#    plt.plot(y1, linewidth=2, color='r')  
    k2_ = K2
    theta2_ = Theta2
    mean2_ = Mean2   
    SD1,SD2,SR1,SR2 = evaluate([bins_, k_, theta_, mean_], [bins2_, k2_, theta2_, mean2_], [bins, k, theta, mean], ' with: k: , theta, mean being: ', 'Amplitude comp Gamma and ref' , 'Amplitude comp Gamma and log ')    
#Gaussian prior with log scale       
    N= len(data)
    data = data.real
    accept = 0
    dt_df = pd.DataFrame(data)
    #replace blank cells with mean
    dt_df.fillna(np.mean(dt_df.iloc[:,0]), inplace=True)
    dt_df = np.array(dt_df, dtype = float)  
#amplitude only positive    
    dt_df = abs(dt_df)    

    sigma = np.std(dt_df)
    mu = np.mean(dt_df)
    count, bins, ignored = plt.hist(dt_df+mu,50,density = 'TRUE')
#    bins = np.log2(bins)
#    mean = np.log2(mean)
#    k = np.log2(k)
#    y = (abs(bins - 2.1*mean)**k)**2/abs(bins - mean)**k * np.exp(-(bins - mean) / theta**2)* np.exp(theta) / (sps.gamma(2**(k-10))**22* (theta**k)**2*4*(10**62))

    y = np.exp(-(bins -mean)**2 / (2*sigma**2))/(np.sqrt(2*np.pi*sigma**2))   
    plt.plot(y,linewidth=2, color='r') 
     
    Mu = []
    Sigma = []
    Bins = [[]]
    Y = [[]]
    for i in range(iters):
        pi,mui,sigmai,binsi,yi = Gaussian(mu, sigma, bins)
        q = np.exp(np.log2(binsi) - np.log2(bins))
        A = min(1, np.mean(q))
        if (A > np.random.rand()): 
            accept += 1.0 
            Mu.append(mui)
            Sigma.append(sigmai)
            Bins.append(binsi)
            Y.append(yi)
        else:
            break
    Mu1 = np.mean(Mu)
    Sigma1 = np.mean(Sigma)
#    Bins1 = np.mean(Bins[Bins != []],axis = 0)
    Bins[np.shape(Bins) == 0] = Bins[-1]
    Bins1 = np.array(np.mean(Bins, axis = 0))
    Mean1 = np.mean(Bins1) 
#    y = abs(Bins1 - Mean1)**(K1 - 1)*np.exp(-(Bins1 - Mean1) / Theta1)/(sps.gamma(K1) * Theta1**K1)
    y = np.exp(-(Bins1 -Mean1)**2 / (2*Sigma1**2))/(np.sqrt(2*np.pi*Sigma1**2)) 
    Gs = np.random.normal(Mu1, Sigma1, 50)
    count_, bins_, ignored_ = plt.hist(Gs+ Mean1, 50, density = 'TRUE')
    plt.plot(y, linewidth=2, color='r') 
    mu_ = Mu1
    sigma_ = Sigma1
    mean_ = Mean1    
#
#use log scale 
#    if sum(abs(dt_df) > 1000) > 0:
#        dt_df[dt_df > 0] = np.log2(dt_df[dt_df > 0])
#        dt_df[dt_df <= 0] = -np.log2(-dt_df[dt_df <= 0])
        
    if sum(abs(dt_df) > 1000) > 0:
        dt_df[dt_df > 0] = np.log2(dt_df[dt_df > 0])
        dt_df[dt_df <= 0] = np.log2(-dt_df[dt_df <= 0])
        
    sigma1 = np.std(dt_df)
    mu1 = np.mean(dt_df)
    count1, bins1, ignored1 = plt.hist(dt_df+mu1,50,density = 'TRUE')
#    y1 = abs(bins1 - mean1)**(k1 - 1)*np.exp(-(bins1 - mean1) / theta1)/(sps.gamma(k1) * theta1**k1)
#    plt.plot(y1,linewidth=2, color='r')  
    Mu2 = []
    Sigma2 = []
    Bins2 = [[]]
    Y = [[]]
    accept = 0
    for i in range(iters):
        pi1,mui1,sigmai1,bins1,yi1 = Gaussian(mu1, sigma1, bins1)
        q = yi1/y
        A = min(1, np.mean(q))
        if (A > np.random.rand()): 
            accept += 1.0 
            Mu2.append(mui)
            Sigma2.append(sigmai1)
            Bins2.append(bins1)
        else:
            break
    Mu2 = np.mean(Mu2)
    Sigma2 = np.mean(Sigma2)
#    Bins1 = np.mean(Bins[Bins != []],axis = 0)
    Bins2[np.shape(Bins2) == 0] = Bins2[-1]
    Bins2 = np.array(np.mean(Bins2, axis = 0))
    Mean2 = np.mean(Bins2) 
    Sigma2 = 2**(Sigma2)
    Mu2 = 2**(Mu2)  
    Bins2 = 2**(Bins2)  
    Mean2 = 2**(Mean2)  
    y2 = np.exp(-(Bins2 -Mean2)**2 / (2*Sigma2**2))/(np.sqrt(2*np.pi*Sigma2**2)) 
    Gs2 = np.random.normal(Mean2, Sigma2, 50)
    count2_, bins2_, ignored2_ = plt.hist(Gs2, 50, density = 'TRUE')
    plt.plot(y2, linewidth=2, color='r')  
    mu2_ = Mu2
    sigma2_ = Sigma2
    mean2_ = Mean2   
    SD1N,SD2N,SR1N,SR2N = evaluate([abs(bins_), mu_, sigma_], [y2*10**8, mu2_, sigma2_], [bins, mu, sigma], ' with: mu, sigma being: ', 'Amplitude comp Gaussian and ref' , 'Amplitude comp Gaussian and log ')    


#rupture data
import obspy
from obspy import read
from obspy import *
import matplotlib.pyplot as plt
import numpy as np
from  scipy.stats import * 
import scipy.special as sps
import pandas as pd
from obspy.signal.trigger import ar_pick
from obspy.signal.tf_misfit import plot_tf_misfits

#st = read('http://examples.obspy.org/RJOB_061005_072159.ehz.new')
#data = st[0]
#print(data)
#print(data.stats) 
#dataset = data.data
#npts = data.stats.npts
#t = np.linspace(1,npts,npts)
#fs = data.stats.sampling_rate
#T = npts/fs


#datf = np.fft.fft(dataset)
#freq = np.fft.fftfreq(t.shape[-1])
#w = np.arctan(np.real(datf)/(np.imag(datf)+0.0001))
#Amax = max(datf)
#plt.plot(datf)

ch1 = read('https://examples.obspy.org/loc_RJOB20050801145719850.z.gse2')[0]
ch2 = read('https://examples.obspy.org/loc_RJOB20050801145719850.n.gse2')[0]
ch3 = read('https://examples.obspy.org/loc_RJOB20050801145719850.e.gse2')[0]
dataz = ch1
datan = ch2
datae = ch3
print(dataz)
print(dataz.stats) 
#w  = 12.80, phi = 47.74/180*np.pi
datasetz = dataz.data
datasetn = datan.data
datasete = datae.data

npts = dataz.stats.npts
t = np.linspace(1,npts,npts)
fs = dataz.stats.sampling_rate
T = npts/fs

# plot the raw and filtered data...
t = np.arange(0, npts / fs, dataz.stats.delta)
plt.figure,
plt.subplot(311)
plt.plot(t, datasetz, 'k')
plt.ylabel('Z wave')
plt.title('Raw Data(time domain)')
plt.subplot(312)
plt.plot(t, datasetn, 'k')
plt.ylabel('N wave')
plt.subplot(313)
plt.plot(t, datasete, 'k')
plt.ylabel('E wave')

dataz.filter('highpass', freq=1.6, corners=5, zerophase=True)

datf = np.fft.fft(dataz)
freq = np.fft.fftfreq(t.shape[-1])
datf1 = datf
freq1 = freq
datf = datf[0:int(0.5*len(freq)-1)]
freq = freq[0:int(0.5*len(freq)-1)]
w = np.arctan(np.real(datf)/(np.imag(datf)+0.0001))
Amax = max(datf)
plt.plot(datf)

plt.subplot(211)
plt.plot(t, datasetz, 'k')
plt.ylabel('Highpassed Data(Baseline correction)')
plt.xlabel('Time [s]')
plt.subplot(212)
plt.plot(freq, datf.real, freq, datf.imag)
#plt.plot(t, datS_filt.data, 'k')
plt.ylabel(['real', 'imag'])
plt.title('specV.s.frequency')
plt.show()

#characteristics(Amax in mm, D in km, Krishna River
#p, s = ar_pick(a = dataset, samp+rate= fs, f1=1.0, f2=20.0, lta_p = 1, sta_p=0.1 ,lta_s =4 , sta_s=1 , m_p=1, m_s=1, l_p=2, l_s=8)
p, s = ar_pick(datasetz, datasetn, datasete, fs,
                         1.0, 20.0, 1.0, 0.1, 4.0, 1.0, 2, 8, 0.1, 0.2)


#plt.figure,
#plt.subplot(411)
##datf.filter('lowpass', freq=4.0, corners=2, zerophase=True)
#plt.plot(freq, datf)
#plt.ylabel('Baseline correction ')
#plt.title('Highpassed ZWave Amplitude(frequency domain)')

#corner:5: fcp = 0.32*vs/rho  3 fcp = 1.53*vp/2pi*rho
TP = 5
FP =1/TP
CP=FP/1.53
TS = 20
FS = 1/TS
CS=FS/0.21
#datP_filt = dataz.copy()
#datP_filt.filter('bandpass', freqmin=1, freqmax=3.6,corners=6, zerophase=True)

plt.figure,
plt.subplot(311)
#plt.plot(freq[0:int((s-p)*fs)],dat_filt[int(p*fs):int(p*fs+(s-p)*fs)])
plt.title('Amplitude(time domain)')
Pwave = dataz[int(p*fs):int(p*fs+(s-p)*fs)]
#Pw = w[int(p*fs):int(p*fs+(s-p)*fs)]
plt.plot(t[int(p*fs):int(p*fs+(s-p)*fs)],dataz[int(p*fs):int(p*fs+(s-p)*fs)])
plt.ylabel('P wave only')

#corner:6: fcs = 0.21*vs/rho
#datS_filt = dataz.copy()
#datS_filt.filter('highpass', freq=1.8, corners=5, zerophase=True)

plt.subplot(312)
#plt.plot(t[int((s-p)*fs)+1:int((s-p)*fs)+int(20*fs)],dat_filt[int((s-p)*fs)+1:int((s-p)*fs)+int(20*fs)])
#plt.ylabel('S wave Amplitude(time domain)')
Swave = dataz[int(s*fs+1):int(s*fs+3*(s-p)*fs)]
#Sw = w[int(s*fs+1):int(s*fs+3*(s-p)*fs)]
plt.plot(t[int(s*fs+1):int(s*fs+3*(s-p)*fs)],dataz[int(s*fs+1):int(s*fs+3*(s-p)*fs)])
plt.ylabel('S wave')

plt.subplot(313)
Sfwave = dataz[int(s*fs+3*(s-p)*fs+1):int(s*fs+(3*(s-p)+6)*fs+1)]
#Sfw = w[int(s*fs+3*(s-p)*fs+1):int(s*fs+(3*(s-p)+6)*fs+1)]
plt.plot(t[int(s*fs+3*(s-p)*fs+1):int(s*fs+(3*(s-p)+6)*fs+1)],dataz[int(s*fs+3*(s-p)*fs+1):int(s*fs+(3*(s-p)+6)*fs+1)])
plt.ylabel('Surface wave')

#########################################frequency domain
plt.figure,
plt.subplot(311)
#plt.plot(freq[0:int((s-p)*fs)],dat_filt[int(p*fs):int(p*fs+(s-p)*fs)])
plt.title('Spectrum(frequency domain)')
Pwave = dataz[int(p*fs):int(p*fs+(s-p)*fs)]
Pfreq = abs(np.fft.fftfreq(int(p*fs+(s-p)*fs))[int(p*fs+(s-p)*fs):int(p*fs):-1])
Pw = np.arctan(datf[2:2+len(Pfreq)].real/datf[2:2+len(Pfreq)].imag)
plt.plot(Pfreq, datf[2:2+len(Pfreq)].real)
plt.ylabel('P wave only')

#corner:6: fcs = 0.21*vs/rho
#datS_filt = dataz.copy()
#datS_filt.filter('highpass', freq=1.8, corners=5, zerophase=True)

plt.subplot(312)
#plt.plot(t[int((s-p)*fs)+1:int((s-p)*fs)+int(20*fs)],dat_filt[int((s-p)*fs)+1:int((s-p)*fs)+int(20*fs)])
#plt.ylabel('S wave Amplitude(time domain)')
Swave = dataz[int(s*fs+1):int(s*fs+3*(s-p)*fs)]
Sfreq = abs(np.fft.fftfreq(int(s*fs+3*(s-p)*fs))[int(s*fs+3*(s-p)*fs):int(s*fs+1):-1])
Sw = np.arctan(datf[2:2+len(Sfreq)].real/datf[2:2+len(Sfreq)].imag)
plt.plot(Sfreq, datf[2:2+len(Sfreq)].real)
plt.ylabel('S wave')

plt.subplot(313)
Sfwave = dataz[int(s*fs+3*(s-p)*fs+1):int(s*fs+(3*(s-p)+6)*fs+1)]
Sffreq = abs(np.fft.fftfreq(int(s*fs+(3*(s-p)+6)*fs+1))[int(s*fs+(3*(s-p)+6)*fs+1):int(s*fs+3*(s-p)*fs+1):-1])
Sfw = np.arctan(datf[2:2+len(Sffreq)].real/datf[2:2+len(Sffreq)].imag)
plt.plot(Sffreq, datf[2:2+len(Sffreq)].real)
plt.ylabel('Surface wave')

#parameters: 
#regional S and surfacewaves:85-105,1Hz/0.05Hz
D0 = np.linspace(1,200,2000)
Ml0 = np.log10(Amax/0.01)-0.15+1.6*np.log10(D0)
Ms0 = np.log10(Amax/20) +1.66*np.log10(D0)+2
M00 = abs(np.exp(1.5*Ml0+9))
MAX0 = max(M00)
MIN0 = min(M00)
#np.sqrt(MAX1/
#    Mb1 = np.log10(Amax/1) + Qp
#scaling law
du0 = (M00-1.6*10**15)/(5*10**16-1.6*10**15)*(0.1-0.04)+0.04 
Dref0 = (M00-1.6*10**15)/(5*10**16-1.6*10**15)*(4.7-1.5)+1.5
TR0 = (M00-1.6*10**15)/(5*10**16-1.6*10**15)*(1.4-0.4)+0.4

D1 = np.linspace(201,600,4000)
Ml1 = np.log10(Amax*10**6)-3.38+3*np.log10(D1)
Ms1 = np.log10(Amax/20) +1.66*np.log10(D1)+2
M01 = abs(np.exp(1.5*Ml1+9))
MAX1 = max(M01)
MIN1 = min(M01)

#    Mb1 = np.log10(Amax/1) + Qp
#scaling law
du1 = (M01-1.6*10**15)/(5*10**16-1.6*10**15)*(0.1-0.04)+0.04 
Dref1 = (M01-1.6*10**15)/(5*10**16-1.6*10**15)*(4.7-1.5)+1.5
TR1 = (M01-1.6*10**15)/(5*10**16-1.6*10**15)*(1.4-0.4)+0.4

D2 = np.linspace(601,2750,20000)
Ml2 = np.log10(Amax*10**6)-2.48+2.76*np.log10(D2)
Ms2 = np.log10(Amax/20) +1.66*np.log10(D2)+2
M02 = abs(np.exp(1.5*Ml2+9)) 
MAX2 = max(M02)
MIN2 = min(M02)
#du2 = (M02-5*10**16)/(1.6*10**18-5*10**16)*(0.4-0.1)+0.1
#Dref2 = (M02-5*10**16)/(1.6*10**18-5*10**16)*(15-4.7)+4.7
#TR2 = (M02-5*10**16)/(1.6*10**18-5*10**16)*(4.5-1.4)+1.4
du2 = (M01-1.6*10**15)/(5*10**16-1.6*10**15)*(0.1-0.04)+0.04 
Dref2 = (M01-1.6*10**15)/(5*10**16-1.6*10**15)*(4.7-1.5)+1.5
TR2 = (M01-1.6*10**15)/(5*10**16-1.6*10**15)*(1.4-0.4)+0.4


D3 = np.linspace(2751,5000,22500)
Ml3 = np.log10(Amax*10**6)-2.48+2.76*np.log10(D3)
Ms3 = np.log10(Amax/20) +1.66*np.log10(D3)+2    
M03 = abs(np.exp(1.5*Ml3+9)) 
MAX3 = max(M03)
MIN3 = min(M03)
#du3 = (M03-1.6*10**18)/(5*10**19-1.6*10**18)*(1.3-0.4)+0.4 
#Dref3 = (M03-1.6*10**18)/(5*10**19-1.6*10**18)*(47.4-15)+15
#TR3 = (M03-1.6*10**18)/(5*10**19-1.6*10**18)*(14.4-4.5)+4.5
du3 = (M01-1.6*10**15)/(5*10**16-1.6*10**15)*(0.1-0.04)+0.04 
Dref3 = (M01-1.6*10**15)/(5*10**16-1.6*10**15)*(4.7-1.5)+1.5
TR3 = (M01-1.6*10**15)/(5*10**16-1.6*10**15)*(1.4-0.4)+0.4

D4 = np.linspace(5001,6400,14000)
Ml4 = np.log10(Amax*10**6)-2.48+2.76*np.log10(D4)
Ms4 = np.log10(Amax/20) +1.66*np.log10(D4)+2   
M04 = abs(np.exp(1.5*Ml4+9)) 
MAX4 = max(M04)
MIN4 = min(M04)
#du4 = abs((M04-5*10**19)/(1.6*10**21-5*10**19)*(4.3-1.3)+1.3).astype(float) 
#Dref4 = abs((M04-5*10**19)/(1.6*10**21-5*10**19)*(150-47.4)+47.4).astype(float) 
#TR4 = abs((M04-5*10**19)/(1.6*10**21-5*10**19)*(45.4-14.4)+14.4).astype(float) 
du4 = (M01-1.6*10**15)/(5*10**16-1.6*10**15)*(0.1-0.04)+0.04 
Dref4 = (M01-1.6*10**15)/(5*10**16-1.6*10**15)*(4.7-1.5)+1.5
TR4 = (M01-1.6*10**15)/(5*10**16-1.6*10**15)*(1.4-0.4)+0.4

plt.plot((du4/TR4).real)
    #bodywave:80-85,0.1-1Hz
    
    #Coda: 105-140
    
    #surfacewave(slowness model)
    
    #T=range()
    #Amplitude   
MCMC_Amp(datf, 1000)   
rho = [np.linspace(2.21*10**15, 3*10**15, 30).ravel(), np.linspace(3.4*10**15, 4.4*10**15, 720).ravel(), np.linspace(4.4*10**15, 5.6*10**15, 2171).ravel(), np.linspace(9.9*10**15, 12.2*10**15, 2259).ravel(), np.linspace(12.8*10**15, 13.1*10**15, 1221).ravel()]
T = range(int(p*fs),int(p*fs+(s-p)*fs))
winL = int(len(T)/24)
#PT,v,k,Vg,x,a,TP,Phii,varg=P(Pwave, rho, winL, len(T))
PT, xp, vp, kp, Vgp, vavgp, TP  = PS(dataz[T], rho, winL, len(T),'P')

T = range(int(s*fs+1), int(s*fs+3*(s-p)*fs))
ref = Swave - PT[(s-p)*fs+1: int(s*fs+3*(s-p)*fs)]
#freq = np.fft.fft(ref) 
#w = np.arctan(np.real(freq) / (np.imag(freq)+0.0001))
#ST, xs, vs, ks, vgs, vavgs, TS = PS(ref, rho, winL, len(T),'S')
#D = max(DP,DS)
Tbody = int((s-p)*fs) + TS
Sf(datasete[int(s*fs+3*(s-p)*fs+1):int(s*fs+(3*(s-p)+6)*fs+1)], datasetn[int(s*fs+3*(s-p)*fs+1):int(s*fs+(3*(s-p)+6)*fs+1)], datasetz[int(s*fs+3*(s-p)*fs+1):int(s*fs+(3*(s-p)+6)*fs+1)], rho, 5, 3*fs, 3*fs, 3*fs, 25)
SPECP()

refx = datasete.filter('highpass', freq=1.6, corners=5, zerophase=True) 
freqx = np.fft.fft(refx) 
wx = np.arctan(np.real(freqx) / (np.imag(freqx)+0.0001))
refy = datasete.filter('highpass', freq=1.6, corners=5, zerophase=True) 
freqy = np.fft.fft(refy) 
wy = np.arctan(np.real(freqy) / (np.imag(freqy)+0.0001))
refz = datasete.filter('highpass', freq=1.6, corners=5, zerophase=True) 
freqz = np.fft.fft(refz) 
wz = np.arctan(np.real(freqz) / (np.imag(freqz)+0.0001))



