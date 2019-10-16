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
    fi = 1/winL# change for more spectra
    
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
#amplitude
        if abs(A0 - At)/abs(Amax - Amin)<1.0:
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


def PS(ref, winL, T, prop):
#corner:5: fcp = 0.32*vs/rho  3 fcp = 1.53*vp/2pi*rho   
#    n = np.linspace(0, T, len(ref))
#    ref = y, T = Ty, prop = ''
    ref1 = ref
    freq = []
    if len(T) % winL > 0:
        T = int(len(T)/winL)*winL
#        ref = ref[0:T-1]
        ref = ref[0:int(T)]
    freq = np.fft.fftfreq(int(T))[0:int(T/2)]
    freq = freq[freq >= 0]
    t0 = np.linspace(0, T-winL, int(int(T)/winL))
    fmin = min(freq)
    fmax = max(freq)
    
    #Damping ratio kse
    kse = 0.05    
    wii = 2*np.pi/winL
    m = 30/32.2 #mass of the frame kips per foot
    c = 2*m*wii*kse
    wd =wii*(1-kse**2)**0.5
    PT = []
    v = []
    vi = []
    v2 = []#displcement
    k = []
    Vg = []
    x = []
    a = []#speed
    dv = []
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
            ti = range(int(t), int(t0[i+1]))
#        else: break
        p0 = ref[ti]
        v0 = p0[0]
        dv0 = p0[1] - p0[0]
        wd = wii*(1-kse**2)**0.5
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
        vi2 = (v0*np.cos(wd*(ti-t))+((dv0+v0)*kse*wii/wd)*np.sin(wd*(ti-t)))*np.exp(-kse*wii*(ti-t))
        print(vi- vi2)
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
                v2.append(vi2)
                k.append(phii/vi)
                x.append(pti)
                PT.append(pti)
#                if i > 1 & (len(Phii[i]) == len(Phii[i-1])):
            if i > 1 and (len(Phii[i]) == len(Phii[i-1])):  
                Phii1 = np.array(Phii)
                k1 = np.array(k)
                vg=(Phii1[i]-Phii1[i-1])/(k1[i]-k1[i-1])
#                ai = [vi[0]-v[i-1][-1],vi[1]-vi[0],vi[2]-vi[1],vi[3]-vi[2],vi[4]-vi[3]]
                ai = [vi2[0]-v[i-1][-1],vi2[1]-vi2[0],vi2[2]-vi2[1],vi2[3]-vi2[2],vi2[4]-vi2[3]]
            else:
                vg=np.array(Phii[i])/np.array(k[i])  
                ai = []
                ai.append(vi2[0])
            if len(vi2) > 1:
                for i in range(len(vi2)-1):
 #                   ai.append(vi[i+1]-vi[i])
                     ai.append(vi2[i+1]-vi2[i])
#                    D.append(x)
                     Vg.append(vg)
                     TP.append(ti)
                     a.append(ai)
        else:
            PT1 = PT
            v1 = v2
            k1 = k
            Vg1 = Vg
            x1 = x
            a1 = a
            TP1 = TP
            Phii1 = Phii
            vavg1 = vavg
            print([PT1,v1,k1,Vg1,x1,a1,TP1,Phii1,vavg1])
            pti, vi, ki, vgi, xi, ai, ti, phii, vavg = PS(p0, int(len(p0)/2), int(len(p0)), prop)
            print([pti, vi, ki, vgi, xi, ai, ti, phii, vavg])
                
            Phii.append(phii)
            v.append(vi)
            k.append(ki)
            x.append(pti)
            PT.append(pti)
            Vg.append(vgi)
            TP.append(ti)
            a.append(ai)
            an = np.unique(a)
        if i == len(t0) - 1: break
#    PT.merge(method = 1)
    if prop == 'P':
        PT = np.array(PT[:]).ravel()
        plt.figure(),
        plt.subplot(311)
        plt.plot(np.array(PT).ravel())
        plt.ylabel('simulated Pwave')
        plt.title(['residule:', np.mean(np.array(PT).ravel() - ref1[0:len(np.array(np.array(PT[:]).ravel()))])]) 
        plt.subplot(312)
        plt.plot(ref1[0:len(np.array(np.array(PT[:]).ravel()))])
        plt.ylabel('original Pwave')
        plt.subplot(313)
        plt.plot(np.array(PT).ravel() - ref1[0:len(np.array(np.array(PT[:]).ravel()))])  
        plt.ylabel('residule')
        plt.tight_layout()
        plt.figure(),
        plt.subplot(311)
        plt.plot(np.array(np.sum(a,axis = 0))[5:-1]-np.array(np.sum(a,axis = 0))[0:-6])
        plt.title('Pwave accelerate')
        plt.subplot(312)
        plt.plot(np.array(v2).ravel())
        plt.title('Pwave displacement')
        plt.subplot(313)
        plt.plot(np.sum(a,axis = 0))
        plt.title('Pwave speed')
        plt.tight_layout()
    #    plot_tf_misfits(np.array(np.array(PT[:]).ravel()), ref1[0:len(np.array(np.array(PT[:]).ravel()))], dt=1, fmin=fmin, fmax=fmax)
    #    vavg = freq*2*np.pi*rho[sum(np.array(np.array(x).ravel()))]/1.53
        vavg = freq*2*np.pi*sum(np.array(np.array(x).ravel()))/1.53
    elif prop == 'S':
        PT = np.array(PT[:]).ravel()
        plt.figure(),
        plt.subplot(311)
        plt.plot(np.array(PT).ravel())
        plt.ylabel('simulated Swave')
        plt.title(['residule:', np.mean(np.array(PT).ravel() - ref1[0:len(np.array(np.array(PT[:]).ravel()))])]) 
        plt.subplot(312)
        plt.plot(ref1[0:len(np.array(np.array(PT[:]).ravel()))])
        plt.ylabel('original Swave')
        plt.subplot(313)
        plt.plot(np.array(PT).ravel() - ref1[0:len(np.array(np.array(PT[:]).ravel()))])  
        plt.ylabel('residule')
        plt.tight_layout()
        plt.figure(),
        plt.subplot(311)
        plt.plot(np.array(np.sum(a,axis = 0))[4:-1]-np.array(np.sum(a,axis = 0))[0:-5])
        plt.title('Swave accelerate')
        plt.subplot(312)
        plt.plot(np.array(v2).ravel())
        plt.title('Swave displacement')
        plt.subplot(313)
        plt.plot(np.sum(a,axis = 0))
        plt.title('Swave speed')
        plt.tight_layout()
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
        RT_3D, vl, kl, Vgl, xy, al, TL, PhiLi, vavgl = SPECPS(xx, xy, rho, winL, Tx, 'L')
       
       
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


def MCMC_1D_improved(x, winL, T, rho, fs, domain):
#corner:5: fcp = 0.32*vs/rho  3 fcp = 1.53*vp/2pi*rho   
#    n = np.linspace(0, T, len(ref))
    if domain == 'time':
        mode = 1
        xx1 = x
        T = int(len(xx1))
        if T % winL > 0:
            T = np.ceil(T/winL)*winL
            xx = np.repeat(0.0, int(T))
            xx[0:len(x)] = xx1
            xx[int(len(x)+1):-1] = xx1[0:int(T-len(x)-1)]
        else:
            xx = x
#    else:
#        mode = 0
#        T = int(len(xx1)/2)
#        if T % winL > 0:
#            T = np.ceil(T/winL)*winL
#        freq = np.fft.fftfreq(int(T))
#        px = np.fft.fft(x)[-len(freq)-1:-1]
#        vx = x[np.array(np.linspace(0, len(xx1)-2, T),dtype = int)]

#        Re = 6371
#    WW = []
#    for dw in dW:
#    winL = dw/360*2*np.pi*Re               
    t0 = np.linspace(0, T, np.ceil(T/winL)+1)
    #    fmin = min(freq)
    #    fmax = max(freq)
    # Housner's design
#   fundamental mode
#    for kse in [0, 0.02, 0.05, 0.1]:
#    S0 = 0.0063
#    Tau = np.repeat(0,(len(t0)-1)* winL)
#    Alpha =np.repeat(0.0,(len(t0)-1)* winL)
    ACC = np.repeat(0.0,(len(t0)-1)* winL)
#    Pti = np.repeat(0.0,(len(t0)-1)* winL)
    Ve = np.repeat(0.0,(len(t0)-1)* winL)
    Di =np.repeat(0.0,(len(t0)-1)* winL)
#    AC = np.repeat(0,(len(t0)-1)* winL)
    Amp = np.repeat(0.0,(len(t0)-1)* winL)
    wii = 2*np.pi/winL 
    w0 = 0
    Spd = np.repeat(0.0, winL)
    Spv= np.repeat(0.0, winL)
    Spa = np.repeat(0.0, winL)
    res = np.repeat(0.0, winL)    
     #   fi = 1/winL# change for more spectra
     #   fii = np.random.normal(fi, 0.01*fi, 5)
    for tau in range(len(t0)-1): 
        temp = 1000
        print(tau)
        t = t0[tau]
        ti = np.array(range(int(t), int(np.ceil(t0[tau+1]))))
#        if mode == 0:
#            p0 = px[ti]
        if mode == 1:
            p0 = xx[ti]

        dt = 1/fs

        s0 = p0[0]
        v0 = p0[1]-p0[0]- v0 
        gamma = 1
        fi = 1/winL
        Amax = max(abs(p0))
        Amin = min(p0)
        if mode == 1:
            wi= -np.arccos(p0/Amax)
#        if mode == 0:
#            wi = -np.arctan(p0.imag/p0.real)
        vgi = wi[0]-w0
        w0 = wi[0]
        kse0 = 0
        wd = wi*(1-kse0**2)**0.5
#        if mode == 0:
#            thetai =np.arcsin(vx[ti]/Amax)[-1]- np.arcsin(vx[ti]/Amax)[0]
#            B = Amax**2/2*(1-(np.sin(2*thetai)*np.cos(2*w0*ti)-(1-np.cos(2*thetai))*np.sin(2*w0*ti))/(2*thetai))
#            B_ =  B**0.5
#            B_[B_*p0<0] = -B_[B_*p0<0]
#            
#            spv = np.repeat(0.0, winL)
#            for n in range(len(ti)):
#                pv = 0
#                for w in wi:
#                    tau = np.linspace(0, ti[n], (ti[n]+1)/dt)
#                    if pv < abs(sum(np.sin(w*(ti[n]-tau))*dt*np.exp(-kse*(w*(ti[n]-dt))))):
#                        pv = sum(np.sin(w*(ti[n]-tau))*dt*np.exp(-kse*(w*(ti[n]-dt))))
#                spv[n] = np.max(pv*B_).real
#    #            print(['Spv:', Spv])
#
#            multi = wd**2+(2*wd*kse0)**2/2*(np.cos(2*wd*ti))**2+kse0*wi*wd*np.cos(2*wd*ti)
#            multi *= 1-np.exp(-2*wi*kse0*ti)/wd**2
#            B_cor = abs(np.pi*wi*B_**0.5/(2*kse0*ki**2+0.000001)*multi)**0.5
#    #            B_cor = B_cor**0.5
#            B_cor[B_cor*p0 < 0] = -B_cor[B_cor*p0 < 0]
#            
#            Sum_B_cor = B_cor
#            Spd = spv/wi
#            Sum_Spd = spv/wi
#            Spv = spv
#            Sum_Spv = spv
#            Spa = spv*wi
#            Sum_Spa = spv*wi
#            diff = sum(abs(B_cor - p0))
        if mode == 1:
#             alpha = p0/dt
#            A0 = v0/wi**2 - 2*kse0*alpha/wi**3
#            A1 = p0/dt/wi**2
#            A2 = v0 - A0
#            A3 = 1/wd*(w0+kse0*wi*A2-alpha/wi**2)
#            spd = A0 + A1*ti + A2*np.exp(-kse0*wi*ti)*np.cos(wd*ti) + A3*np.exp(-kse0*wi*ti)*np.sin(wd*ti)
#            spv = A1 +(wd*A3 - kse0*wi*A2)*np.exp(-kse0*wi*ti)*np.cos(wd*ti) - (wd*A2 + kse0*wi*A3)*np.exp(-kse0*wi*ti)*np.sin(wd*ti) 
#            spa = A1 +(wd*A3 - kse0*wi*A2)*(-kse0*wi*np.exp(-kse0*wi*ti)*np.cos(wd*ti)-wd*ti*np.sin(wd*ti)*np.cos(wd*ti)) - (wd*A2 + kse0*wi*A3)*(-kse0*wi*np.exp(-kse0*wi*ti)*np.cos(wd*ti)-wd*ti*np.sin(wd*ti)*np.cos(wd*ti))
#            Sum_Spd = spd
#            Sum_Spv = spv
#            Sum_Spa = spa
#            Spd = spd
#            Spv = spv
#            Spa = spa
#            diff = sum(abs(Spd - p0))
             k0 = abs(p0[1]-p0[0])/dt
             c0 = 60*4.45
             fD = 0
             fS = 0
             m = 30*4.45
             a = (p0[0] - fD - fS)/m
             kd = k0 + 3*c0/dt +6*m/dt**2
#             pd = m*a0 + k0*s0 +c0*v0 +m*(6*s0/dt**2 +6*v0/dt+2*a0) + c0*(3*v0/dt +2*v0+ a0*dt/2)
             dpd = m*((6*ds/dt+3*dv)+c0*(3*ds+dt*dv/2))/kd 
             ds = dpd/kd 
             dv = 3*ds/dt -3*v0 -dt*a/2
             Spd[0] = s0+ds
             Spv[0] = v0+dv
             Spa[0] = a
             res[0] = abs(s0+ds-p0[1])/abs(p0[1])+abs(v0+dv-(p0[1]-p0[0])/dt)/abs((p0[1]-p0[0])/dt)            
             for i in range(1, len(ti)):
                 fD += c0*dv
                 fS += k0*ds
                 kd += 3*c0/dt +6*m/dt**2
                 dpd += m*(6*s0/dt**2 +6*v0/dt+2*a0) + c0*(3*v0/dt +2*v0+ a0*dt/2)
                 a = (p0[i] - fD - fS)/m 
                 ds = dpd/kd 
                 dv = 3*ds/dt -3*v0 -dt*a/2
                 Spd[i] = p0[i-1]+ds
                 Spv[i] = Spv[i-1]+dv
                 Spa[i] = a
                 res[i] = abs(p0[i-1]+ds-p0[i])/abs(p0[i])+abs(Spv[i-1]+dv-(p0[i]-p0[i-1])/dt)/abs((p0[i]-p0[i-1])/dt)

        ksei = np.linspace(0, 1, 100)

        temp = 1000
        spd = np.repeat(0.0, winL)
        spv= np.repeat(0.0, winL)
        spa = np.repeat(0.0, winL)
        for count in range(1, len(ksei)):
            kse = ksei[len(ksei)-count]
            print(count)
            sigma1 = np.std(p0)
#            if mode == 0:
#                wi = -np.arctan(B_.imag/(_B.real+0.00000001))
            if mode == 1:
                wi = -np.arccos(Spd/max(abs(Spd)))
                wi = np.nan_to_num(wi)+0.00000001
                FD = fD
                FS = fS
                KD = kd
                DS = ds
                DV = dv
                DPD = dpd 
                Sum_Spd = Spd
                Sum_Spa = Spa
                Sum_Spv = Spv
            wd = wi*(1-kse**2)**0.5
#            print(sigma1)
                #amplitude
#            Ki = wii**2/rho[ti]
#            Di = rho[ti]/Ki
 #           vi = v0*wii**2*Di
#            if mode == 0:
#                a = a0*gamma
#                alpha = wi*dt/vgi
#                b_ = np.fft.ifft(B)/alpha 
#
##               Spv = -kse*b*np.sin(wi*(t-dt))*dt*np.exp(-kse*(wi*(t-dt)))+b*np.cos(wi*(t-dt))*dt*np.exp(-kse*(wi*(t-dt)))
##                Sav = -wi*(2*kse**2-1)*b*np.sin(wi*(t-dt))*dt*np.exp(-kse*(wi*(t-dt)))-2*wi*kse*b*np.cos(wi*(t-dt))*dt*np.exp(-kse*(wi*(t-dt)))
#                #pseudo
#                spv = np.repeat(0.0, winL)
#                for n in range(len(ti)):
#                    pv = 0
#                    for w in wi:
#                        tau = np.linspace(0, ti[n], (ti[n]+1)/dt)
#                        if pv < abs(sum(np.sin(w*(ti[n]-tau))*dt*np.exp(-kse*(w*(ti[n]-dt))))):
#                            pv = sum(np.sin(w*(ti[n]-tau))*dt*np.exp(-kse*(w*(ti[n]-dt))))
#                        spv[n] = np.max(pv*b_).real
##            print(['Spv:', Spv])
#                ki = wi/spv
#                wd = wi*(1-kse**2)**0.5
#                multi = wd**2+(2*wd*kse)**2/2*(np.cos(2*wd*ti))**2+kse*wi*wd*np.cos(2*wd*ti)
#                multi *= 1-np.exp(-2*wi*kse*ti)//wd**2
#                B_cor = abs(np.pi*wi*B**0.5/(2*kse*ki**2+0.000001)*multi)**0.5
#                #            B_cor = B_cor**0.5
#                B_cor[B_cor*p0 < 0] = -B_cor[B_cor*p0 < 0]
#                if diff > sum(abs(B_cor - p0)):
#                    B_ = B_cor
#                    diff = sum(abs(B_cor - p0))
#                    Spd = spv/(wi+0.0000000001)
#                    Spa = wi*spv 
#                    Spv = spv
#                    Sum_B_cor += B_cor
#                    Sum_Spd += spv/wi
#                    Sum_Spa += wi*spv
#                    Sum_Spv += spv
#                if diff > sum(abs(Sum_B_cor/(count+1) - p0)):
#                    B_ = Sum_B_cor/(count+1)
#                    diff = sum(abs(Sum_B_cor/(count+1) - p0))
#                    Spd = Sum_Spd/(count+1)
#                    Spa = Sum_Spa/(count+1) 
#                    Spv = Sum_Spv/(count+1)
#
#                Ampi = np.fft.fft(np.fft.ifft(B_))
##            Ak = np.exp(-np.mean(abs(Ampi - p0)**2)/(2*sigma1**2+0.00001))
#                Ak = np.mean(abs(Ampi - p0)/(sigma1**2+0.00001))
#                print(Ak)
#            #frequency
#        #            Af = 0.3*abs(Amax)
#        #            f0i = np.zeros((winL))
#        #            if ti[-1] < len(freq):
#        #                f0i = freq[ti]
#        #            elif ti[0] <= len(freq):
#        #                f0i[0:len(freq)-ti[0]] = np.array(freq[ti[0]:len(freq)])
#        #                f0i[len(freq)-ti[0]:int(ti[-1]-len(freq)+2)] = np.array(freq[0: ti[-1]-len(freq)])
#        #            else: break
#        #            f0 = f0i[pf >= Af]
#        #            f0mu = np.mean(f0)
#        #            fki = f0i[ pti >= Af]
##                if sum(Di*vi!=0) > 0:
##                    fki = (v0/(Di*vi))**0.5       
##                else:
##                    fki = Di*vi
##                fk = sum((fki - wii)**2) 
#                phi = -np.arctan(Ampi.imag/(Ampi.real+0.0000001))
#                sigma2 = np.std(wi)
#                fki = abs(phi - wi)
##            fk = np.exp(-np.mean(fki)/(2*sigma2**2+0.00001))
#                fk = np.mean(fki/(sigma2+0.00001))               
##            print(fk)
#
##            pf = acc/count
##            print(pf)
##            print(abs(Ak + fk))
#                if abs(Ak + fk) < temp:
#                    temp = abs(Ak + fk)
#                    print(['temp:',temp])
#                if temp <= 1:
#                    print(['range:',range(i*winL,i*winL+winL)])
#                    ACC[i*winL:i*winL+winL] = Spa
#                    Alpha[i*winL:i*winL+winL] = alpha
##                Tau[i*winL:i*winL+winL] = tau/count
#                    Ve[i*winL:i*winL+winL] = Spv
#                    Di[i*winL:i*winL+winL] = Spd
#                    Amp[i*winL:i*winL+winL] = spd
#                    print(['Spd:',Spd, 'D', Di])
#                    print(['Spa:',Spa, 'Accelerte', ACC])
#                    print(['alpha:',alpha, 'Alpha', Alpha])
#                    break
            if mode == 1:
                for i in range(winL):
    #                 alpha = Spd/dt
    #                 A0 = v0/wi**2 - 2*kse*p0/dt/wi**3
    #                 A1 = Spd/dt/wi**2
    #                 A2 = v0 - A0
    #                 A3 = 1/wd*(w0+kse*wi*A2-alpha/wi**2)
    #                 spd = A0 + A1*ti + A2*np.exp(-kse*wi*ti)*np.cos(wd*ti) + A3*np.exp(-kse*wi*ti)*np.sin(wd*ti)
    #                 spv = A1 +(wd*A3 - kse*wi*A2)*np.exp(-kse*wi*ti)*np.cos(wd*ti) - (wd*A2 + kse*wi*A3)*np.exp(-kse*wi*ti)*np.sin(wd*ti) 
    #                 spa = A1 +(wd*A3 - kse*wi*A2)*(-kse*wi*np.exp(-kse*wi*ti)*np.cos(wd*ti)-wd*ti*np.sin(wd*ti)*np.cos(wd*ti)) - (wd*A2 + kse*wi*A3)*(-kse*wi*np.exp(-kse*wi*ti)*np.cos(wd*ti)-wd*ti*np.sin(wd*ti)*np.cos(wd*ti))
                    if i ==0:
                        ds = DS
                        dv = DV
                    FD -= c0*kse*ds
                    FS -= k0*dv
                    KD -= 3*c0*kse/dt +6*m/dt**2
#                    DPD -= (1-(count-1))*m*(6*s0/dt**2 +6*kse*v0/dt+2*a0) + c0*kse*(3*v0/dt +2*v0+ a0*dt/2)
                    DPD -= m*(6*s0/dt**2 +6*v0/dt+2*a0) + c0*kse*(3*v0/dt +2*v0+ a0*dt/2)
                    a = kse*(p0[i] - FD - FS)/m 
                    ds = DPD/KD 
                    dv = kse*ds/dt -3*v0 -dt*a/2
#                    dv = 3*ds/dt -3*v0 -dt*a0/2
            
                    spd[i] = p0[i-1]+ds
                    spv[i] = Spv[i-1]+dv
                    spa[i] = a
                    if res[i] > abs(p0[i-1]+DS-p0[i])/abs(p0[i])+abs(spv[i-1]+dv-(p0[i]-p0[i-1])/dt)/abs((p0[i]-p0[i-1])/dt):
                        res[i] = abs(p0[i-1]+DS-p0[i])/abs(p0[i])+abs(spv[i-1]+dv-(p0[i]-p0[i-1])/dt)/abs((p0[i]-p0[i-1])/dt)
                        Sum_Spd[i] += spd[i]
                        Sum_Spa[i] += spa[i]
                        Sum_Spv[i] += spv[i]
                    if res[i] > abs(Sum_Spd[i]/(count+1) - p0[i]):
                        res[i] = abs(Sum_Spd[i]/(count+1) - p0[i])
                        Spd[i] = Sum_Spd[i]/(count+1)
                        Spa[i] = Sum_Spa[i]/(count+1) 
                        Spv[i] = Sum_Spv[i]/(count+1)
                        Sum_Spd[i] += Spd[i]-spd[i]
                        Sum_Spa[i] += Spa[i]-spa[i]
                        Sum_Spv[i] += Spv[i]-spv[i]
                        
            Ampi = Spd
    #            Ak = np.exp(-np.mean(abs(Ampi - p0)**2)/(2*sigma1**2+0.00001))
            Ak = np.mean(abs(Ampi - p0)/(sigma1**2+0.00001))
            print(Ak)
                #frequency
            #            Af = 0.3*abs(Amax)
            #            f0i = np.zeros((winL))
            #            if ti[-1] < len(freq):
            #                f0i = freq[ti]
            #            elif ti[0] <= len(freq):
            #                f0i[0:len(freq)-ti[0]] = np.array(freq[ti[0]:len(freq)])
            #                f0i[len(freq)-ti[0]:int(ti[-1]-len(freq)+2)] = np.array(freq[0: ti[-1]-len(freq)])
            #            else: break
            #            f0 = f0i[pf >= Af]
            #            f0mu = np.mean(f0)
            #            fki = f0i[ pti >= Af]
    #                if sum(Di*vi!=0) > 0:
    #                    fki = (v0/(Di*vi))**0.5       
    #                else:
    #                    fki = Di*vi
    #                fk = sum((fki - wii)**2) 
            phi = -np.arccos(Ampi/max(abs(Ampi)+0.0000001))
            sigma2 = np.std(wi)
            fki = abs(phi - wi)
    #            fk = np.exp(-np.mean(fki)/(2*sigma2**2+0.00001))
            fk = np.mean(fki/(sigma2+0.00001))
    #            print(fk)
    
    #            pf = acc/count
    #            print(pf)
    #            print(abs(Ak + fk))
            if abs(Ak + fk) < temp:
                temp = abs(Ak + fk)
                print(['temp:',temp])
            if temp <= 1:
                print(['range:',range(tau*winL,tau*winL+winL)])
                ACC[tau*winL:tau*winL+winL] = Spa
#                Alpha[tau*winL:tau*winL+winL] = alpha
    #                Tau[i*winL:i*winL+winL] = tau/count
                Ve[tau*winL:tau*winL+winL] = Spv
                Di[tau*winL:tau*winL+winL] = Spd
                Amp[tau*winL:tau*winL+winL] = Spd
                print(['Spd:',Spd, 'D', Di])
                print(['Spa:',Spa, 'Accelerte', ACC])
#                print(['alpha:',alpha, 'Alpha', Alpha])
                break
            
#    if mode == 0:
#        plt.figure(),
#        plt.subplot(311)
#        plt.plot(Amp[0:len(px)])
#        plt.ylabel('simulated spectrum density')
#        plt.ylim(min(px)-10, max(px)+10) 
#        plt.title(['residule:', np.mean(Amp[0:len(px)] - px)]) 
#        plt.subplot(312)
#        plt.plot(px)
#        plt.ylabel('original wave')
#        plt.subplot(313)
#        plt.plot(Amp[0:len(px)] - px)
#        plt.ylabel('residule')
#        #plt.ylabel('residule')
    if mode == 1:
        plt.figure(),
        plt.subplot(311)
        plt.plot(Amp[0:len(xx)])
        plt.ylabel('simulated spectrum density')
        plt.ylim(min(xx)-3*np.std(xx), max(xx)+3*np.std(xx)) 
        plt.title(['residule:', np.mean(Amp[0:len(xx)] - xx)]) 
        plt.subplot(312)
        plt.plot(xx)
        plt.ylabel('original wave')
        plt.subplot(313)
        plt.plot(Amp[0:len(xx)] - xx)
        plt.ylabel('residule')
    
    plt.figure(),    
    plt.subplot(311)
    plt.plot(ACC)
    plt.ylabel('accelerate of spectrum')

        #        plt.ylim(min(px)-10, max(px)+10) 
    plt.title('spectrum deterministic quntities:') 
    plt.subplot(312)
    plt.plot(Ve)
    plt.ylabel('velocity')

    plt.subplot(313)
    plt.plot(Di)
    plt.ylabel('displacement')

        
    return Amp[0:len(px)], Amp[0:len(px)] -  px, ACC, Ve, Di         
                
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

def v_rho_r(ri, r, u, du, v, dv, w, vp_, vs_,a,l):
    data = r[:, ri]
    Pv = []
    Sv = []
    RR = []
    EE = []
    W2 = []
    r0 = 0
    for i in range(np.shape(data)[1]):
        rr = data[:, i]
        rr = rr[rr != 0,]
#        rho = rho_[r0:r0+len(rr)]
        print(rr)
#        r0 = len(rr)
#        C_ = du[rr!=0,]*du[rr!=0,]
        C_ = du[r0:r0+len(rr),]*du[r0:r0+len(rr),]
        A_ = np.repeat(rr**(-2), np.shape(u)[1]).reshape(len(rr), np.shape(u)[1])*2*du[r0:r0+len(rr),] - l*(l+1)*(dv[r0:r0+len(rr),]+u[r0:r0+len(rr),]-v[r0:r0+len(rr),])/(v[r0:r0+len(rr),]+0.00001)**2
        F_ = np.repeat(2*rr**(-1),np.shape(u)[1]).reshape(len(rr), np.shape(u)[1])*du[r0:r0+len(rr),]*(2*u[r0:r0+len(rr),]-l*(l+1)*v[r0:r0+len(rr),])
        L_ = l*(l+1)*np.repeat(rr**(-1), np.shape(u)[1]).reshape(len(rr), np.shape(u)[1])*(dv[r0:r0+len(rr),] + u[r0:r0+len(rr),] - v[r0:r0+len(rr),])**2
        N_ = np.repeat(rr**(-2),np.shape(u)[1]).reshape(len(rr), np.shape(u)[1])*(l+2)*(l+1)*l*(l-1)*v[r0:r0+len(rr),]**2-(2*u[r0:r0+len(rr),]-l*(l+1)*v[r0:r0+len(rr),])**2
        eta = F_/(A_ - 2*L_+0.000001)
        
#        C_ = C_[C_!= 0]
#        A_ = A_[A_!= 0]
#        F_ = F_[F_!= 0]
#        L_ = L_[L_!= 0]
#        N_ = N_[N_!= 0]
#        eta = eta[eta != 0]
        dA_ = A_
        dA_[1:-1] = A_[1:-1]-A_[0:-2]
        dC_ = C_
        dC_[1:-1] = C_[1:-1]-C_[0:-2]
        dF_ = F_
        dF_[1:-1] = F_[1:-1]-F_[0:-2]
        dL_ = L_
        dL_[1:-1] = L_[1:-1]-L_[0:-2]
        dN_ = N_
        dN_[1:-1] = N_[1:-1]-N_[0:-2]
        if i == 12:
            print(i)
            Vpv1 = 11.2622 - 6.3640*rr/a+5.5242*rr**2/a**2
            Vph1 = 11.2622 - 6.3640*rr/a+5.5242*rr**2/a**2
            Vsv1 = 3.6678 - 4.4475*rr**2/a**2
            Vsh1 = 3.6678 - 4.4475*rr**2/a**2
            rho = 13.0885-8.8381**rr**2/a**2
        if i == 11:
            print(i)
            Vpv1 = 12.5815 - 1.2638*rr/a-3.6426*rr**2/a**2-5.5281*rr**3/a**3
            Vph1 = 12.5815 - 1.2638*rr/a-3.6426*rr**2/a**2-5.5281*rr**3/a**3
            Vsv1 = np.repeat(0, len(rr))
            Vsh1 = np.repeat(0, len(rr))
            rho = 12.5815-1.2638*rr/a-3.6426*rr**2/a**2-5.5281*rr**3/a**3        
        if i == 10:
            print(i)
            Vpv1 = 15.3891 - 5.3181*rr/a+5.5242*rr**2/a**2-2.5514*rr**3/a**3
            Vph1 = 15.3891 - 5.3181*rr/a+5.5242*rr**2/a**2-2.5514*rr**3/a**3
            Vsv1 = 6.9254 + 1.4672*rr/a-2.0834*rr**2/a**2+0.9783*rr**3/a**3
            Vsh1 = 6.9254 + 1.4672*rr/a-2.0834*rr**2/a**2+0.9783*rr**3/a**3
            rho = 7.9656-6.4761*rr/a+5.5283*rr**2/a**2-3.0807*rr**3/a**3
        if i == 9:
            print(i)
            Vpv1 = 24.9520 - 40.4673*rr/a+51.4832*rr**2/a**2-26.6419*rr**3/a**3
            Vph1 = 24.9520 - 40.4673*rr/a+51.4832*rr**2/a**2-26.6419*rr**3/a**3
            Vsv1 = 11.1671 - 13.7818*rr/a+17.4575*rr**2/a**2-9.2777*rr**3/a**3
            Vsh1 = 11.1671 - 13.7818*rr/a+17.4575*rr**2/a**2-9.2777*rr**3/a**3
            rho = 7.9656-6.4761*rr/a+5.5283*rr**2/a**2-3.0807*rr**3/a**3
        if i == 8:
            print(i)
            Vpv1 = 29.2766 - 23.6027*rr/a+5.5242*rr**2/a**2-2.5514*rr**3/a**3
            Vph1 = 29.2766 - 23.6027*rr/a+5.5242*rr**2/a**2-2.5514*rr**3/a**3
            Vsv1 = 22.3459 - 17.2473*rr/a-2.0834*rr**2/a**2+0.9783*rr**3/a**3
            Vsh1 = 22.3459 - 17.2473*rr/a-2.0834*rr**2/a**2+0.9783*rr**3/a**3
            rho = 7.9656-6.4761*rr/a+5.5283*rr**2/a**2-3.0807*rr**3/a**3
        if i == 7:
            print(i)
            Vpv1 = 19.0957 - 9.8672*rr/a
            Vph1 = 19.0957 - 9.8672*rr/a
            Vsv1 = 9.9839 - 4.9324*rr/a
            Vsh1 = 9.9839 - 4.9324*rr/a
            rho = 5.3197-1.4836*rr/a
        if i == 6:
            print(i)
            Vpv1 = 39.7027 - 32.6166*rr/a
            Vph1 = 39.7027 - 32.6166*rr/a
            Vsv1 = 22.3512 - 18.5856*rr/a
            Vsh1 = 22.3512 - 18.5856*rr/a
            rho = 11.2494-8.0298*rr/a
        if i == 5:
            print(i)
            Vpv1 = 20.3926 + 8.9496*rr/a
            Vph1 = 20.3926 + 8.9496*rr/a
            Vsv1 = 8.9496 - 4.4597*rr/a
            Vsh1 = 8.9496 - 4.4597*rr/a
            rho = 7.1089+3.8045*rr/a
        if i in [3, 4]:# with symmetry axis
            print(i)
            Vpv1 = 4.1875 + 3.9382*rr/a
            Vph1 = 4.1875 + 3.9382*rr/a
            Vsv1 = 2.1519 + 2.3481*rr/a
            Vsh1 = 2.1519 + 2.3481*rr/a
            rho = 2.691+0.6924*rr/a
        if i in [1, 2]:
            print(i)
            Vpv1 = 4.1875 + 3.9382*rr/a
            Vph1 = 4.1875 + 3.9382*rr/a
            Vsv1 = 2.1519 + 2.3481*rr/a
            Vsh1 = 2.1519 + 2.3481*rr/a
            if i == 1:
                rho = np.repeat(2.6, len(rr))
            else:
                rho = np.repeat(2.9, len(rr))
        if i == 0:
            Vpv1 = np.repeat(1.45, len(rr))
            Vsv1 = np.repeat(0, len(rr))
            Vph1 = np.repeat(1.45, len(rr))
            Vsh1 = np.repeat(0, len(rr))
            rho = np.repeat(1.020, len(rr))
#        rho = np.repeat(rho, np.shape(u)[1]).reshape(len(rr), np.shape(u)[1])
        r0 = len(rr)
        print(rho)
        print(rr)
        print(Vph1)
        A = rho*Vph1**2
        C = rho*Vpv1**2 
        N = rho*Vsh1**2
        L = rho*Vsv1**2 
        
#        rho_ = np.repeat(rho, np.shape(u)[1]).reshape(len(rr), np.shape(u)[1])
        drho_ = rho/10**15
        drho_[1:-1] = rho[1:-1] - rho[0:-2]
        drho_ = np.repeat(drho_, np.shape(u)[1]).reshape(len(rr), np.shape(u)[1])
#        A = np.repeat(A, np.shape(u)[1]).reshape(len(rr), np.shape(u)[1])        
#        C = np.repeat(C, np.shape(u)[1]).reshape(len(rr), np.shape(u)[1])
#        N = np.repeat(N, np.shape(u)[1]).reshape(len(rr), np.shape(u)[1])        
#        L = np.repeat(L, np.shape(u)[1]).reshape(len(rr), np.shape(u)[1])
        print(np.shape(A))
        print(A)
        A = np.repeat(A, np.shape(u)[1]).reshape(len(rr), np.shape(u)[1])        
        C = np.repeat(C, np.shape(u)[1]).reshape(len(rr), np.shape(u)[1])
        N = np.repeat(N, np.shape(u)[1]).reshape(len(rr), np.shape(u)[1])        
        L = np.repeat(L, np.shape(u)[1]).reshape(len(rr), np.shape(u)[1])
        F = eta*(A-2*L)
        if sum(sum(A!=0)) != 0:
            A = A[A != 0]
        if sum(sum(C!=0)) != 0:
            C = C[C != 0]
        if sum(sum(N!=0)) != 0:
            N = N[N != 0]
        if sum(sum(L!=0)) != 0:
            L = L[L != 0]
        dA = A
        dA[1:-1] = A[1:-1]-A[0:-2]
        dC = C
        dC[1:-1] = C[1:-1]-C[0:-2]
        dN = N
        dN[1:-1] = N[1:-1]-N[0:-2]
        dL = L
        dL[1:-1] = L[1:-1]-L[0:-2]
        dF = F
        dF[1:-1] = F[1:-1]-F[0:-2]
        
        R_ = -0.5*np.repeat(rr**2, np.shape(u)[1]).reshape(len(rr), np.shape(u)[1])* np.repeat( a + np.dot(Vpv1**2, C_) + np.dot(Vph1**2,(A_ + eta*F_)) + np.dot(Vsv1**2,(L_ - 2*eta*F_)) + np.dot(Vsh1**2, N_), len(rr)).reshape(len(rr), np.shape(u)[1])
        pv = -np.repeat(rr**2*rho, np.shape(u)[1]).reshape(len(rr), np.shape(u)[1])* np.repeat(np.dot(Vpv1, C_), len(rr)).reshape(len(rr), np.shape(v)[1])
        ph = -np.repeat(rr**2*rho,  np.shape(u)[1]).reshape(len(rr), np.shape(u)[1])* np.repeat(np.dot(Vph1, A_ + eta*F_), len(rr)).reshape(len(rr), np.shape(v)[1])
        sv = -np.repeat(rr**2*rho, np.shape(u)[1]).reshape(len(rr), np.shape(u)[1])*np.repeat(np.dot(Vsv1,(L_ - 2*eta*F_)), len(rr)).reshape(len(rr), np.shape(v)[1])
        sh = -np.repeat(rr**2*rho, np.shape(u)[1]).reshape(len(rr), np.shape(u)[1])*np.repeat(np.dot(Vsh1, N_),len(rr)).reshape(len(rr), np.shape(v)[1])
        pv[pv==0] += vp_[0]*10**21
        ph[ph==0] += vp_[0]*10**21
        sv[sv==0] += vs_[0]*10**21
        sh[sh==0] += vs_[0]*10**21
        E_ = -np.repeat(0.5*rr**2*rho*(Vph1**2 - 2*Vsv1**2), np.shape(u)[1]).reshape(len(rr), np.shape(u)[1])* F_
#        R_ = R_[R_ != 0]
#        pv = pv[pv != 0]
#        ph = ph[ph != 0]
#        sv = sv[sv != 0]
#        sh = sh[sh != 0]
#        E_ = E_[E_ != 0]
    
        dw2w2 = np.repeat(rr**3, np.shape(u)[1]).reshape(len(rr), np.shape(u)[1])*(dA_*A.reshape(len(rr), np.shape(u)[1])+dC.reshape(len(rr), np.shape(u)[1])*C_+dF.reshape(len(rr), np.shape(u)[1])*F_+dL.reshape(len(rr), np.shape(u)[1])*L_+dN.reshape(len(rr), np.shape(u)[1])*N_+drho_*R_.reshape(len(rr), np.shape(u)[1]))
        
        Pv.append(np.nan_to_num((pv**2+ph**2)**(0.5)))
        Sv.append(np.nan_to_num((sv**2+sh**2)**(0.5)))
        RR.append(R_)
        EE.append(E_)
        W2.append(dw2w2)   
        if i == np.shape(data)[1]-1:
            print(i)
            break
    return Pv, Sv, W2, RR, EE  


#rupture data
#iport obspy
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
T = range(int(p*fs),int(s*fs+3*(s-p)*fs))
#winL = int(len(T)/24)
winL = int(fs)
#PT,v,k,Vg,x,a,TP,Phii,varg=P(Pwave, rho, winL, len(T))    m
PT, xp, vp, kp, Vgp, vavgp, TP  = PS(dataz[T], winL, len(T),'P')

T = range(int(s*fs+1), int(s*fs+3*(s-p)*fs))
ref = Swave - PT[(s-p)*fs+1: int(s*fs+3*(s-p)*fs)]
MCMC_1D_improved(dataz[T], 6, T, rho_, fs)


#freq = np.fft.fft(ref) 
#w = np.arctan(np.real(freq) / (np.imag(freq)+0.0001))
#ST, xs, vs, ks, vgs, vavgs, TS = PS(ref, rho, winL, len(T),'S')
#D = max(DP,DS)
Tbody = int((s-p)*fs) + TS
Sf(datasete[int(s*fs+3*(s-p)*fs+1):int(s*fs+(3*(s-p)+6)*fs+1)], datasetn[int(s*fs+3*(s-p)*fs+1):int(s*fs+(3*(s-p)+6)*fs+1)], datasetz[int(s*fs+3*(s-p)*fs+1):int(s*fs+(3*(s-p)+6)*fs+1)], rho, 5, 3*fs, 3*fs, 3*fs, 25)


refx = datasete.filter('highpass', freq=1.6, corners=5, zerophase=True) 
freqx = np.fft.fft(refx) 
wx = np.arctan(np.real(freqx) / (np.imag(freqx)+0.0001))
refy = datasete.filter('highpass', freq=1.6, corners=5, zerophase=True) 
freqy = np.fft.fft(refy) 
wy = np.arctan(np.real(freqy) / (np.imag(freqy)+0.0001))
refz = datasete.filter('highpass', freq=1.6, corners=5, zerophase=True) 
freqz = np.fft.fft(refz) 
wz = np.arctan(np.real(freqz) / (np.imag(freqz)+0.0001))

######################normal mode
a = 6371
l = 11

r = np.zeros((7504, 13))
u = np.zeros((len(datasetz), 3))
u[:,0] = datasetz
v = np.zeros((len(datasetz), 3))
v[:,1] = datasetn
w = np.zeros((len(datasetz), 3))
w[:,2] = datasete
#url: ds.iris.edu/ds/products/emc-earthmodels
#u = sp.Symbol('u')
#v = sp.Symbol('v')

du= u
dv= v
dw= w

du[1:-1] = u[1:-1]-u[0:-2]
dv[1:-1] = v[1:-1]-v[0:-2]
dw[1:-1] = w[1:-1]-w[0:-2]

rin = np.linspace(0, 1221.5, 1222)
r[range(0, 1222), 0] = rin
rout = np.linspace(1221.5, 3480, 3480-1222)
r[range(1222-1222, 3480-1222), 1] = rout
rlow1 = np.linspace(3480, 3630, 3630-3480)
r[range(3480-3480, 3630-3480), 2] = rlow1
rlow2 = np.linspace(3630, 5600, 5600-3630)
r[range(3630-3630, 5600-3630), 3] = rlow2
rlow3 = np.linspace(5600, 5701, 5701-5600)
r[range(5600-5600, 5701-5600), 4] = rlow3 
rtr1 = np.linspace(5701, 5771, 5771-5701)
r[range(5701-5701, 5771-5701), 5] = rtr1
rtr2 = np.linspace(5771, 5971, 5971-5771)
r[range(5771-5771, 5971-5771), 6] = rtr2
rtr3 = np.linspace(5971, 6151, 6151-5971)
r[range(5971-5971, 6151-5971), 7] = rtr3
rlvz = np.linspace(6151, 6291, 6291-6151)
r[range(6151-6151, 6291-6151), 8] = rlvz 
rlid = np.linspace(6292, 6347, 6347-6292)
r[range(6292-6292, 6347-6292), 9] = rlid 
rcrust1 = np.linspace(6347, 6356, 6356-6347)
r[range(6347-6347, 6356-6347), 10] = rcrust1 
rcrust2 = np.linspace(6356, 6368, 6368-6356)
r[range(6356-6356, 6368-6356), 11] = rcrust2 
rocean = np.linspace(6368, 6371, 6371-6368)
r[range(6368-6368, 6371-6368), 12] = rocean


ri = np.array(np.linspace(12, 8, 5), dtype = 'int')# 9.LVZ, 10.LID, 11.Crust, 12.Ocean 
ri_ = np.array(np.linspace(12, 0, 13), dtype = 'int')# 9.LVZ, 10.LID, 11.Crust, 12.Ocean 

SAW642ANB = pd.read_csv('D:/TUT/Medical/biophysics/experiment/bayes/Bayes/Bayes Inference/samples/SAW642ANB_kmps.csv', sep = '|',header =1, skiprows = 67)
data = np.array(SAW642ANB)
rho_ = data[data[:,1]==0, 5]
depth_ = data[data[:,1]==0, 2]
vs_= data[data[:,1]==0, 3]
vp_= data[data[:,1]==0, 4]
Pv, Sv, dw2w2, RR, EE = v_rho_r(ri_, r, u, du, v, dv, w, vp_, vs_, a, l)
plt.figure,
plt.subplot(2,2,1)
plt.plot(depth_[depth_<=(len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4]))], vp_[depth_<= (len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4]))], 'r')
plt.subplot(2,2,2)
plt.plot(depth_[depth_<=(len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])+len(Sv[4]))], vs_[depth_<= (len(Sv[0])+len(Pv[1])+len(Sv[2])+len(Sv[3])+len(Pv[4]))], 'k')
plt.subplot(2,2,3)
plt.plot(range(len(Pv[0])), Pv[0], 'r-', range(len(Pv[0]),len(Pv[0])+len(Pv[1])), Pv[1], 'r-', range(len(Pv[0])+len(Pv[1]), len(Pv[0])+len(Pv[1])+len(Pv[2])), Pv[2], 'r-', range( len(Pv[0])+len(Pv[1])+len(Pv[2]),  len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])), Pv[3], 'r-', range( len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3]), len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])), Pv[4], 'r-')
plt.subplot(2,2,4)
plt.plot(range(len(Sv[0])), Sv[0], 'k-', range(len(Sv[0]),len(Sv[0])+len(Sv[1])), Sv[1], 'k-', range(len(Sv[0])+len(Sv[1]), len(Sv[0])+len(Sv[1])+len(Sv[2])), Sv[2], 'k-', range(len(Sv[0])+len(Sv[1])+len(Sv[2]), len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])), Sv[3], 'k-', range(len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3]), len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])+len(Sv[4])), Sv[4], 'k-')

plt.figure,
plt.subplot(2,2,1)
plt.plot(depth_[depth_<=(len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])+len(Pv[7])+len(Pv[8])+len(Pv[9])+len(Pv[10])+len(Pv[11])+len(Pv[12]))], vp_[depth_<= (len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])+len(Pv[7])+len(Pv[8])+len(Pv[9])+len(Pv[10])+len(Pv[11])+len(Pv[12]))], 'r')
plt.subplot(2,2,2)
plt.plot(depth_[depth_<=(len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])+len(Sv[4])+len(Sv[5])+len(Sv[6])+len(Sv[7])+len(Sv[8])+len(Sv[9])+len(Sv[10])+len(Sv[11])+len(Sv[12]))], vs_[depth_<= (len(Sv[0])+len(Pv[1])+len(Sv[2])+len(Sv[3])+len(Sv[4])+len(Sv[5])+len(Sv[6])+len(Sv[7])+len(Sv[8])+len(Sv[9])+len(Sv[10])+len(Sv[11])+len(Sv[12]))], 'k')
plt.subplot(2,2,3)
plt.plot(range(len(Pv[0])), Pv[0][:,0], 'ro', range(len(Pv[0]),len(Pv[0])+len(Pv[1])), Pv[1][:,0], 'ro', range(len(Pv[0])+len(Pv[1]), len(Pv[0])+len(Pv[1])+len(Pv[2])), Pv[2][:,0], 'ro', range( len(Pv[0])+len(Pv[1])+len(Pv[2]),  len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])), Pv[3][:,0], 'ro', range( len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3]), len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])), Pv[4][:,0], 'ro',  range(len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4]), len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])), Pv[5][:,0], 'ro',  range(len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5]), len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])), Pv[6][:,0], 'ro', range(len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6]), len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])+len(Pv[7])), Pv[7][:,0], 'ro', range(len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])+len(Pv[7]), len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])+len(Pv[7])+len(Pv[8])), Pv[8][:,0], 'ro', range(len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])+len(Pv[7])+len(Pv[8]), len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])+len(Pv[7])+len(Pv[8])+len(Pv[9])), Pv[9][:,0], 'ro', range(len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])+len(Pv[7])+len(Pv[8])+len(Pv[9]), len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])+len(Pv[7])+len(Pv[8])+len(Pv[9])+len(Pv[10])), Pv[10][:,0], 'ro', range(len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])+len(Pv[7])+len(Pv[8])+len(Pv[9])+len(Pv[10]), len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])+len(Pv[7])+len(Pv[8])+len(Pv[9])+len(Pv[10])+len(Pv[11])), Pv[11][:,0], 'ro', range(len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])+len(Pv[7])+len(Pv[8])+len(Pv[9])+len(Pv[10])+len(Pv[11]), len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])+len(Pv[7])+len(Pv[8])+len(Pv[9])+len(Pv[10])+len(Pv[11])+len(Pv[12])), Pv[12][:,0], 'ro')
plt.subplot(2,2,4)
plt.plot(range(len(Sv[0])), Sv[0][:,0], 'ko', range(len(Sv[0]),len(Sv[0])+len(Sv[1])), Sv[1][:,0], 'ko', range(len(Sv[0])+len(Sv[1]), len(Sv[0])+len(Sv[1])+len(Sv[2])), Sv[2][:,0], 'ko', range( len(Sv[0])+len(Sv[1])+len(Sv[2]),  len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])), Sv[3][:,0], 'ko', range( len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3]), len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])+len(Sv[4])), Sv[4][:,0], 'ko',  range(len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4]), len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])+len(Sv[4])+len(Sv[5])), Sv[5][:,0], 'ko',  range(len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5]), len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])+len(Sv[4])+len(Sv[5])+len(Sv[6])), Sv[6][:,0], 'ko', range(len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6]), len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])+len(Pv[7])), Sv[7][:,0], 'ko', range(len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])+len(Sv[4])+len(Sv[5])+len(Sv[6])+len(Sv[7]), len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])+len(Pv[7])+len(Pv[8])), Sv[8][:,0], 'ko', range(len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])+len(Sv[4])+len(Sv[5])+len(Sv[6])+len(Sv[7])+len(Sv[8]), len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])+len(Sv[4])+len(Sv[5])+len(Sv[6])+len(Sv[7])+len(Sv[8])+len(Sv[9])), Sv[9][:,0], 'ko', range(len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])+len(Sv[4])+len(Sv[5])+len(Sv[6])+len(Sv[7])+len(Sv[8])+len(Sv[9]), len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])+len(Sv[4])+len(Sv[5])+len(Sv[6])+len(Sv[7])+len(Sv[8])+len(Sv[9])+len(Sv[10])), Sv[10][:,0], 'ko', range(len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])+len(Sv[4])+len(Sv[5])+len(Sv[6])+len(Sv[7])+len(Sv[8])+len(Sv[9])+len(Sv[10]), len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])+len(Sv[4])+len(Sv[5])+len(Sv[6])+len(Sv[7])+len(Sv[8])+len(Sv[9])+len(Sv[10])+len(Sv[11])), Sv[11][:,0], 'ko', range(len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])+len(Sv[4])+len(Sv[5])+len(Sv[6])+len(Sv[7])+len(Sv[8])+len(Sv[9])+len(Sv[10])+len(Sv[11]), len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])+len(Sv[4])+len(Sv[5])+len(Sv[6])+len(Sv[7])+len(Sv[8])+len(Sv[9])+len(Sv[10])+len(Sv[11])+len(Sv[12])), Sv[12][:,0], 'ko')


plt.figure,
plt.subplot(2,2,1)
plt.plot(depth_[depth_<=(len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])+len(Pv[7])+len(Pv[8])+len(Pv[9])+len(Pv[10])+len(Pv[11])+len(Pv[12]))], vp_[depth_<= (len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])+len(Pv[7])+len(Pv[8])+len(Pv[9])+len(Pv[10])+len(Pv[11])+len(Pv[12]))], 'r')
plt.subplot(2,2,2)
plt.plot(depth_[depth_<=(len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])+len(Sv[4])+len(Sv[5])+len(Sv[6])+len(Sv[7])+len(Sv[8])+len(Sv[9])+len(Sv[10])+len(Sv[11])+len(Sv[12]))], vs_[depth_<= (len(Sv[0])+len(Pv[1])+len(Sv[2])+len(Sv[3])+len(Sv[4])+len(Sv[5])+len(Sv[6])+len(Sv[7])+len(Sv[8])+len(Sv[9])+len(Sv[10])+len(Sv[11])+len(Sv[12]))], 'k')
plt.subplot(2,2,3)
plt.plot(range(len(Pv[0])), Pv[0][:,0], 'ro')
plt.plot(range(len(Pv[0]),len(Pv[0])+len(Pv[1])), Pv[1][:,0], 'ro')
plt.plot(range(len(Pv[0])+len(Pv[1]), len(Pv[0])+len(Pv[1])+len(Pv[2])), Pv[2][:,0], 'ro')
plt.plot(range( len(Pv[0])+len(Pv[1])+len(Pv[2]),  len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])), Pv[3][:,0], 'ro')
plt.plot(range( len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3]), len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])), Pv[4][:,0], 'ro')
plt.plot(range(len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4]), len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])), Pv[5][:,0], 'ro')
plt.plot(range(len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5]), len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])), Pv[6][:,0], 'ro')
plt.plot(range(len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6]), len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])+len(Pv[7])), Pv[7][:,0], 'ro')
plt.plot(range(len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])+len(Pv[7]), len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])+len(Pv[7])+len(Pv[8])), Pv[8][:,0], 'ro')
plt.plot(range(len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])+len(Pv[7])+len(Pv[8]), len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])+len(Pv[7])+len(Pv[8])+len(Pv[9])), Pv[9][:,0], 'ro')
plt.plot(range(len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])+len(Pv[7])+len(Pv[8])+len(Pv[9]), len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])+len(Pv[7])+len(Pv[8])+len(Pv[9])+len(Pv[10])), Pv[10][:,0], 'ro')
plt.plot(range(len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])+len(Pv[7])+len(Pv[8])+len(Pv[9])+len(Pv[10]), len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])+len(Pv[7])+len(Pv[8])+len(Pv[9])+len(Pv[10])+len(Pv[11])), Pv[11][:,0], 'ro')
#plt.plot(range(len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])+len(Pv[7])+len(Pv[8])+len(Pv[9])+len(Pv[10])+len(Pv[11]), len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])+len(Pv[7])+len(Pv[8])+len(Pv[9])+len(Pv[10])+len(Pv[11])+len(Pv[12])),  Pv[12][:,0], 'ro')                                                                                                                                                                                                                                                                                                                         , 'ro')
plt.subplot(2,2,4)
plt.plot(range(len(Sv[0])), Sv[0][:,0], 'ko')
plt.plot(range(len(Sv[0]),len(Sv[0])+len(Sv[1])), Sv[1][:,0], 'ko')
plt.plot(range(len(Sv[0])+len(Sv[1]), len(Sv[0])+len(Sv[1])+len(Sv[2])), Sv[2][:,0], 'ko')
plt.plot(range( len(Sv[0])+len(Sv[1])+len(Sv[2]),  len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])), Sv[3][:,0], 'ko')
plt.plot(range( len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3]), len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])+len(Sv[4])), Sv[4][:,0], 'ko')
plt.plot(range(len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4]), len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])+len(Sv[4])+len(Sv[5])), Sv[5][:,0], 'ko')
plt.plot(range(len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5]), len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])+len(Sv[4])+len(Sv[5])+len(Sv[6])), Sv[6][:,0], 'ko')
plt.plot(range(len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6]), len(Pv[0])+len(Pv[1])+len(Pv[2])+len(Pv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])+len(Pv[7])), Sv[7][:,0], 'ko')
plt.plot(range(len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])+len(Sv[4])+len(Sv[5])+len(Sv[6])+len(Sv[7]), len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])+len(Pv[4])+len(Pv[5])+len(Pv[6])+len(Pv[7])+len(Pv[8])), Sv[8][:,0], 'ko')
plt.plot(range(len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])+len(Sv[4])+len(Sv[5])+len(Sv[6])+len(Sv[7])+len(Sv[8]), len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])+len(Sv[4])+len(Sv[5])+len(Sv[6])+len(Sv[7])+len(Sv[8])+len(Sv[9])), Sv[9][:,0], 'ko')
plt.plot(range(len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])+len(Sv[4])+len(Sv[5])+len(Sv[6])+len(Sv[7])+len(Sv[8])+len(Sv[9]), len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])+len(Sv[4])+len(Sv[5])+len(Sv[6])+len(Sv[7])+len(Sv[8])+len(Sv[9])+len(Sv[10])), Sv[10][:,0], 'ko')
plt.plot(range(len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])+len(Sv[4])+len(Sv[5])+len(Sv[6])+len(Sv[7])+len(Sv[8])+len(Sv[9])+len(Sv[10]), len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])+len(Sv[4])+len(Sv[5])+len(Sv[6])+len(Sv[7])+len(Sv[8])+len(Sv[9])+len(Sv[10])+len(Sv[11])), Sv[11][:,0], 'ko')
#plt.plot(range(len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])+len(Sv[4])+len(Sv[5])+len(Sv[6])+len(Sv[7])+len(Sv[8])+len(Sv[9])+len(Sv[10])+len(Sv[11]), len(Sv[0])+len(Sv[1])+len(Sv[2])+len(Sv[3])+len(Sv[4])+len(Sv[5])+len(Sv[6])+len(Sv[7])+len(Sv[8])+len(Sv[9])+len(Sv[10])+len(Sv[11])+len(Sv[12])), Sv[12][:,0], 'ko')




#beam forming