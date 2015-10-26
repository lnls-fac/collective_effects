#!/usr/bin/env python3

import numpy as np
import scipy.special as scy

c   = 299792458
mu0 = 4*np.pi*1e-7
ep0 = 1/c**2/mu0
Z0  = mu0*c

def longitudinal_resonator(Rs, Q, wr, w):
    # Modelagem de impedância por soma de ressonadores.
    # Inputs:
    # Rs    = Vetor com as resistências Shunts [Ohm]
    # Q     = Vetor com os fatores de qualidade
    # wr    = Vetor com as frequências angulares de ressonância de cada oscilador [rad/s]
    # w     = Frequências angulares em que sera calculada a impedância [rad/s]
    #
    # Outputs:
    # Zl    = impedância longitudinal (real e imaginária) calculada para w
    #
    # Todos os outputs tem as dimensoes do SI
    #
    # Referência orginal:
    #
    # Chao, A., Physics of Collective Beam Instabilities in High Energy
    # Accelerators, Wiley 1993.
    if len(Rs)>1:
        Rs = Rs[:,None] # I am using broadcasting
        Q  = Q[:,None]
        wr = wr[:,None]
        Zl = w*Rs / (w+1j*Q*(wr - w**2/wr))
        return Zl.sum(0).flatten()

    Zl = w*Rs / (w+1j*Q*(wr - w**2/wr))
    return Zl

def transverse_resonator(Rs, Q, wr, w):
    # Modelagem de impedância por soma de ressonadores.
    # Inputs:
    # Rs    = Vetor com as resist??ncias Shunts [Ohm]
    # Q     = Vetor com os fatores de qualidade
    # wr    = Vetor com as frequ??ncias angulares de resson??ncia de cada oscilador [rad/s]
    # w     = Frequ??ncias angulares em que sera calculada a imped??ncia [rad/s]
    #
    # Outputs:
    # Zt    = imped??ncia transversal (real e imagin??ria) calculada para w
    #
    # Todos os outputs tem as dimens???es do SI
    #
    # Referência orginal:
    #
    # Chao, A., Physics of Collective Beam Instabilities in High Energy
    # Accelerators, Wiley 1993.
    if len(Rs)>1:
        Rs = Rs[:,None] # I am using broadcasting
        Q  = Q[:,None]
        wr = wr[:,None]
        Zt = wr*Rs/(w + 1j*Q*(wr - w**2/wr))
        return Zt.sum(0).flatten()

    Zt = wr*Rs/(w + 1j*Q*(wr - w**2/wr))
    return Zt

def resistive_multilayer_round_pipe(w,epr,mur,b,L,E):

    def Mtil(m, epr, mur, bet, nu, b):
        def produto(A,B):
            C = np.zeros((A.shape[0],B.shape[1],A.shape[2]),dtype=complex)
            for i in range(A.shape[0]):
                for j in range(B.shape[1]):
                    for k in range(A.shape[1]):
                        C[i,j,:] = C[i,j,:] + A[i,k,:]*B[k,j,:]
            return C

        for i in range(len(b)): # lembrando que length(b) = número de camadas - 1
            x = nu[i+1,:] * b[i]
            y = nu[i  ,:] * b[i]
            Mt = np.zeros((4,4,w.shape[0]),dtype=complex)

            if i<len(b)-1:
                D = np.zeros((4,4,nu.shape[1]),dtype=complex)
                z = nu[i+1,:]*b[i+1]
                if not any(z.real<0):
                    ind = (z.real<60)

                    A = scy.iv(m,z[ind])
                    B = scy.kv(m,z[ind])
                    C = scy.iv(m,x[ind])
                    E = scy.kv(m,x[ind])

                    D[0,0,:]    =  1
                    D[1,1,ind]  = - B*C/(A*E)
                    D[1,1,~ind] = - np.exp(-2*(z[~ind]-x[~ind]))
                    D[2,2,:]    =  1
                    D[3,3,ind]  = - B*C/(A*E)
                    D[3,3,~ind] = - np.exp(-2*(z[~ind]-x[~ind]))

            Mt[0,0,:] = -nu[i+1,:]**2*b[i]/epr[i+1,:]*(
                    epr[i+1,:]/nu[i+1,:]*(-scy.kve(m-1,x)/scy.kve(m,x) - m/x)
                    - epr[i,:]/nu[i,:]  *( scy.ive(m-1,y)/scy.ive(m,y) - m/y))
            Mt[0,1,:] = -nu[i+1,:]**2*b[i]/epr[i+1,:]*(
                    epr[i+1,:]/nu[i+1,:]*(-scy.kve(m-1,x)/scy.kve(m,x) - m/x)
                    - epr[i,:]/nu[i,:]  *(-scy.kve(m-1,y)/scy.kve(m,y) - m/y))
            Mt[1,0,:] = -nu[i+1,:]**2*b[i]/epr[i+1,:]*(
                    epr[i+1,:]/nu[i+1,:]*( scy.ive(m-1,x)/scy.ive(m,x) - m/x)
                    - epr[i,:]/nu[i,:]  *( scy.ive(m-1,y)/scy.ive(m,y) - m/y))
            Mt[1,1,:] = -nu[i+1,:]**2*b[i]/epr[i+1,:]*(
                    epr[i+1,:]/nu[i+1,:]*( scy.ive(m-1,x)/scy.ive(m,x) - m/x)
                    - epr[i,:]/nu[i,:]  *(-scy.kve(m-1,y)/scy.kve(m,y) - m/y))

            Mt[0,2,:] = (nu[i+1,:]**2/nu[i,:]**2 - 1)*m/(bet*epr[i+1,:])
            Mt[0,3,:] = Mt[0,2,:]
            Mt[1,2,:] = Mt[0,2,:]
            Mt[1,3,:] = Mt[0,2,:]
            Mt[2,2,:] = -nu[i+1,:]**2*b[i]/mur[i+1,:]*(
                    mur[i+1,:]/nu[i+1,:]*(-scy.kve(m-1,x)/scy.kve(m,x) - m/x)
                    - mur[i,:]/nu[i,:]  *( scy.ive(m-1,y)/scy.ive(m,y) - m/y))
            Mt[2,3,:] = -nu[i+1,:]**2*b[i]/mur[i+1,:]*(
                    mur[i+1,:]/nu[i+1,:]*(-scy.kve(m-1,x)/scy.kve(m,x) - m/x)
                    - mur[i,:]/nu[i,:]  *(-scy.kve(m-1,y)/scy.kve(m,y) - m/y))
            Mt[3,2,:] = -nu[i+1,:]**2*b[i]/mur[i+1,:]*(
                    mur[i+1,:]/nu[i+1,:]*( scy.ive(m-1,x)/scy.ive(m,x) - m/x)
                    - mur[i,:]/nu[i,:]  *( scy.ive(m-1,y)/scy.ive(m,y) - m/y))
            Mt[3,3,:] = -nu[i+1,:]**2*b[i]/mur[i+1,:]*(
                    mur[i+1,:]/nu[i+1,:]*( scy.ive(m-1,x)/scy.ive(m,x) - m/x)
                    - mur[i,:]/nu[i,:]  *(-scy.kve(m-1,y)/scy.kve(m,y) - m/y))
            Mt[2,0,:] = (nu[i+1,:]**2/nu[i,:]**2 - 1)*m/(bet*mur[i+1,:])
            Mt[2,1,:] = Mt[2,0,:]
            Mt[3,0,:] = Mt[2,0,:]
            Mt[3,1,:] = Mt[2,0,:]

            if len(b) == 1:
                M = Mt
            else:
                if (i ==0):
                    M = produto(D,Mt)
                elif i < len(b)-1:
                    M = produto(D,produto(Mt,M))
                else:
                    M = produto(Mt,M)
        return M

    def alphaTM(m, epr, mur, bet, nu, b):
        M = Mtil(m, epr, mur, bet, nu, b)

        B = (M[0,1,:]*M[2,2,:] - M[2,1,:]*M[0,2,:]) / (M[0,0,:]*M[2,2,:] - M[0,2,:]*M[2,0,:])
        alphaTM = scy.kv(m,nu[0,:]*b[0])/scy.iv(m,nu[0,:]*b[0]) * B
        return alphaTM

    ####################
    gam = E/511e3
    bet = np.sqrt(1-1/gam**2)
    nu  = np.ones((epr.shape[0],1))*abs(w/c)*np.sqrt(1 - bet**2*epr*mur)

    Zl = 1j*L*w   /(2*np.pi*ep0 * (bet*c)**2*gam**2)*alphaTM(0, epr, mur, bet, nu, b)
    Zv = 1j*L*w**2/(4*np.pi*ep0*c**2*(bet*c)*gam**4)*alphaTM(1, epr, mur, bet, nu, b)

    # The code cant handle w = 0;
    ind0, = np.where(w == 0)
    if ind0:
        Zl[ind0] = (Zl[ind0+1] + Zl[ind0-1]).real/2
        Zv[ind0] = 0 + 1j*Zv[ind0+1].imag


    # n=10;
    # fil   = exp(-((-n:n)/(n/5)).^2)/sqrt(pi)/n*5;
    # Zv = conv(Zv,fil,'same');
    Zh = Zv

    return Zl.conj(), Zv.conj(), Zh.conj()

def ferrite_kicker_impedance(w,a,b,d,L,epr,mur,Zg, model, coupled):
    # Calculates Impedances for a ferrite kicker:
    #       PIOR CASO:
    #       multi-layer cilindrica com vacuo-coating(2microns)-ceramica-ferrite-PEC
    #
    #       MELHOR CASO:
    #       multi-layer cilindrica com vacuo-coating(2microns)-ceramica-PEC
    #
    # Inputs:
    #
    # w   = vector of angular frequencies to evaluate impedances [rad/s]
    # epr = vector with real and imaginary electric permeability of ferrite for
    #       the frequency range of interest
    # mur = vector with real and imaginary magnetic permissivity of ferrite for
    #       the frequency range of interest
    # L   = length of the structure [m]
    # Zg  = Vector with generator impedance in function of w [Ohm]
    #
    # Outputs:
    #
    # Zl = Impedancia Longitudinal [Ohm]
    # Zh = Impedancia Horizontal [Ohm/m]
    # Zv = Impedancia Vertical   [Ohm/m]
    #
    #

    epb     = np.array([1, 1, 9.3, 1, 12, 1],dtype=float)
    mub     = np.array([1, 1, 1, 1, 1, 1],dtype=float)
    ange    = np.array([0, 0, 0, 0, 0, 0],dtype=float)
    angm    = np.array([0, 0, 0, 0, 0, 0],dtype=float)
    sigmadc = np.array([0, 2.4e6,1,1,1, 5.9e7],dtype=float)
    tau     = np.array([0, 0, 0, 0, 0, 0],dtype=float)*27e-15
    b1      = np.array([(b - 2.0e-3 - 10e-6), (b - 2.0e-3), (b-1.0e-3), b , d],dtype=float)

    epr1 = np.zeros((epb.shape[0],w.shape[0]),dtype=complex)
    mur1 = np.zeros((epb.shape[0],w.shape[0]),dtype=complex)
    for j in range(len(epb)):
        epr1[j,:] = epb[j]*(1-1j*np.sign(w)*np.tan(ange[j])) + sigmadc[j]/(1+1j*w*tau[j])/(1j*w*ep0)
        mur1[j,:] = mub[j]*(1-1j*np.sign(w)*np.tan(angm[j]))
    epr1[5,:] = epr
    mur1[5,:] = mur

    if model.startswith('tsutsui'):
        Zl, Zh, Zv = kicker_tsutsui_model(w, epr, mur, a, b, d, L, 10)
    elif model.startswith('pior'):
        Zl, Zv, Zh = resistive_multilayer_round_pipe(w, epr1, mur1, b1, L, 3)
        Zv = np.pi**2/12*Zv
        Zh = np.pi**2/24*Zh
        Zqh = -Zh
        Zqv = Zh
    else:
        indx = [0, 1, 2, 3, 5]
        mur1 = mur1[indx,:]
        epr1 = epr1[indx,:]
        b1    = b1([0, 1, 2, 3])
        Zl, Zv, Zh = resistive_multilayer_round_pipe(w, epr1, mur1, b1, L, 3)
        Zv = np.pi**2/12*Zv
        Zh = np.pi**2/24*Zh
        Zqh = -Zh
        Zqv = Zh

def kicker_coupled_flux(w,h,W,t,L,mur,Zg):
    # Calculates Impedances for a ferrite kicker:
    #   - For the Coupled Flux, it uses Davino-Hahn model.
    #
    #   DAVINO-HAHN MODEL:
    #
    #  #######################################    |
    #  ###############FERRITE#################    t
    #  #######################################    |
    #  ###**                             **###  |
    #  ###**  VACUUM                     **###  |      ______|  |_________|
    #  ###**                             **###  |            |  |         |
    #  ###**             +   .           **###  w            |  |         |
    #  ###**                             **###  |            )||(         \
    #  ###**             |_D_|           **###  |      Zk  L1)||(L2     Zg/
    #  ###**                             **###  |            )||(         \
    #  #######################################               | M|         /
    #  #######################################               |  |         |
    #  #######################################         ______|  |_________|
    #      |______________h______________|
    #
    # Bibliografias:
    #
    # - Davino_D Hahn_H - Improved Analytical Model of the transverse coupling
    #   impedance of ferrite kicker magnets - Phys. Rev. ST-AB v6 012001 2003
    #
    # - Nassibian G Sacherer F - Methods for measuring tranverse coupling
    #   impedances in circular Accelerators - Nucl Inst and Meth. 159 21-27 1979

    # Equivalent Circuit model.
    D = 0.5e-3
    M  = L*D*mu0/W
    #     L2 = L*2*a*mu0/2/b
    L2 = L*h*mu0/W*(mur*t/(mur*t+h*(h/W+1)))

    Zk =      w * (M/L2)**2 * Zg*L2*1j/(1j*w*L2 + Zg)
    Zx = c/D**2 * (M/L2)**2 * Zg*L2*1j/(1j*w*L2 + Zg)

    return Zk.conj(), Zx.conj()  # take the conjugate to adapt impedance convention

def kicker_tsutsui_model(w, epr, mur, a, b, d, L, n):
    #   - For the Uncoupled Flux, we can choose between three models:
    #
    #       TSUTSUI MODEL:
    #
    #  ******************PEC**********************
    #  *******************************************
    #  **#######################################**      |
    #  **################FERRITE################**      |
    #  **#######################################**      |
    #  **                                       **  |   d
    #  **                                       **  b   |
    #  **                                       **  |   |
    #  **     VACUUM        .                   **  |   |
    #  **                                       **
    #  **                                       **
    #  **                                       **
    #  **#######################################**
    #  **#######################################**
    #  **#######################################**
    #  *******************************************
    #  *******************************************
    #                       |__________a________|
    #
    # Inputs:
    #
    # w   = vector of angular frequencies to evaluate impedances [rad/s]
    # epr = vector with real and imaginary electric permeability of ferrite for
    #       the frequency range of interest
    # mur = vector with real and imaginary magnetic permissivity of ferrite for
    #       the frequency range of interest
    # n   = max order of terms to sum
    # L   = length of the structure [m]
    #
    # Outputs:
    #
    # Zl = Impedancia Longitudinal [Ohm]
    # Zh = Impedancia Horizontal [Ohm/m]
    # Zv = Impedancia Vertical   [Ohm/m]
    #
    # Bibliografias:
    #
    # - Tsutsui_H - Some Simplified Models of Ferrite Kicker Magnet for
    #   Calculation of longitudinal Coupling Impedance - CERN-SL-2000-004
    #
    # - Tsutsui_H - Transverse Coupling Impedance of a Simplified Ferrite
    #   Kicker Magnet Model - LHC Project Note 234 - 2000

    # Valores do artigo do Wang et al para testar a implementacao das formulas
    # do modelo do Tsutui.
    # a = 103e-3
    # b = 67e-3
    # d = 122e-3
    # L = 1.0
    # Valores do artigo do Tsutsui para testar a implementacao das formulas
    # a = 20e-3
    # b = 16e-3
    # d = 66e-3
    # L = 0.6

    # Terms for the infinite sum:
    n = np.arange(0,n+1)

    k = np.ones(n.shape)*w/c
    epr = np.ones(n.shape)*epr
    mur = np.ones(n.shape)*mur


    kxn = np.repeats((2*n[:,None]+1)*np.pi/2/a, w.shape[0], axis=1)
    kyn = np.sqrt((epr*mur-1)*k**2 - kxn**2)
    sh  = np.sinh(kxn*b)
    ch  = np.cosh(kxn*b)
    tn  = np.tan(kyn*(b-d))
    ct  = 1/np.tan(kyn*(b-d))

    Zl = 1j*Z0*L/2/a / (
        (kxn/k*sh*ch*(1+epr*mur) + kyn/k*(mur*sh**2*tn - epr*ch**2*ct)
        )/(epr*mur-1) - k/kxn*sh*ch)
    Zl = Zl.sum(0)

    Zv = 1j*Z0*L/2/a * kxn**2/k/(
        (kxn/k*sh*ch*(1+epr*mur) + kyn/k*(mur*ch**2*tn - epr*sh**2*ct)
        )/(epr*mur-1) - k/kxn*sh*ch)
    Zv = Zv.sum(0)

    Zh = 1j*Z0*L/2/a * kxn**2/k/(
        (kxn/k*sh*ch*(1+epr*mur) + kyn/k*(mur*sh**2*tn - epr*ch**2*ct)
        )/(epr*mur-1) - k/kxn*sh*ch)
    Zh = Zh.sum(0)

    return Zl.conj(), Zh.conj(), Zv.conj() # impedance convention
