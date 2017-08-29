#!/usr/bin/env python3
import os as _os
import numpy as _np
import matplotlib.pyplot as _plt
import scipy.special as _scy
import mathphys as _mp

_c   = _mp.constants.light_speed
_mu0 = _mp.constants.vacuum_permeability
_ep0 = _mp.constants.vacuum_permitticity
_Z0  = _mp.constants.vacuum_impedance
E0   = _mp.constants.electron_rest_energy * _mp.units.joule_2_eV

_IMPS  = {'Zll','Zdy','Zdx','Zqy','Zqx'}
_WAKES = {'Wll','Wdy','Wdx','Wqy','Wqx'}
_TITLE = {'Zll':'Longitudinal Impedance',
          'Zdy':'Driving Vertical Impedance',
          'Zdx':'Driving Horizontal Impedance',
          'Zqy':'Detuning Vertical Impedance',
          'Zqx':'Detuning Horizontal Impedance'}
_FACTOR ={'Zll':1e-3, 'Zdy':1e-3, 'Zdx':1e-3, 'Zqy':1e-3, 'Zqx':1e-3,
          'Wll':1e-3, 'Wdy':1e-6, 'Wdx':1e-6, 'Wqy':1e-6, 'Wqx':1e-6}
_BETA   ={'Zll' :lambda x:1,
          'Zdy':lambda x:x.betay,
          'Zdx':lambda x:x.betax,
          'Zqy':lambda x:x.betay,
          'Zqx':lambda x:x.betax}

RW_default_w = _np.logspace(5,12.4,3000)
RW_default_w = _np.array([-_np.flipud(RW_default_w),RW_default_w]).ravel()


def _plotlog(x, y, color=None, label=None, ax=None,linewidth=1.5):
    if ax is None: ax = _plt.gca()

    if any(y > 0):
        ax.loglog(x,y,color=color,label=label,linewidth=linewidth)
        ax.loglog(x,-y,'--',color=color,linewidth=linewidth)
    else:
        ax.loglog(x,-y,'--',color=color,linewidth=linewidth)
        ax.loglog(x,y,color=color,label=label,linewidth=linewidth)

def _prepare_props(props):
    if isinstance(props,str):
        if props.lower() == 'all':
            props = sorted(list(_IMPS))
        elif props in _IMPS:
            props = [props]
        else:
            raise AttributeError(props,' not supported for plot.')
    elif isinstance(props,(list,tuple)):
        wrong = set(props) - _IMPS
        if wrong:
            raise AttributeError(wrong,' not supported for plot.')
    else:
        raise TypeError("Type '"+type(props)+"' not supported for 'props'.")
    return props

class Element:

    _YLABEL ={'Zll' :r'$Z_l [\Omega]$',
              'Zdy':r'$Z_y^D [k\Omega/m]$',
              'Zdx':r'$Z_x^D [k\Omega/m]$',
              'Zqy':r'$Z_y^Q [k\Omega/m]$',
              'Zqx':r'$Z_x^Q [k\Omega/m]$'}

    def __init__(self,name=None, path=None, betax=None,betay=None,quantity=None):
        self.name     = name or 'element'
        if path is not None: path = _os.path.abspath(path)
        self.path     = path or _os.path.abspath('.')
        self.quantity = quantity or 0 # this field shall only be used in Budget
        self.betax    = betax or 0.0  # this field shall only be used in Budget
        self.betay    = betay or 0.0  # this field shall only be used in Budget
        self.w        = _np.array([],dtype=float)
        self.Zll      = _np.array([],dtype=complex)
        self.Zdy      = _np.array([],dtype=complex)
        self.Zdx      = _np.array([],dtype=complex)
        self.Zqy      = _np.array([],dtype=complex)
        self.Zqx      = _np.array([],dtype=complex)
        self.s        = _np.array([],dtype=float)
        self.Wll      = _np.array([],dtype=float)
        self.Wdy      = _np.array([],dtype=float)
        self.Wdx      = _np.array([],dtype=float)
        self.Wqy      = _np.array([],dtype=float)
        self.Wqx      = _np.array([],dtype=float)

    def copy(self):
        other = Element(name =self.name,path=self.path,betax=self.betax,
                        betay=self.betay,quantity=self.quantity)
        other.w   = self.w.copy()
        other.Zll = self.Zll.copy()
        other.Zdy = self.Zdy.copy()
        other.Zdx = self.Zdx.copy()
        other.Zqy = self.Zqy.copy()
        other.Zqx = self.Zqx.copy()
        other.s   = self.s.copy()
        other.Wll = self.Wll.copy()
        other.Wdy = self.Wdy.copy()
        other.Wdx = self.Wdx.copy()
        other.Wqy = self.Wqy.copy()
        other.Wqx = self.Wqx.copy()

        return other

    def save(self):
        name = self.name.replace(' ','_').lower()
        _mp.utils.save_pickle(_os.path.sep.join([self.path, name]),element=self)

    def load(self):
        name = self.name.replace(' ','_').lower()
        data = _mp.utils.load_pickle(_os.path.sep.join([self.path, name]))
        return data['element']

    def plot(self, props='all', logscale=True, show = True, save = False, name='',figsize=(8,4)):

        if name: name = '_'+name
        props = _prepare_props(props)

        for prop in props:
            Imp2 = getattr(self,prop)
            if Imp2 is None or len(Imp2)==0: continue
            _plt.figure(figsize=figsize)
            Imp = Imp2*_FACTOR[prop]
            w = self.w
            if logscale:
                _plotlog(w, Imp.real,color='b',label='Real')
                _plotlog(w, Imp.imag,color='r',label='Imag')
            else:
                _plt.plot(w,Imp.real,'b',label='Real')
                _plt.plot(w,Imp.imag,'r',label='Imag')
            _plt.legend(loc='best')
            _plt.grid(True)
            _plt.xlabel(r'$\omega [rad/s]$')
            _plt.ylabel(Element._YLABEL[prop])
            _plt.title(self.name+': '+_TITLE[prop])
            if save: _plt.savefig(_os.path.sep.join((self.path, prop + name + '.svg')))
        if show: _plt.show()

class Budget(list):

    _YLABEL ={'Zll':r'$Z_l [k\Omega]$',
              'Zdy':r'$\beta_yZ_y^D [k\Omega]$',
              'Zdx':r'$\beta_xZ_x^D [k\Omega]$',
              'Zqy':r'$\beta_yZ_y^Q [k\Omega]$',
              'Zqx':r'$\beta_xZ_x^Q [k\Omega]$'}

    def __init__(self, lista=None, name=None, path = None):
        lista = lista or []
        if lista and not isinstance(lista[0],Element):
            assert 'Input must be a sequence of Element objects.'
        super().__init__(lista)
        self.name = name or 'Budget'
        if path is not None: path = _os.path.abspath(path)
        self.path = path or _os.path.abspath('.')

    def __setitem__(self,k,v):
        assert isinstance(v,Element)
        super().__setitem__(k,v)

    def __setattr__(self,name,value):
        if name in {'name','path'}:
            self.__dict__[name] = str(value)
        else:
            raise AttributeError('Attribute '+name+' is read only.')

    def __getattr__(self,name):
        if name not in _IMPS | _WAKES | {'w','s'}:
            return [getattr(x,name) for x in self]

        w = _np.unique(_np.concatenate([getattr(x,'w') for x in self]))
        if name == 'w': return w
        if name in _IMPS:
            temp = _np.zeros(w.shape,dtype=complex)
            for el in self:
                attr = getattr(el,name)
                if attr is None or len(attr) == 0: continue
                temp += 1j*_np.interp(w,el.w,attr.imag,left=0.0,right=0.0)*el.quantity*_BETA[name](el)
                temp +=    _np.interp(w,el.w,attr.real,left=0.0,right=0.0)*el.quantity*_BETA[name](el)
            return temp

        s = _np.unique(_np.concatenate([getattr(x,'s') for x in self]))
        if name == 's': return s
        if name in _WAKES:
            temp = _np.zeros(s.shape,dtype=float)
            for el in self:
                attr = getattr(el,name)
                if attr is None or len(attr) == 0: continue
                temp += _np.interp(s,el.s,attr,left=0.0,right=0.0)*el.quantity*_BETA[name](el)
            return temp
        raise AttributeError("'"+self.__class__.__name__+ "' object has no attribute '"+name+"'" )

    def __str__(self):
        string  = '{0:^48s}\n'.format(self.name)
        string += '{0:^15s}: {1:^10s} {2:^10s} {3:^10s}\n'.format('Element','Quantity','Betax','Betay')
        for el in self:
            string += '{0:<15s}: {1:^10d} {2:^10.1f} {3:^10.1f}\n'.format(el.name,el.quantity,el.betax,el.betay)
        string +='\n'
        return string

    def copy(self):
        other = Budget(name=self.name,path=self.path)
        for el in self:
            other.append(el.copy())
        return other

    def budget2element(self,name=None,path=None):
        ele = Element(name=name,path=path)
        for prop in _IMPS | _WAKES | {'w','s'}:
            Imp2 = getattr(self,prop)
            if not _np.isclose(Imp2,0).all():
                setattr(ele,prop,Imp2.copy())
        ele.betax = 1.0
        ele.betay = 1.0
        ele.quantity = 1
        return ele

    def save(self):
        name = self.name.replace(' ','_').lower()
        _mp.utils.save_pickle(_os.path.sep.join([self.path,name]),budget=self)

    def load(self):
        name = self.name.replace(' ','_').lower()
        data = _mp.utils.load_pickle(_os.path.sep.join([self.path,name]))
        return data['budget']

    def plot(self, props='all', logscale=True, show = True, save = False,
            name='',figsize=(8,6),fontsize=14,linewidth=1.5):

        color_map = _plt.get_cmap('nipy_spectral')
        if name: name = '_'+name
        props = _prepare_props(props)

        for prop in props:
            a = True
            for el in self:
                Imp3 = getattr(el,prop)
                a &= Imp3 is None or len(Imp3)==0
                if not a: break
            if a: continue
            f,ax = _plt.subplots(2,1,sharex=True,figsize=figsize)
            N = len(self)
            for i,el in enumerate(self):
                Imp2 = getattr(el,prop)
                if Imp2 is None or len(Imp2)==0: continue
                Imp = Imp2*_FACTOR[prop] * el.quantity * _BETA[prop](el)
                w = el.w
                cor = color_map(i/N)
                if logscale:
                    _plotlog(w, Imp.real, color=cor, ax=ax[0],linewidth=linewidth)
                    _plotlog(w, Imp.imag, color=cor, label=el.name, ax=ax[1],linewidth=linewidth)
                else:
                    ax[0].plot(w,Imp.real,color=cor,linewidth=linewidth)
                    ax[1].plot(w,Imp.imag,color=cor,label=el.name,linewidth=linewidth)
            ax[1].legend(loc='best',fontsize=10)
            ax[0].grid(True)
            ax[0].tick_params(labelsize=fontsize)
            ax[1].grid(True)
            ax[1].tick_params(labelsize=fontsize)
            ax[1].set_xlabel(r'$\omega [rad/s]$',fontsize=fontsize)
            ax[0].set_ylabel(r'Re'+Budget._YLABEL[prop],fontsize=fontsize)
            ax[1].set_ylabel(r'Im'+Budget._YLABEL[prop],fontsize=fontsize)
            ax[0].set_title(self.name+': '+_TITLE[prop],fontsize=fontsize)
            if save: f.savefig(_os.path.sep.join((self.path, prop + name + '.svg')))
        if show: _plt.show()

def load_budget(fname):
        data = _mp.utils.load_pickle(fname)
        return data['budget']

def load_element(fname):
        data = _mp.utils.load_pickle(fname)
        return data['element']

def longitudinal_resonator(Rs, Q, wr, w):
    """Returns the longitudinal resonator impedance for w.

    Inputs:
    Rs    = Shunt impedance [Ohm]
    Q     = Quality Factor
    wr    = Angular resonance frequency [rad/s]
    w     = numpy array of frequencies to calculate the impedance [rad/s]

    Rs, Q and wr may be a float, list, tuple or a numpy array.

    Outputs:
    Zl    = Longitudinal Impedance (real and imaginary) [Ohm]

    Sign convention followed:
        Chao, A., Physics of Collective Beam Instabilities in High Energy
        Accelerators, Wiley 1993.

    Examples:
    >>> w = _np.linspace(-5e9,5e9,1000)
    >>> Rs, Q, wr = 1000, 1, 2*_np.pi*1e9
    >>> Zl = longitudinal_resonator(Rs,Q,wr,w)
    >>> Zl.shape
    (1000,)
    >>> Rs, Q, wr = [1000,2000], [1,10], [2*_np.pi*1e9,2*_np.pi*5e8]
    >>> Zl = longitudinal_resonator(Rs,Q,wr,w)
    >>> Zl.shape
    (1000,)
    >>> Rs, Q, wr = _np.array(Rs), _np.array(Q), _np.array(wr)
    >>> Zl = longitudinal_resonator(Rs,Q,wr,w)
    >>> Zl.shape
    (1000,)
    """
    Rs = _np.array(Rs,ndmin=1,dtype=float)[:,None] # I am using broadcasting
    Q  = _np.array(Q, ndmin=1,dtype=float)[:,None]
    wr = _np.array(wr,ndmin=1,dtype=float)[:,None]
    Zl = w*Rs / (w+1j*Q*(wr - w**2/wr))
    return Zl.sum(0).flatten()

def transverse_resonator(Rs, Q, wr, w):
    """Returns the transverse resonator impedance for w.

    Inputs:
    Rs    = Shunt impedance [Ohm]
    Q     = Quality Factor
    wr    = Angular resonance frequency [rad/s]
    w     = numpy array of frequencies to calculate the impedance [rad/s]

    Rs, Q and wr may be a float, list, tuple or a numpy array.

    Outputs:
    Zl    = Longitudinal Impedance (real and imaginary) [Ohm]

    Sign convention followed:
        Chao, A., Physics of Collective Beam Instabilities in High Energy
        Accelerators, Wiley 1993.

    Examples:
    >>> w = _np.linspace(-5e9,5e9,1000)
    >>> Rs, Q, wr = 1000, 1, 2*_np.pi*1e9
    >>> Zt = longitudinal_resonator(Rs,Q,wr,w)
    >>> Zt.shape
    (1000,)
    >>> Rs, Q, wr = [1000,2000], [1,10], [2*_np.pi*1e9,2*_np.pi*5e8]
    >>> Zt = longitudinal_resonator(Rs,Q,wr,w)
    >>> Zt.shape
    (1000,)
    >>> Rs, Q, wr = _np.array(Rs), _np.array(Q), _np.array(wr)
    >>> Zt = longitudinal_resonator(Rs,Q,wr,w)
    >>> Zt.shape
    (1000,)
    """
    Rs = _np.array(Rs,ndmin=1,dtype=float)[:,None] # I am using broadcasting
    Q  = _np.array(Q, ndmin=1,dtype=float)[:,None]
    wr = _np.array(wr,ndmin=1,dtype=float)[:,None]
    Zt = wr*Rs/(w + 1j*Q*(wr - w**2/wr))
    return Zt.sum(0).flatten()

def wake_longitudinal_resonator(Rs, Q, wr, s):
    """Returns the longitudinal resonator wake-function for s.

    Inputs:
    Rs    = Shunt impedance [Ohm]
    Q     = Quality Factor
    wr    = Angular resonance frequency [rad/s]
    w     = numpy array of frequencies to calculate the impedance [rad/s]

    Rs, Q and wr may be a float, list, tuple or a numpy array.

    Outputs:
    Zl    = Longitudinal Impedance (real and imaginary) [Ohm]

    Sign convention followed:
        Chao, A., Physics of Collective Beam Instabilities in High Energy
        Accelerators, Wiley 1993.

    Examples:
    >>> w = _np.linspace(-5e9,5e9,1000)
    >>> Rs, Q, wr = 1000, 1, 2*_np.pi*1e9
    >>> Zl = longitudinal_resonator(Rs,Q,wr,w)
    >>> Zl.shape
    (1000,)
    >>> Rs, Q, wr = [1000,2000], [1,10], [2*_np.pi*1e9,2*_np.pi*5e8]
    >>> Zl = longitudinal_resonator(Rs,Q,wr,w)
    >>> Zl.shape
    (1000,)
    >>> Rs, Q, wr = _np.array(Rs), _np.array(Q), _np.array(wr)
    >>> Zl = longitudinal_resonator(Rs,Q,wr,w)
    >>> Zl.shape
    (1000,)
    """
    Rs    = _np.array(Rs,ndmin=1,dtype=float)[:,None] # I am using broadcasting
    Q     = _np.array(Q, ndmin=1,dtype=float)[:,None]
    wr    = _np.array(wr,ndmin=1,dtype=float)[:,None]
    alpha   = wr / (2*Q)
    wrl     = _np.sqrt(wr*wr - alpha*alpha)
    Wl      = _np.zeros(len(s))
    sel     = s > 0.0
    Wl[sel] = (2*alpha * Rs * _np.exp(-alpha*s[sel]/_c)*(_np.cos(wrl*s[sel]/_c) -
                                               alpha/wrl*_np.sin(wrl*s[sel]/_c))
            ).sum(0).flatten()
    Wl[s==0] = (alpha * Rs) .sum(0).flatten()
    return Wl

def wake_transverse_resonator(Rs, Q, wr, s):
    """Returns the Transverse resonator wake-function for w.

    Inputs:
    Rs    = Shunt impedance [Ohm]
    Q     = Quality Factor
    wr    = Angular resonance frequency [rad/s]
    w     = numpy array of frequencies to calculate the impedance [rad/s]

    Rs, Q and wr may be a float, list, tuple or a numpy array.

    Outputs:
    Zl    = Longitudinal Impedance (real and imaginary) [Ohm]

    Sign convention followed:
        Chao, A., Physics of Collective Beam Instabilities in High Energy
        Accelerators, Wiley 1993.

    Examples:
    >>> w = _np.linspace(-5e9,5e9,1000)
    >>> Rs, Q, wr = 1000, 1, 2*_np.pi*1e9
    >>> Zt = longitudinal_resonator(Rs,Q,wr,w)
    >>> Zt.shape
    (1000,)
    >>> Rs, Q, wr = [1000,2000], [1,10], [2*_np.pi*1e9,2*_np.pi*5e8]
    >>> Zt = longitudinal_resonator(Rs,Q,wr,w)
    >>> Zt.shape
    (1000,)
    >>> Rs, Q, wr = _np.array(Rs), _np.array(Q), _np.array(wr)
    >>> Zt = longitudinal_resonator(Rs,Q,wr,w)
    >>> Zt.shape
    (1000,)
    """
    Rs = _np.array(Rs,ndmin=1,dtype=float)[:,None] # I am using broadcasting
    Q  = _np.array(Q, ndmin=1,dtype=float)[:,None]
    wr = _np.array(wr,ndmin=1,dtype=float)[:,None]
    alpha   = wr / (2*Q)
    wrl     = _np.sqrt(wr*wr - alpha*alpha)
    Wt      = _np.zeros(len(s))
    sel     = s > 0.0
    Wt[sel] = (Rs * wr / (Q * wrl) * _np.exp(-alpha*s[sel]/_c)*_np.sin(wrl*s[sel]/_c)
              ).sum(0).flatten()
    return Wt


def _prepare_input_epr_mur(w,epb,mub,ange,angm,sigmadc,tau):
    epr = _np.zeros((len(epb),len(w)),dtype=complex)
    mur = _np.zeros((len(epb),len(w)),dtype=complex)
    for j in range(len(epb)):
        epr[j,:] = epb[j]*(1-1j*_np.sign(w)*_np.tan(ange[j])) + sigmadc[j]/(1+1j*w*tau[j])/(1j*w*_ep0)
        mur[j,:] = mub[j]*(1-1j*_np.sign(w)*_np.tan(angm[j]))
    return epr, mur

def resistive_multilayer_round_pipe(w,epr,mur,b,L,E):

    def Mtil(m, epr, mur, bet, nu, b):
        def produto(A,B):
            C = _np.zeros((A.shape[0],B.shape[1],A.shape[2]),dtype=complex)
            for i in range(A.shape[0]):
                for j in range(B.shape[1]):
                    for k in range(A.shape[1]):
                        C[i,j,:] = C[i,j,:] + A[i,k,:]*B[k,j,:]
            return C

        for i in range(len(b)): # lembrando que length(b) = número de camadas - 1
            x = nu[i+1,:] * b[i]
            y = nu[i  ,:] * b[i]
            Mt = _np.zeros((4,4,w.shape[0]),dtype=complex)

            if i<len(b)-1:
                D = _np.zeros((4,4,nu.shape[1]),dtype=complex)
                z = nu[i+1,:]*b[i+1]
                if not any(z.real<0):
                    ind = (z.real<60)

                    A = _scy.iv(m,z[ind])
                    B = _scy.kv(m,z[ind])
                    C = _scy.iv(m,x[ind])
                    E = _scy.kv(m,x[ind])

                    D[0,0,:]    =  1
                    D[1,1,ind]  = - B*C/(A*E)
                    D[1,1,~ind] = - _np.exp(-2*(z[~ind]-x[~ind]))
                    D[2,2,:]    =  1
                    D[3,3,ind]  = - B*C/(A*E)
                    D[3,3,~ind] = - _np.exp(-2*(z[~ind]-x[~ind]))

            Mt[0,0,:] = -nu[i+1,:]**2*b[i]/epr[i+1,:]*(
                    epr[i+1,:]/nu[i+1,:]*(-_scy.kve(m-1,x)/_scy.kve(m,x) - m/x)
                    - epr[i,:]/nu[i,:]  *( _scy.ive(m-1,y)/_scy.ive(m,y) - m/y))
            Mt[0,1,:] = -nu[i+1,:]**2*b[i]/epr[i+1,:]*(
                    epr[i+1,:]/nu[i+1,:]*(-_scy.kve(m-1,x)/_scy.kve(m,x) - m/x)
                    - epr[i,:]/nu[i,:]  *(-_scy.kve(m-1,y)/_scy.kve(m,y) - m/y))
            Mt[1,0,:] = -nu[i+1,:]**2*b[i]/epr[i+1,:]*(
                    epr[i+1,:]/nu[i+1,:]*( _scy.ive(m-1,x)/_scy.ive(m,x) - m/x)
                    - epr[i,:]/nu[i,:]  *( _scy.ive(m-1,y)/_scy.ive(m,y) - m/y))
            Mt[1,1,:] = -nu[i+1,:]**2*b[i]/epr[i+1,:]*(
                    epr[i+1,:]/nu[i+1,:]*( _scy.ive(m-1,x)/_scy.ive(m,x) - m/x)
                    - epr[i,:]/nu[i,:]  *(-_scy.kve(m-1,y)/_scy.kve(m,y) - m/y))

            Mt[0,2,:] = (nu[i+1,:]**2/nu[i,:]**2 - 1)*m/(bet*epr[i+1,:])
            Mt[0,3,:] = Mt[0,2,:]
            Mt[1,2,:] = Mt[0,2,:]
            Mt[1,3,:] = Mt[0,2,:]
            Mt[2,2,:] = -nu[i+1,:]**2*b[i]/mur[i+1,:]*(
                    mur[i+1,:]/nu[i+1,:]*(-_scy.kve(m-1,x)/_scy.kve(m,x) - m/x)
                    - mur[i,:]/nu[i,:]  *( _scy.ive(m-1,y)/_scy.ive(m,y) - m/y))
            Mt[2,3,:] = -nu[i+1,:]**2*b[i]/mur[i+1,:]*(
                    mur[i+1,:]/nu[i+1,:]*(-_scy.kve(m-1,x)/_scy.kve(m,x) - m/x)
                    - mur[i,:]/nu[i,:]  *(-_scy.kve(m-1,y)/_scy.kve(m,y) - m/y))
            Mt[3,2,:] = -nu[i+1,:]**2*b[i]/mur[i+1,:]*(
                    mur[i+1,:]/nu[i+1,:]*( _scy.ive(m-1,x)/_scy.ive(m,x) - m/x)
                    - mur[i,:]/nu[i,:]  *( _scy.ive(m-1,y)/_scy.ive(m,y) - m/y))
            Mt[3,3,:] = -nu[i+1,:]**2*b[i]/mur[i+1,:]*(
                    mur[i+1,:]/nu[i+1,:]*( _scy.ive(m-1,x)/_scy.ive(m,x) - m/x)
                    - mur[i,:]/nu[i,:]  *(-_scy.kve(m-1,y)/_scy.kve(m,y) - m/y))
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
        alphaTM = _scy.kv(m,nu[0,:]*b[0])/_scy.iv(m,nu[0,:]*b[0]) * B
        return alphaTM

    ####################
    gam = E/E0
    bet = _np.sqrt(1-1/gam**2)
    nu  = _np.ones((epr.shape[0],1))*abs(w/_c)*_np.sqrt(1 - bet**2*epr*mur)

    Zl = 1j*L*w   /(2*_np.pi*_ep0 * (bet*_c)**2*gam**2)*alphaTM(0, epr, mur, bet, nu, b)
    Zv = 1j*L*w**2/(4*_np.pi*_ep0*_c**2*(bet*_c)*gam**4)*alphaTM(1, epr, mur, bet, nu, b)

    # The code cant handle w = 0;
    ind0, = _np.where(w == 0)
    if ind0:
        Zl[ind0] = (Zl[ind0+1] + Zl[ind0-1]).real/2
        Zv[ind0] = 0 + 1j*Zv[ind0+1].imag


    # n=10;
    # fil   = exp(-((-n:n)/(n/5)).^2)/sqrt(pi)/n*5;
    # Zv = conv(Zv,fil,'same');
    Zh = Zv.copy()

    return Zl.conj(), Zv.conj(), Zh.conj()

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
    M  = L*D*_mu0/W
    #     L2 = L*2*a*_mu0/2/b
    L2 = L*h*_mu0/W*(mur*t/(mur*t+h*(h/W+1)))

    Zk =      w * (M/L2)**2 * Zg*L2*1j/(1j*w*L2 + Zg)
    Zx = _c/D**2 * (M/L2)**2 * Zg*L2*1j/(1j*w*L2 + Zg)

    return Zk.conj(), Zx.conj()  # take the conjugate to adapt impedance convention

def kicker_tsutsui_model(w, epr, mur, a, b, d, L, n):
    '''
       - For the Uncoupled Flux, we can choose between three models:

           TSUTSUI MODEL:

      ******************PEC**********************
      *******************************************
      **#######################################**      |
      **################FERRITE################**      |
      **#######################################**      |
      **                                       **  |   d
      **                                       **  b   |
      **                                       **  |   |
      **     VACUUM        .                   **  |   |
      **                                       **
      **                                       **
      **                                       **
      **#######################################**
      **#######################################**
      **#######################################**
      *******************************************
      *******************************************
                           |__________a________|

     Inputs:

     w   = vector of angular frequencies to evaluate impedances [rad/s]
     epr = vector with real and imaginary electric permeability of ferrite for
           the frequency range of interest
     mur = vector with real and imaginary magnetic permissivity of ferrite for
           the frequency range of interest
     n   = max order of terms to sum
     L   = length of the structure [m]

     Outputs:

     Zl = Impedancia Longitudinal [Ohm]
     Zh = Impedancia Horizontal [Ohm/m]
     Zv = Impedancia Vertical   [Ohm/m]

     Bibliografias:

     - Tsutsui_H - Some Simplified Models of Ferrite Kicker Magnet for
       Calculation of longitudinal Coupling Impedance - CERN-SL-2000-004

     - Tsutsui_H - Transverse Coupling Impedance of a Simplified Ferrite
       Kicker Magnet Model - LHC Project Note 234 - 2000
'''

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
    n = _np.arange(0,n+1)[:,None]

    k = _np.ones(n.shape)*w/_c
    epr = _np.ones(n.shape)*epr
    mur = _np.ones(n.shape)*mur


    kxn = _np.repeat((2*n + 1)*_np.pi/2/a, w.shape[0], axis=1)
    kyn = _np.sqrt((epr*mur-1)*k**2 - kxn**2)
    sh  = _np.sinh(kxn*b)
    ch  = _np.cosh(kxn*b)
    tn  = _np.tan(kyn*(b-d))
    ct  = 1/_np.tan(kyn*(b-d))

    Zl = 1j*_Z0*L/2/a / (
        (kxn/k*sh*ch*(1+epr*mur) + kyn/k*(mur*sh**2*tn - epr*ch**2*ct)
        )/(epr*mur-1) - k/kxn*sh*ch)
    Zl = Zl.sum(0)

    Zv = 1j*_Z0*L/2/a * kxn**2/k/(
        (kxn/k*sh*ch*(1+epr*mur) + kyn/k*(mur*ch**2*tn - epr*sh**2*ct)
        )/(epr*mur-1) - k/kxn*sh*ch)
    Zv = Zv.sum(0)

    Zh = 1j*_Z0*L/2/a * kxn**2/k/(
        (kxn/k*sh*ch*(1+epr*mur) + kyn/k*(mur*sh**2*tn - epr*ch**2*ct)
        )/(epr*mur-1) - k/kxn*sh*ch)
    Zh = Zh.sum(0)

    return Zl.conj(), Zh.conj(), Zv.conj() # impedance convention

def yokoya_factors(plane='ll'):
    if plane == 'll':
        return 1
    elif plane == 'dy':
        return _np.pi**2/12
    elif plane in {'qy','dx'}:
        return _np.pi**2/24
    elif plane == 'qx':
        return -_np.pi**2/24
    else:
        raise Exception("Plane not identified. Possible options: 'll','dy','dx','qy','qx' ")


def taper(w,R1,R2,L1,L2):
    '''
                   L2
    _____|- - - - - - - - - - - - |_____
          \                      /    :
           \                    /     :
            \                  /      :
             \_______L1_______/       : R2
                :                     :
                : R1                  :
                :                     :
    - - - - - - - - - - - - - - - - - - -
    '''
    theta = (L2-L1)/2 / (R2-R1)