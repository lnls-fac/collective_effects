# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.2
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_cppcolleff', [dirname(__file__)])
        except ImportError:
            import _cppcolleff
            return _cppcolleff
        if fp is not None:
            try:
                _mod = imp.load_module('_cppcolleff', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _cppcolleff = swig_import_helper()
    del swig_import_helper
else:
    import _cppcolleff
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


TWOPI = _cppcolleff.TWOPI
light_speed = _cppcolleff.light_speed

def interpola(*args):
  return _cppcolleff.interpola(*args)
interpola = _cppcolleff.interpola
class Particle_t(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Particle_t, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Particle_t, name)
    __repr__ = _swig_repr
    __swig_setmethods__["xx"] = _cppcolleff.Particle_t_xx_set
    __swig_getmethods__["xx"] = _cppcolleff.Particle_t_xx_get
    if _newclass:xx = _swig_property(_cppcolleff.Particle_t_xx_get, _cppcolleff.Particle_t_xx_set)
    __swig_setmethods__["xl"] = _cppcolleff.Particle_t_xl_set
    __swig_getmethods__["xl"] = _cppcolleff.Particle_t_xl_get
    if _newclass:xl = _swig_property(_cppcolleff.Particle_t_xl_get, _cppcolleff.Particle_t_xl_set)
    __swig_setmethods__["de"] = _cppcolleff.Particle_t_de_set
    __swig_getmethods__["de"] = _cppcolleff.Particle_t_de_get
    if _newclass:de = _swig_property(_cppcolleff.Particle_t_de_get, _cppcolleff.Particle_t_de_set)
    __swig_setmethods__["ss"] = _cppcolleff.Particle_t_ss_set
    __swig_getmethods__["ss"] = _cppcolleff.Particle_t_ss_get
    if _newclass:ss = _swig_property(_cppcolleff.Particle_t_ss_get, _cppcolleff.Particle_t_ss_set)
    def __init__(self, *args): 
        this = _cppcolleff.new_Particle_t(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _cppcolleff.delete_Particle_t
    __del__ = lambda self : None;
Particle_t_swigregister = _cppcolleff.Particle_t_swigregister
Particle_t_swigregister(Particle_t)

class Bunch_t(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Bunch_t, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Bunch_t, name)
    __repr__ = _swig_repr
    __swig_setmethods__["num_part"] = _cppcolleff.Bunch_t_num_part_set
    __swig_getmethods__["num_part"] = _cppcolleff.Bunch_t_num_part_get
    if _newclass:num_part = _swig_property(_cppcolleff.Bunch_t_num_part_get, _cppcolleff.Bunch_t_num_part_set)
    __swig_setmethods__["Ib"] = _cppcolleff.Bunch_t_Ib_set
    __swig_getmethods__["Ib"] = _cppcolleff.Bunch_t_Ib_get
    if _newclass:Ib = _swig_property(_cppcolleff.Bunch_t_Ib_get, _cppcolleff.Bunch_t_Ib_set)
    __swig_setmethods__["particles"] = _cppcolleff.Bunch_t_particles_set
    __swig_getmethods__["particles"] = _cppcolleff.Bunch_t_particles_get
    if _newclass:particles = _swig_property(_cppcolleff.Bunch_t_particles_get, _cppcolleff.Bunch_t_particles_set)
    def __init__(self, *args): 
        this = _cppcolleff.new_Bunch_t(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _cppcolleff.delete_Bunch_t
    __del__ = lambda self : None;
    def generate_bunch(self): return _cppcolleff.Bunch_t_generate_bunch(self)
    def sort(self): return _cppcolleff.Bunch_t_sort(self)
Bunch_t_swigregister = _cppcolleff.Bunch_t_swigregister
Bunch_t_swigregister(Bunch_t)

class Ring_t(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Ring_t, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Ring_t, name)
    __repr__ = _swig_repr
    __swig_setmethods__["betax"] = _cppcolleff.Ring_t_betax_set
    __swig_getmethods__["betax"] = _cppcolleff.Ring_t_betax_get
    if _newclass:betax = _swig_property(_cppcolleff.Ring_t_betax_get, _cppcolleff.Ring_t_betax_set)
    __swig_setmethods__["alphax"] = _cppcolleff.Ring_t_alphax_set
    __swig_getmethods__["alphax"] = _cppcolleff.Ring_t_alphax_get
    if _newclass:alphax = _swig_property(_cppcolleff.Ring_t_alphax_get, _cppcolleff.Ring_t_alphax_set)
    __swig_setmethods__["etax"] = _cppcolleff.Ring_t_etax_set
    __swig_getmethods__["etax"] = _cppcolleff.Ring_t_etax_get
    if _newclass:etax = _swig_property(_cppcolleff.Ring_t_etax_get, _cppcolleff.Ring_t_etax_set)
    __swig_setmethods__["etaxl"] = _cppcolleff.Ring_t_etaxl_set
    __swig_getmethods__["etaxl"] = _cppcolleff.Ring_t_etaxl_get
    if _newclass:etaxl = _swig_property(_cppcolleff.Ring_t_etaxl_get, _cppcolleff.Ring_t_etaxl_set)
    __swig_setmethods__["gammax"] = _cppcolleff.Ring_t_gammax_set
    __swig_getmethods__["gammax"] = _cppcolleff.Ring_t_gammax_get
    if _newclass:gammax = _swig_property(_cppcolleff.Ring_t_gammax_get, _cppcolleff.Ring_t_gammax_set)
    __swig_setmethods__["tunex"] = _cppcolleff.Ring_t_tunex_set
    __swig_getmethods__["tunex"] = _cppcolleff.Ring_t_tunex_get
    if _newclass:tunex = _swig_property(_cppcolleff.Ring_t_tunex_get, _cppcolleff.Ring_t_tunex_set)
    __swig_setmethods__["chromx"] = _cppcolleff.Ring_t_chromx_set
    __swig_getmethods__["chromx"] = _cppcolleff.Ring_t_chromx_get
    if _newclass:chromx = _swig_property(_cppcolleff.Ring_t_chromx_get, _cppcolleff.Ring_t_chromx_set)
    __swig_setmethods__["tunex_shift"] = _cppcolleff.Ring_t_tunex_shift_set
    __swig_getmethods__["tunex_shift"] = _cppcolleff.Ring_t_tunex_shift_get
    if _newclass:tunex_shift = _swig_property(_cppcolleff.Ring_t_tunex_shift_get, _cppcolleff.Ring_t_tunex_shift_set)
    __swig_setmethods__["circum"] = _cppcolleff.Ring_t_circum_set
    __swig_getmethods__["circum"] = _cppcolleff.Ring_t_circum_get
    if _newclass:circum = _swig_property(_cppcolleff.Ring_t_circum_get, _cppcolleff.Ring_t_circum_set)
    __swig_setmethods__["mom_comp"] = _cppcolleff.Ring_t_mom_comp_set
    __swig_getmethods__["mom_comp"] = _cppcolleff.Ring_t_mom_comp_get
    if _newclass:mom_comp = _swig_property(_cppcolleff.Ring_t_mom_comp_get, _cppcolleff.Ring_t_mom_comp_set)
    __swig_setmethods__["T0"] = _cppcolleff.Ring_t_T0_set
    __swig_getmethods__["T0"] = _cppcolleff.Ring_t_T0_get
    if _newclass:T0 = _swig_property(_cppcolleff.Ring_t_T0_get, _cppcolleff.Ring_t_T0_set)
    __swig_setmethods__["energy"] = _cppcolleff.Ring_t_energy_set
    __swig_getmethods__["energy"] = _cppcolleff.Ring_t_energy_get
    if _newclass:energy = _swig_property(_cppcolleff.Ring_t_energy_get, _cppcolleff.Ring_t_energy_set)
    __swig_setmethods__["cav_s"] = _cppcolleff.Ring_t_cav_s_set
    __swig_getmethods__["cav_s"] = _cppcolleff.Ring_t_cav_s_get
    if _newclass:cav_s = _swig_property(_cppcolleff.Ring_t_cav_s_get, _cppcolleff.Ring_t_cav_s_set)
    __swig_setmethods__["cav_V"] = _cppcolleff.Ring_t_cav_V_set
    __swig_getmethods__["cav_V"] = _cppcolleff.Ring_t_cav_V_get
    if _newclass:cav_V = _swig_property(_cppcolleff.Ring_t_cav_V_get, _cppcolleff.Ring_t_cav_V_set)
    def __init__(self): 
        this = _cppcolleff.new_Ring_t()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _cppcolleff.delete_Ring_t
    __del__ = lambda self : None;
    def track_one_turn(self, *args): return _cppcolleff.Ring_t_track_one_turn(self, *args)
Ring_t_swigregister = _cppcolleff.Ring_t_swigregister
Ring_t_swigregister(Ring_t)

class Results_t(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Results_t, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Results_t, name)
    __repr__ = _swig_repr
    __swig_setmethods__["nturns"] = _cppcolleff.Results_t_nturns_set
    __swig_getmethods__["nturns"] = _cppcolleff.Results_t_nturns_get
    if _newclass:nturns = _swig_property(_cppcolleff.Results_t_nturns_get, _cppcolleff.Results_t_nturns_set)
    __swig_setmethods__["xx_ave"] = _cppcolleff.Results_t_xx_ave_set
    __swig_getmethods__["xx_ave"] = _cppcolleff.Results_t_xx_ave_get
    if _newclass:xx_ave = _swig_property(_cppcolleff.Results_t_xx_ave_get, _cppcolleff.Results_t_xx_ave_set)
    __swig_setmethods__["xx_std"] = _cppcolleff.Results_t_xx_std_set
    __swig_getmethods__["xx_std"] = _cppcolleff.Results_t_xx_std_get
    if _newclass:xx_std = _swig_property(_cppcolleff.Results_t_xx_std_get, _cppcolleff.Results_t_xx_std_set)
    __swig_setmethods__["xl_ave"] = _cppcolleff.Results_t_xl_ave_set
    __swig_getmethods__["xl_ave"] = _cppcolleff.Results_t_xl_ave_get
    if _newclass:xl_ave = _swig_property(_cppcolleff.Results_t_xl_ave_get, _cppcolleff.Results_t_xl_ave_set)
    __swig_setmethods__["xl_std"] = _cppcolleff.Results_t_xl_std_set
    __swig_getmethods__["xl_std"] = _cppcolleff.Results_t_xl_std_get
    if _newclass:xl_std = _swig_property(_cppcolleff.Results_t_xl_std_get, _cppcolleff.Results_t_xl_std_set)
    __swig_setmethods__["de_ave"] = _cppcolleff.Results_t_de_ave_set
    __swig_getmethods__["de_ave"] = _cppcolleff.Results_t_de_ave_get
    if _newclass:de_ave = _swig_property(_cppcolleff.Results_t_de_ave_get, _cppcolleff.Results_t_de_ave_set)
    __swig_setmethods__["de_std"] = _cppcolleff.Results_t_de_std_set
    __swig_getmethods__["de_std"] = _cppcolleff.Results_t_de_std_get
    if _newclass:de_std = _swig_property(_cppcolleff.Results_t_de_std_get, _cppcolleff.Results_t_de_std_set)
    __swig_setmethods__["ss_ave"] = _cppcolleff.Results_t_ss_ave_set
    __swig_getmethods__["ss_ave"] = _cppcolleff.Results_t_ss_ave_get
    if _newclass:ss_ave = _swig_property(_cppcolleff.Results_t_ss_ave_get, _cppcolleff.Results_t_ss_ave_set)
    __swig_setmethods__["ss_std"] = _cppcolleff.Results_t_ss_std_set
    __swig_getmethods__["ss_std"] = _cppcolleff.Results_t_ss_std_get
    if _newclass:ss_std = _swig_property(_cppcolleff.Results_t_ss_std_get, _cppcolleff.Results_t_ss_std_set)
    __swig_setmethods__["Wlkick"] = _cppcolleff.Results_t_Wlkick_set
    __swig_getmethods__["Wlkick"] = _cppcolleff.Results_t_Wlkick_get
    if _newclass:Wlkick = _swig_property(_cppcolleff.Results_t_Wlkick_get, _cppcolleff.Results_t_Wlkick_set)
    __swig_setmethods__["Wdkick"] = _cppcolleff.Results_t_Wdkick_set
    __swig_getmethods__["Wdkick"] = _cppcolleff.Results_t_Wdkick_get
    if _newclass:Wdkick = _swig_property(_cppcolleff.Results_t_Wdkick_get, _cppcolleff.Results_t_Wdkick_set)
    __swig_setmethods__["FBkick"] = _cppcolleff.Results_t_FBkick_set
    __swig_getmethods__["FBkick"] = _cppcolleff.Results_t_FBkick_get
    if _newclass:FBkick = _swig_property(_cppcolleff.Results_t_FBkick_get, _cppcolleff.Results_t_FBkick_set)
    def __init__(self, *args): 
        this = _cppcolleff.new_Results_t(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _cppcolleff.delete_Results_t
    __del__ = lambda self : None;
    def calc_stats(self, *args): return _cppcolleff.Results_t_calc_stats(self, *args)
    def set_Wkicks(self, *args): return _cppcolleff.Results_t_set_Wkicks(self, *args)
    def set_FBkick(self, *args): return _cppcolleff.Results_t_set_FBkick(self, *args)
Results_t_swigregister = _cppcolleff.Results_t_swigregister
Results_t_swigregister(Results_t)

class Feedback_t(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Feedback_t, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Feedback_t, name)
    __repr__ = _swig_repr
    __swig_setmethods__["track"] = _cppcolleff.Feedback_t_track_set
    __swig_getmethods__["track"] = _cppcolleff.Feedback_t_track_get
    if _newclass:track = _swig_property(_cppcolleff.Feedback_t_track_get, _cppcolleff.Feedback_t_track_set)
    __swig_setmethods__["npoints"] = _cppcolleff.Feedback_t_npoints_set
    __swig_getmethods__["npoints"] = _cppcolleff.Feedback_t_npoints_get
    if _newclass:npoints = _swig_property(_cppcolleff.Feedback_t_npoints_get, _cppcolleff.Feedback_t_npoints_set)
    __swig_setmethods__["delay"] = _cppcolleff.Feedback_t_delay_set
    __swig_getmethods__["delay"] = _cppcolleff.Feedback_t_delay_get
    if _newclass:delay = _swig_property(_cppcolleff.Feedback_t_delay_get, _cppcolleff.Feedback_t_delay_set)
    __swig_setmethods__["phase"] = _cppcolleff.Feedback_t_phase_set
    __swig_getmethods__["phase"] = _cppcolleff.Feedback_t_phase_get
    if _newclass:phase = _swig_property(_cppcolleff.Feedback_t_phase_get, _cppcolleff.Feedback_t_phase_set)
    __swig_setmethods__["freq"] = _cppcolleff.Feedback_t_freq_set
    __swig_getmethods__["freq"] = _cppcolleff.Feedback_t_freq_get
    if _newclass:freq = _swig_property(_cppcolleff.Feedback_t_freq_get, _cppcolleff.Feedback_t_freq_set)
    __swig_setmethods__["gain"] = _cppcolleff.Feedback_t_gain_set
    __swig_getmethods__["gain"] = _cppcolleff.Feedback_t_gain_get
    if _newclass:gain = _swig_property(_cppcolleff.Feedback_t_gain_get, _cppcolleff.Feedback_t_gain_set)
    __swig_setmethods__["satur"] = _cppcolleff.Feedback_t_satur_set
    __swig_getmethods__["satur"] = _cppcolleff.Feedback_t_satur_get
    if _newclass:satur = _swig_property(_cppcolleff.Feedback_t_satur_get, _cppcolleff.Feedback_t_satur_set)
    __swig_setmethods__["bpmbeta"] = _cppcolleff.Feedback_t_bpmbeta_set
    __swig_getmethods__["bpmbeta"] = _cppcolleff.Feedback_t_bpmbeta_get
    if _newclass:bpmbeta = _swig_property(_cppcolleff.Feedback_t_bpmbeta_get, _cppcolleff.Feedback_t_bpmbeta_set)
    __swig_setmethods__["kikbeta"] = _cppcolleff.Feedback_t_kikbeta_set
    __swig_getmethods__["kikbeta"] = _cppcolleff.Feedback_t_kikbeta_get
    if _newclass:kikbeta = _swig_property(_cppcolleff.Feedback_t_kikbeta_get, _cppcolleff.Feedback_t_kikbeta_set)
    def __init__(self): 
        this = _cppcolleff.new_Feedback_t()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _cppcolleff.delete_Feedback_t
    __del__ = lambda self : None;
    def apply_kick(self, *args): return _cppcolleff.Feedback_t_apply_kick(self, *args)
Feedback_t_swigregister = _cppcolleff.Feedback_t_swigregister
Feedback_t_swigregister(Feedback_t)

class WakePl(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, WakePl, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, WakePl, name)
    __repr__ = _swig_repr
    __swig_setmethods__["general"] = _cppcolleff.WakePl_general_set
    __swig_getmethods__["general"] = _cppcolleff.WakePl_general_get
    if _newclass:general = _swig_property(_cppcolleff.WakePl_general_get, _cppcolleff.WakePl_general_set)
    __swig_setmethods__["resonator"] = _cppcolleff.WakePl_resonator_set
    __swig_getmethods__["resonator"] = _cppcolleff.WakePl_resonator_get
    if _newclass:resonator = _swig_property(_cppcolleff.WakePl_resonator_get, _cppcolleff.WakePl_resonator_set)
    __swig_setmethods__["s"] = _cppcolleff.WakePl_s_set
    __swig_getmethods__["s"] = _cppcolleff.WakePl_s_get
    if _newclass:s = _swig_property(_cppcolleff.WakePl_s_get, _cppcolleff.WakePl_s_set)
    __swig_setmethods__["W"] = _cppcolleff.WakePl_W_set
    __swig_getmethods__["W"] = _cppcolleff.WakePl_W_get
    if _newclass:W = _swig_property(_cppcolleff.WakePl_W_get, _cppcolleff.WakePl_W_set)
    __swig_setmethods__["wr"] = _cppcolleff.WakePl_wr_set
    __swig_getmethods__["wr"] = _cppcolleff.WakePl_wr_get
    if _newclass:wr = _swig_property(_cppcolleff.WakePl_wr_get, _cppcolleff.WakePl_wr_set)
    __swig_setmethods__["Rs"] = _cppcolleff.WakePl_Rs_set
    __swig_getmethods__["Rs"] = _cppcolleff.WakePl_Rs_get
    if _newclass:Rs = _swig_property(_cppcolleff.WakePl_Rs_get, _cppcolleff.WakePl_Rs_set)
    __swig_setmethods__["Q"] = _cppcolleff.WakePl_Q_set
    __swig_getmethods__["Q"] = _cppcolleff.WakePl_Q_get
    if _newclass:Q = _swig_property(_cppcolleff.WakePl_Q_get, _cppcolleff.WakePl_Q_set)
    def __init__(self): 
        this = _cppcolleff.new_WakePl()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _cppcolleff.delete_WakePl
    __del__ = lambda self : None;
WakePl_swigregister = _cppcolleff.WakePl_swigregister
WakePl_swigregister(WakePl)

class Wake_t(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Wake_t, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Wake_t, name)
    __repr__ = _swig_repr
    __swig_setmethods__["Wd"] = _cppcolleff.Wake_t_Wd_set
    __swig_getmethods__["Wd"] = _cppcolleff.Wake_t_Wd_get
    if _newclass:Wd = _swig_property(_cppcolleff.Wake_t_Wd_get, _cppcolleff.Wake_t_Wd_set)
    __swig_setmethods__["Wq"] = _cppcolleff.Wake_t_Wq_set
    __swig_getmethods__["Wq"] = _cppcolleff.Wake_t_Wq_get
    if _newclass:Wq = _swig_property(_cppcolleff.Wake_t_Wq_get, _cppcolleff.Wake_t_Wq_set)
    __swig_setmethods__["Wl"] = _cppcolleff.Wake_t_Wl_set
    __swig_getmethods__["Wl"] = _cppcolleff.Wake_t_Wl_get
    if _newclass:Wl = _swig_property(_cppcolleff.Wake_t_Wl_get, _cppcolleff.Wake_t_Wl_set)
    def __init__(self): 
        this = _cppcolleff.new_Wake_t()
        try: self.this.append(this)
        except: self.this = this
    def apply_kicks(self, *args): return _cppcolleff.Wake_t_apply_kicks(self, *args)
    __swig_destroy__ = _cppcolleff.delete_Wake_t
    __del__ = lambda self : None;
Wake_t_swigregister = _cppcolleff.Wake_t_swigregister
Wake_t_swigregister(Wake_t)


def do_tracking(*args):
  return _cppcolleff.do_tracking(*args)
do_tracking = _cppcolleff.do_tracking
# This file is compatible with both classic and new-style classes.


