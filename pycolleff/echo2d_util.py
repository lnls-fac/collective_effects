#!/usr/bin/env python3

import numpy as _np
import matplotlib.pyplot as _plt

class EchoObj(_np.ndarray):
    def __new__(cls):
        return _np.ndarray.__new__(cls,10,dtype=float)
    def __init__(self):
        for i in range(len(self)):
            self[i] = 0.0

def line(x1,y1,x2,y2):
    t = EchoObj()
    t[0] = x1
    t[1] = y1
    t[2] = x2
    t[3] = y2
    return t

def circle_out(x,y,R,theta):
    "(x,y) is the position there the tangent lines cross."
    t = EchoObj()
    t[0] = x - R*_np.tan(theta/2)
    t[1] = y
    t[2] = x + R*_np.tan(theta/2)*_np.cos(theta)
    t[3] = y - R*_np.tan(theta/2)*_np.sin(theta)
    t[4] = t[0] - R
    t[5] = y
    t[6] = t[0] + R
    t[7] = y - 2*R
    t[8] = 0
    return t

def circle_inn(x,y,R,theta):
    "(x,y) is the position there the tangent lines cross."
    t = EchoObj()
    t[0] = x - R*_np.tan(theta/2)*_np.cos(theta)
    t[1] = y + R*_np.tan(theta/2)*_np.sin(theta)
    t[2] = x + R*_np.tan(theta/2)
    t[3] = y
    t[4] = t[2] - R
    t[5] = y + 2*R
    t[6] = t[2] + R
    t[7] = y
    t[8] = 1
    return t

def reflect(t_in):
    if isinstance(t_in,(list,tuple)) and isinstance(t_in[0],EchoObj):
        t_out = t_in.copy()
        for i in range(len(t_out)): t_out[-i-1] = reflect(t_in[i])
    elif isinstance(t_in,EchoObj):
        t_out = t_in.copy()
        t_out[0] = -t_in[2]
        t_out[1] =  t_in[3]
        t_out[2] = -t_in[0]
        t_out[3] =  t_in[1]
        t_out[4] = -t_in[6]
        t_out[6] = -t_in[4]
    else:
        raise Exception('error')
    return t_out

def translate(t_in,delta):
    if isinstance(t_in,(list,tuple)) and isinstance(t_in[0],EchoObj):
        t_out = t_in.copy()
        for i in range(len(t_in)): t_out[i] = translate(t_in[i],delta)
    elif isinstance(t_in,EchoObj):
        t_out = t_in.copy()
        t_out[[0,2]] += delta
        if _np.any(t_out[[4,6]] != 0): t_out[[4,6]] += delta
    else:
        raise Exception('error')
    return t_out

def create_linear_taper(fname=None, r_in=0.012, r_out=0.004, t = 20,
                         s_in=0.0025, s_out=0.0025, C=None):
    p2 = _np.array([ 0.0           , r_in])
    p3 = _np.array([t*_np.abs(r_in-r_out) , r_out])

    if C is None:
        t1 = line(p2[0],p2[1],p3[0],p3[1])
        points = [t1,]
    else:
        theta  = _np.arctan(1/t)
        if r_in > r_out:
            t1 = circle_out(p2[0],p2[1],r_out,theta)
            t3 = circle_inn(p3[0],p3[1],r_in,theta)
        else:
            t1 = circle_inn(p2[0],p2[1],r_out,theta)
            t3 = circle_out(p3[0],p3[1],r_in,theta)
        t2 = line(t1[2],t1[3],t3[0],t3[1])
        points = [t1,t2,t3]

    if s_in > 0.0:
        t0 = line(t1[0]-s_in,t1[1],t1[0],t1[1])
        points = [t0,] + points
    if s_out > 0.0:
        t_end = points[-1]
        t_end = line(t_end[2],t_end[3],t_end[2]+s_out,t_end[3])
        points += [t_end,]

    if fname is not None:
        create_geometry_file(fname,points)

    return points

def create_collimator(fname=None,R_in=0.012,R_out=None, r=0.004,
                                 t_in=20,   t_out=None, g=0.02, C_in=None,C_out=None):

    if t_out is None: t_out = t_in
    if R_out is None: R_out = R_in
    init  = 0.0025

    theta_in  = _np.arctan(1/t_in)
    theta_out = _np.arctan(1/t_out)

    points  = create_linear_taper(r_in=R_in, r_out=r,t=t_in, s_in=init,s_out=0.0,C=C_in)
    print(points)
    t_m     = line(points[-1][2],points[-1][3],points[-1][2]+g,points[-1][3])
    print(t_m)
    points2 = create_linear_taper(r_in=r,r_out=R_out,t=t_out,s_in=0.0,s_out=init,C=C_out)
    print(points2)
    points2 = translate(points2,t_m[2]-points2[0][0])
    print(points2)
    points  += [t_m,] + points2

    if fname is not None:
        create_geometry_file(fname,points)

    return points

def create_geometry_file(fname):
    if not fname.endswith('.txt'): fname += '.txt'
    with open(fname,'w') as f:
        f.write('{0:d}'.format(len(points)) + '\n')
        for p in points:  f.write(''.join(['{0:12.5f}'.format(x) for x in p]) + '\n')

def plot_geometry(points):
    _plt.figure()
    for p in points:
        color = 'b' if _np.all(p[[4,5,6,7]]==0) else 'r'
        _plt.plot(p[[0,2]],p[[1,3]]*1e3,color)
    _plt.xlabel('s [m]')
    _plt.ylabel('R [mm]')
    _plt.grid('on')
    _plt.show()

def create_input_file(fname, geo_fname,
    geo_unit='m', geo_type='recta',geo_width=0.024,geo_bound='magn',geo_conv=True,
    beam_sigma=5e-4, beam_offset=-1,
    modes = (0,1),
    mesh_leng = 10001, mesh_step=None):

    if mesh_step is None: mesh_step = beam_sigma/5.0

    with open(fname, 'w') as f:
        f.write('%%%% Generated automatically\n')
        f.write('''
%%%%%%%%%%%%%%%%%%%%% geometry %%%%%%%%%%%%%%%%%%%%%%%

GeometryFile= {0:<21s} % geometry in "Units"
Units= {1:<28s} % m/cm/mm - only for GeometryFile!!!
GeometryType= {2:<21s} % recta/round
Width= {3:<28.5g} % in m
SymmetryCondition= {4:<16s} % magn/elec - boundary condition on axis
Convex= {5:<27d} % 1(convex)/0(no)
'''.format(geo_fname,geo_unit,geo_type,geo_width,geo_bound,1 if geo_conv else 0))


        f.write('''
%%%%%%%%%%%%%%% input beam and field %%%%%%%%%%%%%%%%%%

InPartFile= -                       % -(Gaussian pencil beam)
BunchSigma= {0:<23.5g} % in m
Offset= {1:<27d} % in mesh steps/-1(default, near wall)
InFieldFile= -                      %  -(no output)
'''.format(beam_sigma,beam_offset))


        f.write('''
%%%%%%%%%%%%%%%%%%%%%% model %%%%%%%%%%%%%%%%%%%%%%%%%%%%

WakeIntMethod= ind                  % ind/dir/-
Modes= {0:s}
ParticleMotion= -                   % 1d/2d/3d/-
'''.format(' '.join(['{0:d}'.format(x) for x in modes])))


        f.write('''
%%%%%%%%%%%%%%%%%%%%%% mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%

MeshLength= {0:<23d} % in mesh steps
StepY= {1:<28.5g} %
StepZ= {1:<28.5g} %
NStepsInConductive=0                % 0(default)
TimeSteps=0                         % 0(default)
'''.format(mesh_leng,mesh_step))

        f.write('''
%%%%%%%%%%%%%%%%%%%% monitors %%%%%%%%%%%%%%%%%%%%%%%%%%

OutPartFile= -                      % -(no output)
OutFieldFile= -                     % -(no output)
''')

# def parse_from_input_file(fname):
#
#     if mesh_step is None: mesh_step = beam_sigma/5.0
#
# geo_fname,
#     geo_unit='m', geo_type='recta',geo_width=0.024,geo_bound='magn',geo_conv=True,
#     beam_sigma=5e-4, beam_offset=-1,
#     modes = (0,1),
#     mesh_leng = 10001, mesh_step=None
#
#     with open(fname, 'r') as f:
#         f.write('%%%% Generated automatically\n')
#         f.write('''
# %%%%%%%%%%%%%%%%%%%%% geometry %%%%%%%%%%%%%%%%%%%%%%%
#
# GeometryFile= {0:<21s} % geometry in "Units"
# Units= {1:<28s} % m/cm/mm - only for GeometryFile!!!
# GeometryType= {2:<21s} % recta/round
# Width= {3:<28.5g} % in m
# SymmetryCondition= {4:<16s} % magn/elec - boundary condition on axis
# Convex= {5:<27d} % 1(convex)/0(no)
# '''.format(geo_fname,geo_unit,geo_type,geo_width,geo_bound,1 if geo_conv else 0))
#
#
#         f.write('''
# %%%%%%%%%%%%%%% input beam and field %%%%%%%%%%%%%%%%%%
#
# InPartFile= -                       % -(Gaussian pencil beam)
# BunchSigma= {0:<23.5g} % in m
# Offset= {1:<27d} % in mesh steps/-1(default, near wall)
# InFieldFile= -                      %  -(no output)
# '''.format(beam_sigma,beam_offset))
#
#
#         f.write('''
# %%%%%%%%%%%%%%%%%%%%%% model %%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# WakeIntMethod= ind                  % ind/dir/-
# Modes= {0:s}
# ParticleMotion= -                   % 1d/2d/3d/-
# '''.format(' '.join(['{0:d}'.format(x) for x in modes])))
#
#
#         f.write('''
# %%%%%%%%%%%%%%%%%%%%%% mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
# MeshLength= {0:<23d} % in mesh steps
# StepY= {1:<28.5g} %
# StepZ= {1:<28.5g} %
# NStepsInConductive=0                % 0(default)
# TimeSteps=0                         % 0(default)
# '''.format(mesh_leng,mesh_step))
#
#         f.write('''
# %%%%%%%%%%%%%%%%%%%% monitors %%%%%%%%%%%%%%%%%%%%%%%%%%
#
# OutPartFile= -                      % -(no output)
# OutFieldFile= -                     % -(no output)
# ''')
