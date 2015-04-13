function [Zl, Zv, Zh, Zl2] = lnls_calc_impedance_bellow(w,raio,periodo,corruga, comprimento)
c = 299792458;
mu0 = 4*pi*1e-7;
Z0 = c*mu0;

k = w/c;

% R = 17.2e-3;
% L = 1.156e-3;
% t = 6.8e-3;

R = raio;
L = periodo;
t = corruga;

epsilon = t/(2*R);
eta = 2*pi*R/L;
kappa = R*k;

y = epsilon*eta;
fy = 2*y*besseli(1,y)/(besseli(0,y)+besseli(2,y));
B1 = epsilon/eta*kappa*fy;

% Zv = -1i*comprimento*Z0*2*L*B1./(4*pi*R^3*k);
Zl2 = -1i*comprimento*Z0*2*L*B1./(4*pi*R);

Zv = -comprimento*1i*Z0*L*t/(4*pi*R^3)*ones(size(k));
Zl = -comprimento*1i*Z0*L*t/(4*pi*R)*k;
Zh = Zv;


