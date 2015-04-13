function [kickf, Zt_eff] = lnls_calc_kick_factor(w,Z,sigma,w0,nb)

% c = 299792458;
% 
% h = exp(-(w*sigma/c).^2);
% ImZ = imag(Z);
% 
% kickf = (1/(2*pi))*trapz(w,ImZ.*h);


% Calcula o kick factor para nb pacotes com comprimento longitudinal sigma
% igualmente espacados.

c = 299792458;

pmin = ceil((w(1)+nb*w0)/(w0*nb)); % arredonda em direcao a +infinito
pmax = ceil((w(end)-nb*w0)/(w0*nb))-1; % arredonda em direcao a -infinito

p = pmin:pmax;
wp = w0*p*nb;

h = exp(-(wp*sigma/c).^2);
interpol_Z = interp1(w,imag(Z),wp);

Zt_eff = sum(interpol_Z.*h);
kickf = nb*(w0/2/pi)*Zt_eff;
