function [lossf, Zl_eff] = lnls_calc_loss_factor(w,Z,sigma,w0,nb)
% Calcula o loss factor and effective impedance para nb pacotes com 
% comprimento longitudinal sigma igualmente espacados.
%
% Chamada:
%   lossf = lnls_calc_loss_factor(w,Z,sigma,w0,nb)
%
% Inputs:
%   w = frequencia angular [rad/s]
%   Z = Impedancia longitudinal [Ohm]
%   sigma = tamanho longitudinal do feixe [m]
%   w0 = frequencia angular de revolucao no anel [rad/s]
%   nb = numero de pacotes preenchidos

c = 299792458;

pmin = ceil((w(1)+nb*w0)/(w0*nb)); % arredonda em direcao a +infinito
pmax = ceil((w(end)-nb*w0)/(w0*nb))-1; % arredonda em direcao a -infinito

p = pmin:pmax;
wp = w0*p*nb;

h = exp(-(wp*sigma/c).^2);
interpol_Z = interp1(w,real(Z),wp);

lossf = nb*(w0/2/pi)*sum(interpol_Z.*h);
interpol_Z = interp1(w,imag(Z),wp);
Zl_eff = nb*w0*sum(interpol_Z.*h./(wp+0.001)); % para evitar 0/0