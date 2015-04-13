function [Zv Zh Zl w_val varargout] = lnls_calc_impedance_flat_chamber(thick, cond, perm, h, L, y0, w)
% [Zv Zh Zl w_val] = lnls_calc_impedance_flat_chamber(thick, cond, perm, h, L, y0, w)
% [Zv Zh Zl w_val Zqv Zqh] = lnls_calc_impedance_flat_chamber(thick, cond, perm, h, L, y0, w)
% Imped???ncia de duas placas paralelas.
% Inputs:
% thick = espessura das placas [m]
% cond  = condutividade eletrica do material constituinte [1/(Ohm.m)]
% perm  = permeabilidade magnetica relativa do material
% h     = distancia entre as placas [m]
% L     = comprimento longitudinal das placas [m]
% y0    = posicao do feixe relativa ao centro das duas placas [m]
% w     = frequencias para as quais a impedancia sera calculada [Hz]
%
% Outputs: 
% Zv    = impedancia dipolar vertical calculada para w 
% Zh    = impedancia dipolar vertical calculada para w 
% Zqv   = impedancia quadrupolar vertical calculada para w 
% Zqh   = impedancia quadrupolar vertical calculada para w 
% Zl    = impedancia longitudinal calculada para w 
% w_val = indices das frequwncias para as quais a formula usada para o
%         calculo da impedancia e valida (leva em consideracao apenas o
%         limite de aplicacao 2, descrito abaixo).
%
% Todos os outputs tem as dimensoes do SI
%
% Limites de aplicacao:
% 1-) gamma >> Max(1,h^2/sigma_s^2)
% gamma   = energia relativistica da particula
% sigma_s = tamanho longitudinal do pacote
%
% 2-) Min(1/h/k^2, h/2-y0, thick) >> delta_c
% delta_c = skin depth
% k       = w/c, numero de onda
% 
% A primeira implementacao foi baseada no artigo: 
% A. Piwinski, Impedances in Lossy Elliptical Vacuum Chambers DESY 94-068,
% April 1994. Encontrado em http://www-library.desy.de/report94.html
%
% Em 2012/12/12 descobri a referencia
% Gunzel, R.F., Transverse Coupling Impedance of the storage ring at the
% European Synchrotron Radiation Facility, Phys. Rev. Special topics, 9,
% 114402 (2006)
%
% que me alertou a existencia do wake quadrupolar para camaras sem simetria
% axial, mas com simetria nos dois planos e citou a referencia
% que possuia a derivacao da impedancia para as flat chambers.
%
%Por enquanto nao estou estudando o efeito de o ondulador ter uma espessura
%muito fina de metal. preciso estudar melhor esse tipo de situacao.

% Constantes utilizadas
c   = 299792458; % velocidade da luz
u0  = 4e-7*pi; % permeabilidade magnetica do vacuo
Z0  = c*u0; % Impedancia do vacuo
bigger = 3; %constante multiplicativa de delta_c para indicar ">>"

% % Funcoes auxiliares
% fz = 1+pi*y0/h*tan(pi*y0/h);
% ft = fz/h^3/cos(pi*y0/h)^2;
% 
% w_val = [];
% Zt    = zeros(1,length(w));
% Zl    = zeros(1,length(w));
% for j=1:length(w)
%     k = w(j)/c;
%     delta_c = sqrt(2/Z0/perm/cond/abs(k));
%     if (min([1/h/k^2, h/2-y0, thick])> bigger*delta_c)
%         w_val(end+1) = j;
%     end
%     Zl(j)  =  L*(1-1i*sign(w(j)))/pi/h*sqrt(abs(w(j))*perm*Z0/2/c/cond)*fz;
%     Zt(j)  =  L*pi*(sign(w(j)) - 1i)/sqrt(2*abs(w(j))*cond/c/perm/Z0)*ft;
% end



k = w/c;
delta_c = sqrt(2/Z0/perm/cond./abs(k));

w_val = (min([1/h./k.^2, h/2, thick])> bigger*delta_c);

Zl  =  L*(1-1i*sign(w))/pi/h.*sqrt(abs(w)*perm*Z0/2/c/cond);
Zt  =   L * (pi^2/12) * (sign(w) - 1i) * (Z0/(2*pi*(h/2)^3)).*delta_c;

Zv = Zt;   Zqv = Zv/2;
Zh = Zt/2; Zqh = Zh;

varargout{1} = Zqv;
varargout{2} = Zqh;



