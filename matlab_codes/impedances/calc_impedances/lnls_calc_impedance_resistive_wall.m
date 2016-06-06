function [Zt Zl w_val] = lnls_calc_impedance_resistive_wall(thick, cond, perm, b, L, w)
% Impedância de parede resistiva para câmara de vácuo cil�ndrica.
% Inputs:
% thick = espessura da parede [m]
% cond  = condutividade elétrica do material constituinte [1/(Ohm.m)]
% perm  = permeabilidade magn�tica relativa do material
% b     = raio da c�mara [m]
% L     = comprimento longitudinal da c�mara [m]
% y0    = posição do feixe relativa ao centro da câmara [m]
% w     = frequências para as quais a impedância será calculada [Hz]
%
% Outputs: 
% Zt    = impedância transversal (real e imaginária) calculada para w 
% Zl    = impedância longitudinal (real e imaginária) calculada para w 
% w_val = �ndices das frequências para as quais a fórmula usada para o
%         c�lculo da impedância é válida (leva em consideração apenas os
%         limites de aplicação 1 e 2, descritos abaixo).
%
% Todos os outputs tem as dimens�es do SI
%
% Limites de aplica��o:
% 1-) delta_c << Min(b,thick)
% delta_c = skin depth
%
% 2-) |k| << (b/s0)^3/b
% s0 = comprimento caracter�stico (definido abaixo) 
%
% 3-) Considerei a condutividade DC, pois segundo BANE, para um tubo de
% cobre de 10 mm de raio, a diferença entre as duas (AC e DC) só seria 
% significativa para frequências da ordem de 1e13 Hz.
%
% Referência orginal:
%
% Chao, A., Physics of Collective Beam Instabilities in High Energy
% Accelerators, Wiley 1993.
%
% Referência Alternativa (contém generalizações da fórmula):
% Bane, K.; Sands, M., The Short-Range Resistive Wall Wakefields, SLAC-PUB
% 7074 1995, encontrado em:
% http://www.slac.stanford.edu/cgi-wrap/getdoc/slac-pub-7074.pdf
%
% Referência da qual retirei as fórmulas (já estava no SI):
% Ng, K. Y., Explicit Expressions of Impedances and Wake Functions,
% SLAC-PUB-15078 2010
%

% Constantes utilizadas
c   = 299792458; % velocidade da luz
u0  = 4e-7*pi; % permeabilidade magnética do vácuo
Z0  = c*u0; % Impedância do vácuo
bigger = 3; %constante multiplicativa de delta_c para indicar ">>"

% Parâmetros auxiliares
s0 = (2*b^2*perm/Z0/cond)^(1/3); % ~ 20um

w_val = [];
Zt    = zeros(1,length(w));
Zl    = zeros(1,length(w));
for j=1:length(w)
    k = w(j)/c;
    delta_c = sqrt(2/Z0/perm/cond/abs(k));
    if ((min([b, thick])> bigger*delta_c)&&(abs(k)<bigger*(b/s0)^3/b))
        w_val(end+1) = j;
    end
    Zl(j)  =  L*Z0/pi/b/(2*sqrt(1i*Z0*cond/k/perm) - 1i*b*k);
    Zt(j)  =  (1/k)*L*Z0/pi/b^3/(sqrt(1i*Z0*cond/k/perm) - 1i*b*k/2);
end



%% Fórmulas que consideram espessura finita da câmara (rever).
% for j=1:length(w)
%     k = w(j)/c;
%     delta_c = sqrt(2/Z0/perm/cond/abs(k));
%     if ((b> bigger*abs(delta_c)/2)&&(abs(k)<bigger*(b/s0)^3/b))
%         w_val(end+1) = j;
%     end
%     Zl(j)  =  L*Z0/pi/b/(2*sqrt(1i*Z0*cond/k/perm) - 1i*b*k)* ...
%         (1-exp(-2*(1-sign(k)*1i)*thick/delta_c))/(1+exp(-2*(1-sign(k)*1i)*thick/delta_c));
%     Zt(j)  =  (1/k)*L*Z0/pi/b^3/(sqrt(1i*Z0*cond/k/perm) - 1i*b*k/2)* ...
%         (1-exp(-2*(1-sign(k)*1i)*thick/delta_c))/(1+exp(-2*(1-sign(k)*1i)*thick/delta_c));
% end

