function Zl = lnls_calc_impedance_longitudinal_resonator(Rs, Q, wr, w)
% Modelagem de impedância por soma de ressonadores.
% Inputs:
% Rs    = Vetor com as resistências Shunts [Ohm]
% Q     = Vetor com os fatores de qualidade 
% wr    = Vetor com as frequências angulares de ressonância de cada oscilador [rad/s]
% w     = Frequências angulares em que sera calculada a impedância [rad/s]
%
% Outputs: 
% Zl    = impedância longitudinal (real e imaginária) calculada para w
%
% Todos os outputs tem as dimens�es do SI
%
% Referência orginal:
%
% Chao, A., Physics of Collective Beam Instabilities in High Energy
% Accelerators, Wiley 1993.

Zl = zeros(1,length(w));
for i=1:length(Rs)
    Zl    = Zl + Rs(i)./(1+1i*Q(i)*(wr(i)./w-w/wr(i)));
end