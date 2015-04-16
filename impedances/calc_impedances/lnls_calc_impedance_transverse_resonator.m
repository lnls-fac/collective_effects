function Zt = lnls_calc_impedance_transverse_resonator(Rs, Q, wr, w)
% Modelagem de imped??ncia por soma de ressonadores.
% Inputs:
% Rs    = Vetor com as resist??ncias Shunts [Ohm]
% Q     = Vetor com os fatores de qualidade 
% wr    = Vetor com as frequ??ncias angulares de resson??ncia de cada oscilador [rad/s]
% w     = Frequ??ncias angulares em que sera calculada a imped??ncia [rad/s]
%
% Outputs: 
% Zt    = imped??ncia transversal (real e imagin??ria) calculada para w
%
% Todos os outputs tem as dimens???es do SI
%
% Refer???ncia orginal:
%
% Chao, A., Physics of Collective Beam Instabilities in High Energy
% Accelerators, Wiley 1993.

Zt = zeros(1,length(w));
for i=1:length(Rs)
    Zt    = Zt + wr(i)*Rs(i)./(w + 1i*Q(i)*(wr(i)-w.^2/wr(i)));
end
