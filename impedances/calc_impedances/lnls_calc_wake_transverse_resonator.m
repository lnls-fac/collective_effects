function waket = lnls_calc_wake_transverse_resonator(Rs, Q, wr, tau)
% Modelagem de imped??ncia por soma de ressonadores.
% Inputs:
% Rs    = Vetor com as resist??ncias Shunts [Ohm]
% Q     = Vetor com os fatores de qualidade 
% wr    = Vetor com as frequ??ncias angulares de resson??ncia de cada oscilador [rad/s]
% tau   = Vetor com o delay entre as particulas fonte e teste. DEVE SER
%         NEGATIVO. [s]
%
% Outputs: 
% waket    = wake transversal calculado para tau
%
% Todos os outputs tem as dimensoes do SI
%
% Referencia orginal:
%
% Chao, A., Physics of Collective Beam Instabilities in High Energy
% Accelerators, Wiley 1993.

waket = zeros(size(tau));
for i=1:length(Rs)
    Ql = sqrt(Q(i).^2 - 1/4);
    wrl = wr(i) .* Ql ./ Q(i);
    waket = waket + wr(i)*Rs(i)/Ql*sin(wrl*tau).*exp(wr(i)*tau/(2*Q(i)));
end
