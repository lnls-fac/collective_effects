function wakel = lnls_calc_wake_longitudinal_resonator(Rs, Q, wr, tau)
% Modelagem de imped??ncia por soma de ressonadores.
% Inputs:
% Rs    = Vetor com as resistencias Shunts [Ohm]
% Q     = Vetor com os fatores de qualidade 
% wr    = Vetor com as frequencias angulares de ressonancia de cada oscilador [rad/s]
% tau   = Vetor com o delay entre as particulas fonte e teste. DEVE SER
%         NEGATIVO. [s]
%
% Outputs: 
% wakel    = wake longitudinal calculado para tau
%
% Todos os outputs tem as dimensoes do SI
%
% Referencia orginal:
%
% Chao, A., Physics of Collective Beam Instabilities in High Energy
% Accelerators, Wiley 1993.

wakel = zeros(size(tau));
I = (tau ~= 0);
for i=1:length(Rs)
    Ql = sqrt(Q(i).^2 - 1/4);
    wrl = wr(i) .* Ql ./ Q(i);
    wakel(I) = wakel(I) - wr(i)*Rs(i)/Q(i)*(cos(wrl*tau(I)) +1/(2*Ql)*sin(wrl*tau(I))).*exp(wr(i)*tau(I)/(2*Q(i)));
    wakel(~I) = wakel(~I) - (1/2) * wr(i)*Rs(i)/Q(i);
end
