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

% Percebi que a implementa????o padr??o ?? tal que Ztmax = Rs. Ent??o
% multipliquei a implmenta????o inicial (do Chao) por wr/c. (2012/12/11)
% Separei o c??lculo de frequ??ncia =0 para para n??o divergir.

% Constantes utilizadas
c   = 299792458; % velocidade da luz

Zt = zeros(1,length(w));
ind = w ~= 0;
for i=1:length(Rs)
    Zt(ind)    = Zt(ind) + (wr(i)/c)*(c./w(ind))*Rs(i)./(1+1i*Q(i)*(wr(i)./w(ind)-w(ind)/wr(i)));
    Zt(~ind)   = Zt(~ind) + (wr(i)/c)*c*Rs(i)/(1i*Q(i)*wr(i));
end
