function Zl = lnls_calc_metallic_coating_on_ceramic_pipe(tc, tm, cond, perm, b, L, sigmal, w)
% FALTA TERMINAR DE IMPLEMENTAR: APARENTE DISCORDÂNCIA ENTRE BIBLIOGRAFIAS.
% Impedância tubo cilíndrico de cerâmica com coating metálico.
% Inputs:
% tc = espessura da parede de cerâmica [m]
% tc = espessura do coating metálico [m]
% cond  = condutividade elétrica do metal [1/(Ohm.m)]
% perm  = permeabilidade elétrica relativa da cerâmica
% b     = raio da câmara [m]
% L     = comprimento longitudinal da câmara [m]
% sigmal = comprimento longitudinal do pacote [m]
% w     = frequências para as quais a impedância será calculada [Hz]
%
% Outputs: 
% Zl    = impedância longitudinal (real e imaginária) calculada para w 
%
% Todos os outputs tem as dimensões do SI
%
% Limites de aplicação:
% 1-) b >> Max(tc,tm)
%
% 2-) Max((perm-1)*tc^2/3,(1-1/perm)*b*tc/2) << sigmal^2
%
% Considerações: O problema é resolvido pelo método exato, de casamento dos
% campos nas interfaces, em [1], portanto, as aproximações são oriundas
% dessa referência. Contudo, a fórmula para a impedância só é apresentada
% em [2], onde o problema é resolvido perturbativamente e não há
% especificação clara das aproximações. Não consegui chegar no resultado de
% [2] a partir do de [1], é necessário estudar mais o problema.
% A fórmula empregada foi retirada de [3], que apenas reescreve e
% particulariza para o caso de metal-cerâmica o resultado de [2], que é
% geral e pode ser aplicado para metal-metal também.
% Ainda, [4] modela o problema diferentemente, considerando um comprimento 
% finito para o sistema. 
%
% Referências:
%
% [1] Piwinski A.; Penetration of the field of a bunched beam through a 
%     ceramic vacuum chamber with metallic coating, IEEE 1977.
% [2] Lin, X. E.; RF loss in and leakage through thin metal film,
%     SLAC-PUB-7924 1998.
% [3] Ng, K. Y., Explicit Expressions of Impedances and Wake Functions,
%     SLAC-PUB-15078 2010.
% [4] Danilov V. et al; An Improved Impedance Model of Metallic Coatings,
%     EPAC 2002.
%

% Constantes utilizadas
c   = 299792458; % velocidade da luz
u0  = 4e-7*pi; % permeabilidade magnética do vácuo
Z0  = c*u0; % Impedância do vácuo
bigger = 5; %constante multiplicativa de delta_c para indicar ">>"

% Parâmetros auxiliares
s0 = (2*b^2*perm/Z0/cond)^(1/3); % ~ 20um

x = (b > bigger*max(tc,tm)) && ...
    (max([(perm-1)*tc^2/3,(1-1/perm)*b*tc/2])*bigger < sigmal^2);
if x
    k       = w/c;
    delta_c = sqrt(2/Z0/perm/cond/abs(k));
    nu      = (1-sign(w)*1i)./delta_c;
    A       = (1-1/perm)*nu*tc;
    Z0m     = 
    Zl      = L*Z0m*(A+tanh(nu*tm))/(1+A*tanh(nu*tm));
else
    disp('Expressão não é válida');
end


