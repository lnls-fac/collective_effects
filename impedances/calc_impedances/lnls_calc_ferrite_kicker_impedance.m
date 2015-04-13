function [Zl, Zh, Zv] = lnls_calc_ferrite_kicker_impedance(w,a,b,d,L,epr,mur,Zg,model, coupled)
% Calculates Impedances for a ferrite kicker: 
%   - For the Coupled Flux, it uses Davino-Hahn model.
%
%   DAVINO-HAHN MODEL:
%  
%  #######################################    |
%  ###############FERRITE#################    t
%  #######################################    |
%  ###**                             **###  |  
%  ###**  VACUUM                     **###  |      ______|  |_________|
%  ###**                             **###  |            |  |         |
%  ###**             +   .           **###  w            |  |         |
%  ###**                             **###  |            )||(         \
%  ###**             |_D_|           **###  |      Zk  L1)||(L2     Zg/
%  ###**                             **###  |            )||(         \
%  #######################################               | M|         /
%  #######################################               |  |         |     
%  #######################################         ______|  |_________|
%      |______________h______________|
%
%   - For the Uncoupled Flux, we can choose between three models:
%
%       TSUTSUI MODEL:
%  
%  ******************PEC**********************
%  *******************************************
%  **#######################################**      |
%  **################FERRITE################**      |
%  **#######################################**      |
%  **                                       **  |   d
%  **                                       **  b   |
%  **                                       **  |   |
%  **     VACUUM        .                   **  |   |
%  **                                       **
%  **                                       **
%  **                                       **
%  **#######################################**
%  **#######################################**
%  **#######################################**
%  *******************************************
%  *******************************************
%                       |__________a________| 
%
%       PIOR CASO: 
%       multi-layer cilindrica com vacuo-coating(2microns)-ceramica-ferrite-PEC
%
%       MELHOR CASO:
%       multi-layer cilindrica com vacuo-coating(2microns)-ceramica-PEC
%
% Inputs:
%
% w   = vector of angular frequencies to evaluate impedances [rad/s];
% epr = vector with real and imaginary electric permeability of ferrite for
%       the frequency range of interest;
% mur = vector with real and imaginary magnetic permissivity of ferrite for
%       the frequency range of interest;
% L   = length of the structure [m];
% Zg  = Vector with generator impedance in function of w [Ohm]
%
% Outputs:
%
% Zl = Impedancia Longitudinal [Ohm]
% Zh = Impedancia Horizontal [Ohm/m]
% Zv = Impedancia Vertical   [Ohm/m]
%
% Bibliografias:
%
% - Tsutsui_H - Some Simplified Models of Ferrite Kicker Magnet for 
%   Calculation of longitudinal Coupling Impedance - CERN-SL-2000-004
%
% - Tsutsui_H - Transverse Coupling Impedance of a Simplified Ferrite
%   Kicker Magnet Model - LHC Project Note 234 - 2000
%
% - Davino_D Hahn_H - Improved Analytical Model of the transverse coupling
%   impedance of ferrite kicker magnets - Phys. Rev. ST-AB v6 012001 2003
%
% - Nassibian G Sacherer F - Methods for measuring tranverse coupling 
%   impedances in circular Accelerators - Nucl Inst and Meth. 159 21-27 1979

c   = 299792458;
mu0 = 4*pi*1e-7;
ep0 = 1/c^2/mu0;


epb     = [1 1                9.3            1      12   1];
mub     = [1 1                 1             1       1   1];
ange    = [0 0                 0             0       0   0];
angm    = [0 0                 0             0       0   0];
sigmadc = [0 5.9e7             1             1       1  5.9e7];
tau     = [0 0                 0             0       0   0]*27e-15;
b1       = [(b - 2.0e-3 - 2e-6), (b - 2.0e-3), (b-1.0e-3), b , d];

for j = 1: length(epb)
    epr1(j,:) = epb(j)*(1-1i.*sign(w).*tan(ange(j))) + sigmadc(j)./(1+1i*w*tau(j))./(1i*w*ep0);
    mur1(j,:) = mub(j)*(1-1i.*sign(w).*tan(angm(j)));
end
epr1(5,:) = epr;
mur1(5,:) = mur;

if strcmp(model, 'tsutsui')
    [Zl, Zh, Zv] = lnls_calc_impedance_tsutsui_model(w, epr, mur, a, b, d, L, 10);
elseif strcmp(model,'pior')
    [Zl, Zv, Zh] = lnls_calc_impedance_multilayer_round_pipe(w, epr1, mur1, b1, L, 3);
    Zv = pi^2/12*Zv;
    Zh = pi^2/24*Zh;
else
    indx = [1 2 3 4 6];
    mur1 = mur1(indx,:);
    epr1 = epr1(indx,:);
    b1    = b1([1 2 3 4]);
    [Zl, Zv, Zh] = lnls_calc_impedance_multilayer_round_pipe(w, epr1, mur1, b1, L, 3);
    Zv = pi^2/12*Zv;
    Zh = pi^2/24*Zh;
end

if ~exist('coupled','var'), coupled = true; end;

if coupled
    % Equivalent Circuit model.
    h = 2*a;
    W = 2*b;
    t = d-b;
    D = 0.5e-3;
    M  = L*D*mu0/W;
    %     L2 = L*2*a*mu0/2/b;
    L2 = L*h*mu0/W*(mur*t./(mur*t+h*(h/W+1)));
    
    Zk = conj((M./L2).^2 .* Zg.*L2*1i.*w./(1i*w.*L2 + Zg));
    Zx = c./w/D^2 .* Zk;
    
    Zl = Zl + Zk;
    Zh = Zh + Zx;
end

    

function [Zl, Zh, Zv] = lnls_calc_impedance_tsutsui_model(w, epr, mur, a, b, d, L, n)
% Inputs:
%
% w   = vector of angular frequencies to evaluate impedances [rad/s];
% epr = vector with real and imaginary electric permeability of ferrite for
%       the frequency range of interest;
% mur = vector with real and imaginary magnetic permissivity of ferrite for
%       the frequency range of interest;
% n   = numbers of terms to sum;
% L   = length of the structure [m];
%
% Outputs:
%
% Zl = Impedancia Longitudinal [Ohm]
% Zh = Impedancia Horizontal [Ohm/m]
% Zv = Impedancia Vertical   [Ohm/m]

% Valores do artigo do Wang et al para testar a implementacao das formulas
% do modelo do Tsutui.
% a = 103e-3;
% b = 67e-3;
% d = 122e-3;
% L = 1.0;
% Valores do artigo do Tsutsui para testar a implementacao das formulas
% a = 20e-3;
% b = 16e-3;
% d = 66e-3;
% L = 0.6;

c   = 299792458;
mu0 = 4*pi*1e-7;
Z0  = mu0*c;

% Terms for the infinit sum:
n = (0:n)';

k = ones(size(n))*w/c;

epr = ones(size(n))*epr;
mur = ones(size(n))*mur;

kxn = (2*n+1)*pi/2/a*ones(size(w));
kyn = sqrt((epr.*mur - 1).*k.^2 - kxn.^2);
sh  = sinh(kxn*b);
ch  = cosh(kxn*b);
tn  = tan(kyn*(b-d));
ct  = cot(kyn*(b-d));

Zl = conj(1i*Z0*L/2/a * sum(1./ ( ...
    ( kxn./k.*sh.*ch.*(1+epr.*mur) + kyn./k.*(mur.*sh.^2.*tn - epr.*ch.^2.*ct) ...
    )./(epr.*mur - 1) - k./kxn.*sh.*ch ...
    ),1));

Zv = conj(1i*Z0*L/2/a * sum(kxn.^2./k./ ( ...
    ( kxn./k.*sh.*ch.*(1+epr.*mur) + kyn./k.*(mur.*ch.^2.*tn - epr.*sh.^2.*ct) ...
    )./(epr.*mur - 1) - k./kxn.*sh.*ch ...
    ),1));

kxn = 2*(n+1)*pi/2/a*ones(size(w));
kyn = sqrt((epr.*mur - 1).*k.^2 - kxn.^2);
sh  = sinh(kxn*b);
ch  = cosh(kxn*b);
tn  = tan(kyn*(b-d));
ct  = cot(kyn*(b-d));

Zh = conj(1i*Z0*L/2/a * sum(kxn.^2./k./ ( ...
    ( kxn./k.*sh.*ch.*(1+epr.*mur) + kyn./k.*(mur.*sh.^2.*tn - epr.*ch.^2.*ct) ...
    )./(epr.*mur - 1) - k./kxn.*sh.*ch ...
    ),1));

