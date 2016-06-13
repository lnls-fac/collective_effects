function imp = multilayer_impedance(w)

mu0 = 4*pi*1e-7;
c   = 299792458;
ep0 = 1/c^2/mu0;
E = 3;

% normal wall
% epb     = [1 1 1];
% mub     = [1 1 1];
% ange    = [0 0 0];
% angm    = [0 0 0];
% sigmadc = [0 5.9e7 1]; % Not sure about NEG resistivity
% relax_time     = [0 1 0]*27e-15;
% b       = [12.000 13.000]*1e-3;
% L       = 480;

% wall with coating:
% epb     = [1 1 1 1];
% mub     = [1 1 1 1];
% ange    = [0 0 0 0];
% angm    = [0 0 0 0];
% sigmadc = [0 4e6 5.9e7 1]; % Not sure about NEG resistivity
% relax_time  = [0 0 1 0]*27e-15;
% b       = [12.000 12.001 13.000]*1e-3;
% L       = 480;

% in vacuum Undulators
% epb     = [1 1 1 1];
% mub     = [1 1 1 100];
% ange    = [0 0 0 0];
% angm    = [0 0 0 0];
% sigmadc = [0 5.9e7 1 6.25e5]; % Copper Sheet
% relax_time = [0 1 0 0]*27e-15;
% b       = [4.5 4.65 4.7]/2*1e-3;
% L       = 2.0;

% EPUs
% epb     = [1 1 1];
% mub     = [1 1 1];
% ange    = [0 0 0];
% angm    = [0 0 0];
% sigmadc = [0 5.9e7 1]; % Copper Sheet
% relax_time = [0 1 0 ]*27e-15;
% b       = [12 14]/2*1e-3;
% L       = 2.7;

% Fast Correctors
% epb     = [1 1 1 1];
% mub     = [1 1 1 1];
% ange    = [0 0 0 0];
% angm    = [0 0 0 0];
% sigmadc = [0 4e6 1.3e6 0.1]; % Not sure about NEG resistivity
% relax_time = [0 0 0 0]*27e-15;
% b       = [12.000 12.001 12.3]*1e-3;
% L       = 0.1 * 80;

% PMM
epb     = [1   1    1    9.3   1];
mub     = [1   1    1     1    1];
ange    = [0   0    0     0    0];
angm    = [0   0    0     0    0];
sigmadc = [0  4e6  5.9e7  1    1]; % Copper Sheet
relax_time = [0   0    1     0    0]*27e-15;

coat    = 2e-3;
neg     = 1e-3;
b       = [(4.5-coat-neg) (4.5-coat) 4.5 5.5]*1e-3;
L       = 0.5;

% n=10;
% fil   = exp(-((-n:n)/(n/5)).^2)/sqrt(pi)/n*5;

epr = zeros(length(epb),length(w));
mur = zeros(length(epb),length(w));
for j = 1: length(epb)
    epr(j,:) = epb(j)*(1-1i.*sign(w).*tan(ange(j))) + sigmadc(j)./(1+1i*w*relax_time(j))./(1i*w*ep0);
    mur(j,:) = mub(j)*(1-1i.*sign(w).*tan(angm(j)));
end
[Zl, Zv, ~] = lnls_calc_impedance_multilayer_round_pipe(w, epr, mur, b, L, E);

imp = Zv;