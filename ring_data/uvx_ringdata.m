function ringdata = uvx_ringdata(phase)

ringdata.lattice_version = 'uvx_bby6t';

w = logspace(5,10,3000);
w2= logspace(10,11,8000);
w = [w w2(2:end)];
ringdata.w = [-fliplr(w) w];

%Tunes and chromaticities
ringdata.nuty = 5.19;
ringdata.nutx = 4.16;
ringdata.chromy = 0.0;
ringdata.chromx =- 0.0;

ringdata.nb = 148;           % number of bunches
ringdata.w0 = 2*pi*299792458/93.19991; % revolution angular frequency [Hz]

ringdata.eta = 8.3e-3;       % momentum compaction factor
ringdata.E = 1.37;              % energy [GeV];

if(strcmp(phase,'BBY6T'))
    ringdata.stage   = 'BBY6T';
    ringdata.sigma   = 10.2e-3;     % longitudinal length [m]
    ringdata.I_tot   = 0.25;       % total current [A];
    ringdata.nus     = 0.0078;    % synchrotron tune
    ringdata.ex      = 100e-9;    % horizontal emittance [m.rad]
    ringdata.ey      = 100e-11;    % vertical emittance [m.rad]
    ringdata.espread = 7e-4;     % energy spread   
    ringdata. Vrf    = 0.5;          % Tensao de RF  [MV]
    % damping times [s]
    ringdata.taux = 8.4e-3;
    ringdata.tauy = 7.5e-3;
    ringdata.taue = 3.5e-3;
end

