function ringdata = booster_ringdata(phase)

ringdata.lattice_version = 'Booster_Sirius_v700-opt_2';

w = logspace(5,9,3000);
w2= logspace(9,11.48,800000);
w = [w w2(2:end)];
ringdata.w = [-fliplr(w) w];


%Tunes and chromaticities
ringdata.nuty = 7.3032;
ringdata.nutx = 19.1904;
ringdata.chromy = 0.0;
ringdata.chromx =- 0.0;

ringdata.nb = 828;           % number of bunches
ringdata.w0 = 2*pi*299792458/496.8; % revolution angular frequency [Hz]

ringdata.eta = 7.01e-4;       % momentum compaction factor
ringdata.I_tot   = 300e-3; %6.1e-4;       % total current [A];

if(strcmp(phase,'high_energy'))
    ringdata.stage   = 'Extraction Energy - 3GeV';
    ringdata.E       = 3;             % energy [GeV];
    ringdata.sigma   = 11.3e-3;     % longitudinal length [m]
    ringdata.nus     = 0.0044;    % synchrotron tune
    ringdata.espread = 8.95e-4;     % energy spread   
    ringdata. Vrf    = 1;          % Tensao de RF  [MV]
    % damping times [s]
    ringdata.taux = 10.6e-3;
    ringdata.tauy = 12.7e-3;
    ringdata.taue =  7e-3;
elseif(strcmp(phase,'low_energy'))
    ringdata.stage   = 'Injection Energy - 150MeV';
    ringdata.E       = 0.15;             % energy [GeV];
    ringdata.espread = 2.5e-3;     % energy spread  
    ringdata.nus     = 0.0096;    % synchrotron tune
    % o valor nao eh baseado no tamanho de equilibrio, mas sim na rampa
    ringdata.sigma   = 9e-3;     % longitudinal length [m]
    ringdata. Vrf    = 0.15;          % Tensao de RF  [MV]
    % damping times [s]
    ringdata.taux    =  125e-3; % 83;
    ringdata.tauy    =  125e-3; %101; %coloquei aqui o tempo de meia rampa
    ringdata.taue    =  125e-3; %57;
elseif(strcmp(phase,'medium_energy'))
    ringdata.stage   = 'Medium Energy - 1.5GeV';
    ringdata.E       = 1.5;             % energy [GeV];
    ringdata.espread = 1e-3;     % energy spread
    ringdata.nus     = 0.0096;    % synchrotron tune
    ringdata.sigma   = 4.5e-3;     % longitudinal length [m]
    ringdata. Vrf    = 0.65;          % Tensao de RF  [MV]
    % damping times [s]
    ringdata.taux    =  83e-3; % 83;
    ringdata.tauy    =  101e-3; %101; %coloquei aqui o tempo de meia rampa
    ringdata.taue    =  57e-3; %57;   
end