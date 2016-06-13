function ringdata = sirius_ringdata(phase)
% Sirius ring data;
c = 299792458;

% damping times [s]
ringdata.taux =15.9e-3;
ringdata.tauy = 21e-3;
ringdata.taue =12.4e-3;

%beta functions, tunes and chromaticities
ringdata.betax = 2; %[m]
ringdata.betay = 10;%[m]
ringdata.nuty = 4.17;
ringdata.nutx = 5.27;
ringdata.chromy = 0;
ringdata.chromx = 0;

ringdata.sigma = 11e-3;     % longitudinal length [m]
ringdata.nb = 148;           % number of bunches
ringdata.w0 = 2*pi*c/93.1999; % revolution angular frequency [Hz]
ringdata.nus = 0.0078;      % synchrotron tune

ringdata.eta = 8.29e-3;       % momentum compaction factor
ringdata.I_tot = 0.250;      % total current [A];
ringdata.E = 1.37;              % energy [GeV];
ringdata.lattice_version = 'Sirius_v403-ac10_5';

w = logspace(5,10,3000);
w2= logspace(10,12,8000);
w = [w w2(2:end)];
ringdata.w = [-fliplr(w) w];


%Tunes and chromaticities
ringdata.nuty = 14.147;
ringdata.nutx = 46.171;
ringdata.chromy = 0.0;
ringdata.chromx =- 0.0;

ringdata.nb = 864;           % number of bunches
ringdata.w0 = 2*pi*299792458/518.25; % revolution angular frequency [Hz]

ringdata.eta = 1.7e-4;       % momentum compaction factor
ringdata.E = 3;              % energy [GeV];

if(strcmp(phase,'commissioning'))
    ringdata.stage   = 'Commissioning';
    ringdata.sigma   = 3.3e-3;     % longitudinal length [m]
    ringdata.I_tot   = 0.01;       % total current [A];
    ringdata.nus     = 0.00393;    % synchrotron tune
    ringdata.ex      = 2.81e-10;    % horizontal emittance [m.rad]
    ringdata.ey      = 2.81e-12;    % vertical emittance [m.rad]
    ringdata.espread = 8.3e-4;     % energy spread   
    ringdata. Vrf    = 2;          % Tensao de RF  [MV]
    % damping times [s]
    ringdata.taux =15.9e-3;
    ringdata.tauy = 21e-3;
    ringdata.taue =12.4e-3;
elseif(strcmp(phase,'ids'))
    ringdata.stage   = 'IDS';
    ringdata.sigma = 3.39e-3;      % longitudinal length [m]
    ringdata.I_tot = 0.10;         % total current [A];
    ringdata.nus = 0.00393;        % synchrotron tune
    ringdata.ex      = 2.88e-10;   % horizontal emittance [m.rad]
    ringdata.ey      = 2.62e-12;   % vertical emittance [m.rad]
    ringdata.espread = 10.7e-4;    % energy spread   
    ringdata. Vrf    = 2.7;        % Tensao de RF  [MV]
    % damping times [s]
    ringdata.taux = 9.4e-3;
    ringdata.tauy = 12.4e-3;
    ringdata.taue = 7.4e-3;
elseif(strcmp(phase,'ids_landau'))
    ringdata.stage   = 'IDS + Landau';
    ringdata.sigma   = 13.4e-3;     % longitudinal length [m]
    ringdata.I_tot   = 0.10;        % total current [A];
    ringdata.nus     = 0.00393;     % synchrotron tune
    ringdata.ex      = 2.66e-10;    % horizontal emittance [m.rad]
    ringdata.ey      = 2.61e-12;    % vertical emittance [m.rad]
    ringdata.espread = 10.1e-4;     % energy spread   
    ringdata. Vrf    = 2.7;         % Tens??o de RF  [MV]
    % damping times [s]
    ringdata.taux = 9.4e-3;
    ringdata.tauy = 12.4e-3;
    ringdata.taue = 7.4e-3;
elseif(strcmp(phase,'full_current'))
    ringdata.stage   = 'Full Current';
    ringdata.sigma   = 13.4e-3;    % longitudinal length [m]
    ringdata.I_tot   = 0.3;        % total current [A];
    ringdata.nus     = 0.00393;    % synchrotron tune
    ringdata.ex      = 2.67e-10;   % horizontal emittance [m.rad]
    ringdata.ey      = 2.61e-12;   % vertical emittance [m.rad]
    ringdata.espread = 10.1e-4;    % energy spread   
    ringdata. Vrf    = 2.7;        % Tensao de RF  [MV]
    % damping times [s]
    ringdata.taux = 9.4e-3;
    ringdata.tauy = 12.4e-3;
    ringdata.taue = 7.4e-3;
elseif(strcmp(phase,'single_bunch'))
    ringdata.stage   = 'Single Bunch';
    ringdata.sigma   = 3.83e-3;     % longitudinal length [m]
    ringdata.I_tot   = 0.002;      % total current [A];
    ringdata.nus     = 0.00393;    % synchrotron tune
    ringdata.ex      = 2.67e-10;   % horizontal emittance [m.rad]
    ringdata.ey      = 2.61e-12;   % vertical emittance [m.rad]
    ringdata.espread = 10.1e-4;    % energy spread   
    ringdata. Vrf    = 2.7;        % Tensao de RF  [MV]
    % damping times [s]
    ringdata.taux = 9.4e-3;
    ringdata.tauy = 12.4e-3;
    ringdata.taue = 7.4e-3;
    
    ringdata.nb = 1;
end



