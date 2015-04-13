function ringdata = sirius_ringdata(phase)

ringdata.lattice_version = 'SI.V07.C05';

w = logspace(5,10,3000);
w2= logspace(10,12.4,8000); % limites determinados atraves de estudo do
w = [w w2(2:end)];
ringdata.w = [-fliplr(w) w];


%Tunes 
ringdata.nuty = 13.116;
ringdata.nutx = 48.131;

ringdata.nb = 864;           % number of bunches
ringdata.w0 = 2*pi*299792458/518.396; % revolution angular frequency [Hz]

ringdata.eta = 1.7e-4;       % momentum compaction factor
ringdata.E = 3;              % energy [GeV];

if(strcmp(phase,'commissioning'))
    ringdata.stage   = 'Commissioning';
    I = linspace(0,4,40)';
    fit = [3.22e-3, 1.75e-3, -1.86e-3]; % fitting dos dados da natalia;
    ringdata.sigma   = [I, fit(1) + I*fit(2) + I.^2*fit(3)];  % longitudinal length [m]
    ringdata.I_tot   = 0.100;       % total current [A];
    ringdata.nus     = 0.00435;    % synchrotron tune
    % damping times [s]
    ringdata.taux =17.0e-3;
    ringdata.tauy = 22.6e-3;
    ringdata.taue =12.4e-3;
elseif(strcmp(phase,'phase_1'))
    ringdata.stage   = 'Phase 1 with IDS';
    I = linspace(0,4,40)';
    fit = [3.83e-3, 8.96e-4, -5.24e-4]; % fitting dos dados da natalia;
    ringdata.sigma   = [I, fit(1) + I*fit(2) + I.^2*fit(3)];  % longitudinal length [m]
    ringdata.I_tot = 0.10;         % total current [A];
    ringdata.nus = 0.00435;        % synchrotron tune
    % damping times [s]
    ringdata.taux = 12.9e-3;
    ringdata.tauy = 14.95e-3;
    ringdata.taue = 8.42e-3;
elseif(strcmp(phase,'phase_2'))
    ringdata.stage   = 'Phase 2 with IDS';
    I = linspace(0,4,40)';
    fit = [14.0e-3, 3.38e-4, -3.0e-5]; % fitting dos dados da natalia;
    ringdata.sigma   = [I, fit(1) + I*fit(2) + I.^2*fit(3)];  % longitudinal length [m]
       
    ringdata.I_tot   = 0.35;        % total current [A];
    ringdata.nus     = 0.00435;    % synchrotron tune
    % damping times [s]
    ringdata.taux = 10.37e-3;
    ringdata.tauy = 12.3e-3;
    ringdata.taue = 6.77e-3;
elseif(strcmp(phase,'phase_2_HC'))
    ringdata.stage   = 'Phase 2 with IDS High Current';
    I = linspace(0,4,40)';
    fit = [14.0e-3, 3.38e-4, -3.0e-5]; % fitting dos dados da natalia;
    ringdata.sigma   = [I, fit(1) + I*fit(2) + I.^2*fit(3)];  % longitudinal length [m]
    
    ringdata.I_tot   = 0.5;        % total current [A];
    ringdata.nus     = 0.00435;    % synchrotron tune
    % damping times [s]
    ringdata.taux = 10.37e-3;
    ringdata.tauy = 12.3e-3;
    ringdata.taue = 6.77e-3;
elseif(strcmp(phase,'single_bunch'))
    ringdata.stage   = 'Single Bunch';
    ringdata.sigma   = 3.83e-3;     % longitudinal length [m]
    ringdata.I_tot   = 0.002;      % total current [A];
    ringdata.nus     = 0.00435;    % synchrotron tune
    % damping times [s]
    ringdata.taux = 10.37e-3;
    ringdata.tauy = 12.3e-3;
    ringdata.taue = 6.77e-3;
    
    ringdata.nb = 1;
end



