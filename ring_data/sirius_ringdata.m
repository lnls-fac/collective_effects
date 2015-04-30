function ringdata = sirius_ringdata(phase)

ringdata.lattice_version = 'SI.V07.C05';

w = logspace(5,12.4,11000); % limites determinados atraves de estudo do
% w = linspace(5e8,3e9,600000)*2*pi;
ringdata.w = [-fliplr(w) w];


%Tunes 
ringdata.nuty = 13.116;
ringdata.nutx = 48.131;
ringdata.nb = 864;           % harmonic Number
ringdata.w0 = 2*pi*299792458/518.396; % revolution angular frequency [Hz]
ringdata.eta = 1.7e-4;       % momentum compaction factor
ringdata.E = 3;              % energy [GeV];

I = linspace(0,4,40);
ringdata.I   = I*1e-3;
if(strcmp(phase,'commissioning'))
    ringdata.stage   = 'Commissioning';
    ringdata.I_tot   = 0.100;       % total current [A];
    ringdata.nus     = 0.00435;    % synchrotron tune
    ringdata.espread = 7.64e-4 + 0*I;
    ringdata.emitx   = 271e-12 + 0*I;
    ringdata.emity   = 2.71e-12 + 0*I;
    % damping times [s]
    ringdata.taux = 17.1e-3;
    ringdata.tauy = 22.7e-3;
    ringdata.taue = 13.6e-3;
    ringdata.U0   = 456740.6; %eV
elseif(strcmp(phase,'phase_1'))
    ringdata.stage   = 'Phase 1 with IDS';
    ringdata.I_tot = 0.10;         % total current [A];
    ringdata.nus = 0.00435;        % synchrotron tune
    ringdata.espread = 1e-02*(0.093995 + 0.038011*I -0.018279*I.^2 + 0.0047843*I.^3 -0.00047294*I.^4);
    ringdata.emitx   = 1e-09*(0.23011  + 0.15699 *I -0.063581*I.^2 + 0.015965* I.^3 -0.0015505* I.^4);
    ringdata.emity   = 1e-12*(2.1496   + 1.8725*  I -0.84932 *I.^2 + 0.22507*  I.^3 -0.022538 * I.^4);
    % damping times [s]
    ringdata.taux = 12.4e-3;
    ringdata.tauy = 15.1e-3;
    ringdata.taue =  8.5e-3;
    ringdata.U0   = 685374.1; %eV
elseif(strcmp(phase,'phase_2'))
    ringdata.stage   = 'Phase 2 with IDS';
    ringdata.I_tot   = 0.35;        % total current [A];
    ringdata.nus     = 0.00435;    % synchrotron tune
    ringdata.espread = 1e-02*(0.088704 + 0.015765*I -0.005477*I.^2 + 0.0012452*I.^3 -0.00011434*I.^4);
    ringdata.emitx   = 1e-09*(0.18859  + 0.056781*I -0.015909*I.^2 + 0.003445* I.^3 -0.00031039*I.^4);
    ringdata.emity   = 1e-12*(1.6497   + 1.0422*  I -0.51454 *I.^2 + 0.14498*  I.^3 -0.015059 * I.^4);
    % damping times [s]
    ringdata.taux = 10.6e-3;
    ringdata.tauy = 12.5e-3;
    ringdata.taue =  6.9e-3;
    ringdata.U0   = 829761.9; %eV
elseif(strcmp(phase,'phase_2_HC'))
    ringdata.stage   = 'Phase 2 with IDS High Current';
    ringdata.I_tot   = 0.5;        % total current [A];
    ringdata.nus     = 0.00435;    % synchrotron tune
    ringdata.espread = 1e-02*(0.088704 + 0.015765*I -0.005477*I.^2 + 0.0012452*I.^3 -0.00011434*I.^4);
    ringdata.emitx   = 1e-09*(0.18859  + 0.056781*I -0.015909*I.^2 + 0.003445* I.^3 -0.00031039*I.^4);
    ringdata.emity   = 1e-12*(1.6497   + 1.0422*  I -0.51454 *I.^2 + 0.14498*  I.^3 -0.015059 * I.^4);
    % damping times [s]
    ringdata.taux = 10.6e-3;
    ringdata.tauy = 12.5e-3;
    ringdata.taue =  6.9e-3;
    ringdata.U0   = 829761.9; %eV
end



