function  globdata = calc_loss_factor(globdata)
% function  alt = calc_loss_factor(globdata)
%This script calculates the loss factor according to methodologies:
% a. using the long. wake data
% b. using the long impedance data

% Extracts and Initialize Needed Variables:
h    = globdata.ringpar.h;
T0   = 2*pi/globdata.ringpar.omega0;
Iavg = globdata.ringpar.Iavg;
sigs = globdata.simpar.sigma;

wake  = globdata.results.W;
saxis = globdata.results.s;
freq  = globdata.results.freq;
ReZ   = globdata.results.ReZlong;

c = 299792458;

omega   = (freq*2*pi);
k       = omega./c;
ksq     = k.^2;

% Calculates klossZ vs. sigma:
sigmax = globdata.ringpar.sigmamax;
sigmin = globdata.simpar.sigma;
sigi   = linspace(sigmin,sigmax,100);

for i=1:length(sigi)
    rhok   = exp(-ksq.*sigi(i)^2);
    kZi(i) = trapz(k, ReZ.*rhok)*c/pi*1e-12;
end

kZ = kZi(1);

sigvec = [2.65 5.3 2.65 4 10 10]*1e-3;  % bunch length scenarios
Ivec   = [500 500 10 110 110 500]*1e-3; % current scenarios

for i=1:length(sigvec)
    rhok     = exp(-ksq.*sigvec(i)^2);
    kZvec(i) = trapz(k, ReZ.*rhok)*c/pi*1e-12;
end
Plossvec = kZvec.*Ivec.^2*T0*1e12/h;

globdata.results.klossZ   = kZi;
globdata.results.sigmak   = sigi;
globdata.results.Plossvec = Plossvec;

%% Calculates klossW

ss    = saxis.^2;
rhos  = (1/(sigs*sqrt(2*pi)))*exp(-ss./(2*sigs^2));
kW    = trapz(wake.*rhos, saxis);
Ploss = kW*Iavg^2*T0*1e12/h;

globdata.results.klossW = kW;
globdata.results.Ploss  = Ploss;

%% Print loss factor calculated in both ways

disp(['klossZ = ' num2str(kZ*1000) ' mV/pC']);
disp(['klossW = ' num2str(kW*1000) ' mV/pC']);
disp(['Ploss = ' num2str(Ploss) ' W (for ' num2str(Iavg*1000) ' mA avg current)']);


