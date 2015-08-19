function  globdata = calc_kick_factor(globdata)
% function  globdata = calc_kick_factor(globdata)
%This function calculates the kick factor according to methodologies:
% a. using the long. wake data
% b. using the long impedance data

% Extracts and Initialize Needed Variables:
sigs  = globdata.simpar.sigma;
wake  = globdata.results.W;
saxis = globdata.results.s;
freq  = globdata.results.freq;
ImZ   = globdata.results.ImZt;

c = 299792458;

sigmasq = sigs^2;
w =(freq*2*pi);
k = w/c;
ksq = k.^2;


% Calculates kickZ vs. sigma:
sigmax = globdata.ringpar.sigmamax;
sigmin = globdata.simpar.sigma;
sigi = linspace(sigmin,sigmax,100);

rhok  = exp(-ksq.*sigs^2);
kickZ = trapz(k, ImZ.*rhok)*c/pi*1e-12;

for i=1:length(sigi)
    rhok = exp(-ksq.*sigi(i)^2);
    kickZi(i) = trapz(k, ImZ.*rhok)*c/pi*1e-12;
end

% Calculates kickW:
ss = saxis.^2;
rhos = (1/(sigs*sqrt(2*pi)))*exp(-ss./(2*sigmasq));
kickW = trapz(wake.*rhos, saxis);

% Assign results to structure:
globdata.results.kickZ = kickZi;
globdata.results.sigmak = sigi;
globdata.results.kickW = kickW;

% Print kick factor calculated in both ways:
disp(['Kick_Z, V/pC/m  ' num2str(kickZ)]);
disp(['Kick_W, V/pC/m  ' num2str(kickW)]);
