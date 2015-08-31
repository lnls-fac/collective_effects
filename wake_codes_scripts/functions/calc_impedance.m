function  globdata = calc_impedance(globdata)
%function  globdata = calc_impedance(globdata)
%This script calculates impedance from Wakepotential data 


%% Extracts Needed Variables
m     = globdata.simpar.m;
sigs  = globdata.simpar.sigma;
wake  = globdata.results.W * 1e12; % rescale W to [V/C]
saxis = globdata.results.s;
f0    = globdata.ringpar.omega0/(2*pi);

c     = 299792458;
sigt  = sigs/c;

%% Impedance calculation based on Wakepotential

% translates scale to zero origin
wakeorigin = saxis-min(saxis);

% frequency scale (Hz)
Np = length(saxis);
ds = max(wakeorigin)/(Np-1);
dt = ds/c;

% Scale shift/zero shift

shift = (saxis(1))/c;

% calculates FFT and frequency scale
fftt = fft(-wake',Np)';
p    = 1:Np;
w    = 2*pi*(p-1)/(Np*dt);
VHat = (dt*exp(1i*w*shift).*fftt);
freq = w/(2*pi);

% Spectrum of a gaussian bunch:
Jwlist = exp(-(w*sigt).^2/2);
Z      = VHat./Jwlist;

%% Exports calculated FFT - there is a 10 MHz difference between the peaks, which corresponds to 1/2 df

%defines the frequency range due to the bunch length
wmax = globdata.simpar.cutoff/sigt;
indcs = find(w <= wmax);

if m==0
    globdata.results.freq = freq(indcs);
    globdata.results.ReZlong =  real(Z(indcs));
    globdata.results.ImZlong = -imag(Z(indcs));
    globdata.results.naxis = globdata.results.freq(2:4)./f0;
    globdata.results.ImZoN = globdata.results.ImZlong(2:4)./globdata.results.naxis;
elseif m>0
    globdata.results.freq =   freq(indcs);
    globdata.results.ReZt = imag(Z(indcs));
    globdata.results.ImZt = real(Z(indcs));
end