c = 299792458;


w = linspace(10^1,10^12,3000000);
% w = [-fliplr(w) w];


Rsl = 5.3e3; % = 3.6*518.25/354.0*1e3;
wrl = 20*1e9*2*pi;
Ql =  1000;
Zl = lnls_calc_impedance_longitudinal_resonator(Rsl, Ql, wrl, w);

sigz= 3.3e-3 ;
h = exp(-(sigz*w/c).^2);

tau = (-20:0.005:20)*1e-3;
q   = tau / sigz;

lamb = zeros(size(q));
int = zeros(size(q));

for jj =1:length(q)
    int(jj) = trapz(Zl.*h.*exp(1i*w*q(jj)*sigz/c));
end

I_b = 1e-4;
re = 2.8179e-15;
T0 = 518.25/c;
Nb = I_b * T0 / 1.6e-19;
nus = 0.0039;
gamma = 5800;
sigd  = 7e-4;

In = Nb* re/(2*pi*nus*gamma*sigd);
for jj =1:length(tau)
    lamb(jj) = exp(-q(jj)^2 + In*trapz(int(1:jj)));
end
 
figure; plot(tau,real(lamb));