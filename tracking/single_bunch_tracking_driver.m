const = lnls_constants;
c = const.c;
el_ch = const.q0;
% ring data
ring.nturns   = 10000;
ring.rev_time = 518.396/c;
ring.E        = 3e9;
ring.mom_comp = 1.7e-4;
ring.beta     = 11;
ring.eta      = 0.0;
ring.etap     = 0;
ring.har_num  = 864;
ring.tune     = 13.117;
ring.dtunedp  = 0.0;
ring.dtunedj  = 200000;

bunch.num_part = 10000;
bunch.I_b      = 3.0e-3;

tau = (-1000:1000)*1e-12;
V = 3.0e6;
wrf = 2*pi*ring.har_num/ring.rev_time;
phi0 = 171.24/180*pi;
Vl = 1.0e6*1;
wl = wrf*3;
phil = 0.30/180*pi;
bunch.potential= V*(sin(phi0-wrf*tau)-sin(phi0)) + Vl*(sin(phil-wl*tau)-sin(phil));
bunch.tau      = tau;
bunch.espread  = 7.64e-4;
bunch.emit     = 271e-12;

% Resistive Wall Impedance
radius = 12e-3;
sigma = 5.9e7;
Z0 = c*4*pi*1e-7;
a = 3/sqrt(2*pi)/4;
p = 2.7;
L = 480;


Zovern = 0.2;
fr  = 2.4* c/(radius*2*pi); 
Rs = Zovern*fr*ring.rev_time;
Q = 1;
wr = fr*2*pi;
Ql = sqrt(Q.^2 - 1/4);
wrl = wr .* Ql ./ Q;

tau = -(0:1000)*1e-12;
clear wake;
wake.long.sim            = true;
wake.long.track          = false;
wake.long.general.tau    = tau;
% wake.long.general.wake = wr*Rs/Q*(cos(wrl*tau) + 1/(2*Ql)*sin(wrl*tau)).*exp(wr*tau/(2*Q));
% wake.long.general.wake(1) = wake.long.wake(1)/2;
wake.long.resonator.wr   = wr;
wake.long.resonator.Rs   = Rs;
wake.long.resonator.Q    = Q;
wake.long.wall.sigma     = 0;
wake.long.wall.fullgap   = radius;
wake.long.wall.length    = L;

beta_imp = 11;
Rs = Zovern*fr*ring.rev_time/radius;

tau = -(0:1000)*1e-12;
s0 = (2*radius^2/Z0/sigma)^(1/3);
W0 = c*Z0/4/pi * 2*s0*L/radius^4;
wake = beta_imp*W0*(8/3*exp(tau*c/s0).*sin(-tau*c*sqrt(3)/s0-pi/6) ...
    + 1/sqrt(2*pi)*1./(a^p + (tau*c/(4*s0)).^(p/2)).^(1/p));


wake.dipo.track          = true;
% wake.dipo.general.tau  = tau;
% wake.dipo.general.wake = beta_imp*wr*Rs/Ql*sin(wrl*tau).*exp(wr*tau/(2*Q));
wake.dipo.resonator.wr   = wr;
wake.dipo.resonator.Rs   = Rs;
wake.dipo.resonator.Q    = Q;
wake.dipo.resonator.beta = beta_imp;
wake.dipo.wall.sigma     = 0;
wake.dipo.wall.fullgap   = radius;
wake.dipo.wall.length    = L;

Rs = -Rs;
wake.quad.track          = false;
% wake.quad.general.tau  = tau;
% wake.quad.general.wake = beta_imp*wr*Rs/Ql*sin(wrl*tau).*exp(wr*tau/(2*Q));
wake.quad.resonator.wr   = wr;
wake.quad.resonator.Rs   = Rs;
wake.quad.resonator.Q    = Q;
wake.quad.resonator.beta = beta_imp;
wake.quad.wall.sigma     = 0;
wake.quad.wall.fullgap   = radius;
wake.quad.wall.length    = L;

wake.feedback.track = false;
wake.feedback.npoints = 8;
wake.feedback.freq   = 0.11;
wake.feedback.phase  = 3/4*pi;
wake.feedback.gain   = 0.1;

[ave_bun,rms_bun, ave_kickx, fdbkx] = single_bunch_tracking(ring, bunch, wake);

% I_b = linspace(0.05,2.0,30)*1e-3;
% fft_ave = zeros(30,ring.nturns);
% rmsx    = zeros(30,ring.nturns);
% for i=1:30
%     bunch.I_b = I_b(i);
%     [ave_bun,rms_bun, ave_kickx, fdbkx] = single_bunch_tracking(ring, bunch, wake);
%     fft_ave(i,:) = 2*abs(fft(ave_bun(1,:)));
%     rmsx(i,:) = rms_bun(1,:);
%     fprintf('%d : %5.3f mA\n',i,I_b(i)*1e3);
% end
% n = ring.nturns;
% tune = (0:n/2)/n;
% ind = tune > 0.11 & tune < 0.125;
% [I,T] = meshgrid(I_b,tune);
% 
% pfft = fft_ave(:,1:n/2+1)';
% pfft = pfft./repmat(max(pfft),n/2+1,1);
% figure;  surface(I(ind,:),T(ind,:),pfft(ind,:),'LineStyle','none');
% xlim([min(I_b),max(I_b)]);ylim([min(tune(ind)),max(tune(ind))]);
% 
% [I,N] = meshgrid(1:n,I_b);
% figure; surface(N,I,log(rmsx/min(min(rmsx))),'LineStyle','none');
% xlim([1,n]);ylim([min(I_b),max(I_b)]);