const = lnls_constants;
c = const.c;
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
ring.dtunedj  = 000000;

bunch.num_part = 600;
bunch.I_b      = 3.0e-3;

tau = (-1000:1000)*1e-12;
V = 3.0e6;
wrf = 2*pi*ring.har_num/ring.rev_time;
phi0 = 171.24/180*pi;
Vl = 1.0e6*0;
wl = wrf*3;
phil = 0.30/180*pi;
bunch.potential= V*(sin(phi0-wrf*tau)-sin(phi0)) + Vl*(sin(phil-wl*tau)-sin(phil));
bunch.tau      = tau;
bunch.espread  = 7.64e-4;
bunch.emit     = 271e-12;


Zovern = 0.2;
radius = 12e-3;
fr  = 2.4* c/(radius*2*pi); 
Rs = Zovern*fr*ring.rev_time;
Q = 1;
wr = fr*2*pi;
Ql = sqrt(Q.^2 - 1/4);
wrl = wr .* Ql ./ Q;

% load impedance;
% w = [0, logspace(4,18,20000)];
% budget = sirius_impedance_budget(w,'ring','phase_1');
% budgetbckup = budget;
% plot_impedances(w,budget,true,'log','log',false);
% budget = budgetbckup(end);



tau  = -[0, logspace(-14,-9,300)];
% [wakel,wakeh,wakev] = getwake_from_budget(w,budget,tau);


clear wake;
wake.long.sim            = true;
wake.long.track          = false;
wake.long.general.tau    = tau;
wake.long.general.wake = wakel;
% wake.long.resonator.wr   = wr;
% wake.long.resonator.Rs   = Rs;
% wake.long.resonator.Q    = Q;


beta_imp = 11;
Rs = Zovern*fr*ring.rev_time/radius;

wake.dipo.track          = true;
wake.dipo.general.tau  = tau;
wake.dipo.general.wake = wakev;
% wake.dipo.resonator.wr   = wr;
% wake.dipo.resonator.Rs   = Rs;
% wake.dipo.resonator.Q    = Q;
% wake.dipo.resonator.beta = beta_imp;


Rs = -Rs;
wake.quad.track          = false;
% wake.quad.general.tau  = tau;
% wake.quad.general.wake = beta_imp*wr*Rs/Ql*sin(wrl*tau).*exp(wr*tau/(2*Q));
% wake.quad.resonator.wr   = wr;
% wake.quad.resonator.Rs   = Rs;
% wake.quad.resonator.Q    = Q;
% wake.quad.resonator.beta = beta_imp;


wake.feedback.track       = false;
wake.feedback.npoints     = 8;
wake.feedback.freq        = 0.11;
wake.feedback.phase       = 3/4*pi;
wake.feedback.gain        = 0.1;
wake.feedback.delay       = 4;
wake.feedback.max_kick    = 3e-7;
wake.feedback.pickup_beta = 10;
wake.feedback.kick_beta   = 17;

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