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

bunch.num_part = 1;
bunch.fillPat = 1:864;
bunch.I_b     = ones(1,length(bunch.fillPat),1)*500e-3/length(bunch.fillPat);

tau = (-1000:1000)*1e-12;
V = 3.0e6;
wrf = 2*pi*ring.har_num/ring.rev_time;
phi0 = 171.24/180*pi;
Vl = 1.0e6*0;
wl = wrf*3;
phil = 0.30/180*pi;
bunch.potential= V*(sin(phi0-wrf*tau)-sin(phi0)) + Vl*(sin(phil-wl*tau)-sin(phil));
bunch.potential = repmat(bunch.potential,length(bunch.fillPat),1);
bunch.tau      = tau;
bunch.espread  = 7.64e-4;
bunch.emit     = 271e-12;


Zovern = 0.2;
radius = 12e-3;
fr  = 3.1e9; 
Rs = Zovern*fr*ring.rev_time;
Q = 1000;
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
wake.long.track          = true;
% wake.long.general.tau    = tau;
% wake.long.general.wake = wakel;
wake.long.resonator.wr   = wr;
wake.long.resonator.Rs   = Rs;
wake.long.resonator.Q    = Q;
wake.long.resonator.Memory = 0;

beta_imp = 11;
Rs = Zovern*fr*ring.rev_time/radius;

wake.dipo.track          = false;
% wake.dipo.general.tau  = tau;
% wake.dipo.general.wake = wakev;
wake.dipo.resonator.wr   = wr;
wake.dipo.resonator.Rs   = Rs;
wake.dipo.resonator.Q    = Q;
wake.dipo.resonator.beta = beta_imp;
wake.dipo.resonator.Memory = 0;

Rs = -Rs;
wake.quad.track          = false;
% wake.quad.general.tau  = tau;
% wake.quad.general.wake = beta_imp*wr*Rs/Ql*sin(wrl*tau).*exp(wr*tau/(2*Q));
wake.quad.resonator.wr   = wr;
wake.quad.resonator.Rs   = Rs;
wake.quad.resonator.Q    = Q;
wake.quad.resonator.beta = beta_imp;
wake.quad.resonator.Memory = 0;

wake.feedback.track = false;
wake.feedback.npoints = 8;
wake.feedback.freq   = 0.11;
wake.feedback.phase  = 3/4*pi;
wake.feedback.gain   = 0.1;

[x, tau] = multi_bunch_tracking(ring, bunch, wake);
