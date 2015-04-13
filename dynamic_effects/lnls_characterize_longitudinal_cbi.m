function deltaw = lnls_characterize_longitudinal_cbi(ringdata, budget, azi_modes, save)

w = ringdata.w;

%% Load parameters

tau  = ringdata.taue;
nb = ringdata.nb;
w0 = ringdata.w0; 
nus = ringdata.nus;
eta = ringdata.eta;
I_tot = ringdata.I_tot;
sigma = interp1(1e-3*ringdata.sigma(:,1)',ringdata.sigma(:,2)',I_tot/nb);
E = ringdata.E;

Zl = zeros(1,length(w));
imped  = 1:length(budget);

for j=imped
    Z     = budget{j}.Zl;
    quant = budget{j}.quantity;
    Zl    = Zl + Z*quant;
end

%% Calc coherent tune-shifts
fprintf('Calculation of Longitudinal Coupled Bunch Instability and Tune Shifts\n');
fprintf('%-20s: %-20.4g\n','Shynchrotron Tune', nus);
fprintf('%-20s: %-20.4g\n\n','Damping Time [ms]', tau*1e3);
fprintf('Number of unstable Modes Modes\n');
fprintf('%-5d ', azi_modes);fprintf('\n');
fprintf('%s',repmat('-',1,30)); fprintf('\n');

deltaw = [];

for m=azi_modes
    deltaw0 = lnls_calc_longitudinal_cbi(w, Zl, sigma, nb, w0, nus, eta, E, I_tot, m);
    deltaw = [ deltaw; deltaw0];
    n = imag(deltaw0) > 1/tau ;
    n = sum(n);
    fprintf('%-5d ',n);
end

fprintf('\n\n\n');

rel_tuneshift  = real(deltaw)/w0/nus;
rel_growth   = imag(deltaw)*tau;

%% Plot Results

f2 = figure('Position',[1 1 850 528]);

axes3 = axes('Parent',f2,'Position',[0.12 0.103 0.850 0.415], 'FontSize',16);
box(axes3,'on'); hold(axes3,'all'); grid(axes3,'on');

plot((1:nb)-1,rel_tuneshift,'Parent',axes3,'LineWidth',2);
xlabel({'Coupled Bunch Mode'},'FontSize',16);
ylabel({'Re(\Omega - \nu_{\beta})/\nu_s'},'FontSize',16);
xlim([0 (nb-1)]);


axes4 = axes('Parent',f2,'Position',[0.12 0.518 0.850 0.415],...
             'FontSize',16, 'XTickLabel',{''});
box(axes4,'on');hold(axes4,'all');grid(axes4,'on');
plot1 = plot(axes4, (1:nb)-1,rel_growth,'LineWidth',2,'DisplayName','Chrom=0, m=0');
plot(axes4, (1:nb)-1,repmat(1/tau,1,nb),'LineWidth',2,...
     'LineStyle','--','Color','k','DisplayName','1/\tau_d');
for i=1:length(azi_modes)
    leg = sprintf('m = %d',azi_modes(i));
    set(plot1(i),'DisplayName',leg);
end
legend(axes4,'show','Location','Best');
ylabel({'1/\tau_{damp} [Hz]'},'FontSize',16);
xlim([0 (nb-1)]);


if save
    saveas(f2,['cbi_l_mode_' ringdata.stage '.fig']);
end
