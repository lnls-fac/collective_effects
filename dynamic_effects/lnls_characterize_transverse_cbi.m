% function lnls_characterize_transverse_cbi(ringdata, budget,plane, n_azi, n_rad, chrom, save)
function deltaw = lnls_characterize_transverse_cbi(ringdata, budget,plane, azi_modes, chrom, save)
% n_rad = 4;
% n_azi = 6;
w = ringdata.w;

%% Load parameters
if strcmp(plane,'v')
    betas = 'betay';
    Zts    = 'Zv';
    tau  = ringdata.tauy;
    nut = ringdata.nuty;
    label = 'Vertical';
else
    betas = 'betax';
    Zts    = 'Zh';
    tau  = ringdata.taux;
    nut = ringdata.nutx;
    label = 'Horizontal';
end

nb = ringdata.nb;
w0 = ringdata.w0; 
nus = ringdata.nus;
eta = ringdata.eta;
I_tot = ringdata.I_tot;
sigma = interp1(1e-3*ringdata.sigma(:,1)',ringdata.sigma(:,2)',I_tot/nb);
E = ringdata.E;

Zt = zeros(1,length(w));
imped  = 1:length(budget);

for j=imped
    Z     = budget{j}.(Zts);
    beta  = budget{j}.(betas);
    quant = budget{j}.quantity;
    Zt    = Zt + Z*beta*quant;
end


%% new method

% fprintf('Calculation of %s Coupled Bunch Instability and Tune Shifts\n', label);
% fprintf('%-20s: %-20.4g\n','Betatron Tune', nut);
% fprintf('%-20s: %-20.4g\n\n','Damping Time [ms]', tau*1e3);
% fprintf('%-7s # of unstable Modes Modes\n',' ');
% fprintf('%-7s ','Chrom');fprintf('%-7.2f ', chrom);fprintf('\n');
% fprintf('%s',repmat('-',1,8*length(chrom))); fprintf('\n');
% fprintf('%-7s ','#');
% 
% 
% n_uns_modes = zeros(1,length(chrom));
% max_rate = zeros(1,length(chrom));
% max_shift = zeros(1,length(chrom));
% g_rate = zeros(1,nb);
% tu_shift = zeros(1,nb);
% for ii=1:length(chrom)
%     for mu=0:(nb-1)
%         deltaw =lnls_calc_transverse_mode_couplingopt(w, Zt, n_rad, n_azi, sigma, ...
%             I_tot/nb, E, w0, nut, nus, eta, chrom(ii), nb, mu);
%         [gmax, ~] = max(imag(deltaw));
%         g_rate(mu + 1) = gmax*nus*w0;
%         tu_shift(mu + 1) =  max(abs(real(deltaw)) - round(abs(real(deltaw))));
%     end
%     if chrom(ii)==0;
%         deltaw0 = 1i*g_rate + tu_shift;
%     end
%     n = g_rate > 1/tau ;
%     n = sum(n);
%     n_uns_modes(ii) = n;
%     fprintf('%-7d ',n);
%     max_rate(ii) = max(g_rate);
%     [~, ind] = max(abs(tu_shift));
%     max_shift(ii) = tu_shift(ind);
% end
% fprintf('\n\n');
% rel_tuneshift = max_shift';
% rel_tuneshift0  = real(deltaw0);
% rel_growth    = max_rate'*tau;
% rel_growth0   = imag(deltaw0)*tau;
% 
% 
% %% Plot Results
% 
% % Create figure
% scrsz = get(0,'ScreenSize');
% figure1 = figure('OuterPosition',[scrsz(1)+100 scrsz(2)+40 scrsz(3)*0.9 scrsz(4)*0.9]);
% 
% h = 0.27;
% v = 0.39;
% hs = 0.055;
% vs = 0.09;
% 
% % Create subplot
% subplot11 = subplot(2,3,1,'Parent',figure1,'FontSize',16,'Position',[hs (2*vs+v) h v]);
% box(subplot11,'on');
% hold(subplot11,'all');
% % Create multiple lines using matrix input to plot
% plot11 = plot(subplot11, chrom,rel_tuneshift);
% % Create xlabel
% xlabel({'Chromaticity'},'FontSize',16);
% % Create ylabel
% ylabel({'Re(\Omega - \nu_b)/\nu_s'},'FontSize',16);
% % Create title
% title({'Tune Shift of the most shifted coupled bunch mode'},'FontSize',16);
% 
% 
% % Create subplot
% subplot21 = subplot(2,3,4,'Parent',figure1,'FontSize',16,'Position',[hs vs*2/3 h v]);
% box(subplot21,'on');
% hold(subplot21,'all');
% % Create multiple lines using matrix input to plot
% plot21 = plot(subplot21, (1:nb)-1,rel_tuneshift0);
% xlim([0 (nb-1)]);
% % Create xlabel
% xlabel({'Coupled Bunch Mode'},'FontSize',16);
% % Create ylabel
% ylabel({'Re(\Omega - \nu_b)/\nu_s'},'FontSize',16);
% % Create title
% title({'Tune shifts @ zero chromaticity and m=0'},'FontSize',16);
% 
% 
% % Create subplot
% subplot13 = subplot(2,3,3,'Parent',figure1,'FontSize',16,'Position',[(3*hs+2*h) (2*vs+v) h v]);
% box(subplot13,'on');
% hold(subplot13,'all');
% % Create multiple lines using matrix input to plot
% plot(subplot13, chrom,rel_growth);
% % Create xlabel
% xlabel({'Chromaticity'},'FontSize',16);
% % Create ylabel
% ylabel({'\tau_{damp}/\tau_g'},'FontSize',16);
% % Create title
% title({'Growth Rates of the most unstable coupled bunch mode'},'FontSize',16);
% % Create legend
% 
% % Create subplot
% subplot23 = subplot(2,3,6,'Parent',figure1,'FontSize',16,'Position',[(3*hs+2*h) vs*2/3 h v]);
% box(subplot23,'on');
% hold(subplot23,'all');
% % Create multiple lines using matrix input to plot
% plot23 = plot(subplot23, (1:nb)-1,rel_growth0);
% xlim([0 (nb-1)]);
% % Create xlabel
% xlabel({'Coupled Bunch Mode'},'FontSize',16);
% % Create ylabel
% ylabel({'\tau_{damp}/\tau_g'},'FontSize',16);
% % Create title
% title({'Growth rates @ zero chromaticity and m=0'},'FontSize',16);
% 
% % Create subplot
% subplot12 = subplot(2,3,2,'Parent',figure1,'FontSize',16,'Position',[(2*hs+h) (2*vs+v) h v]);
% box(subplot12,'on');
% hold(subplot12,'all');
% % Create multiple lines using matrix input to plot
% plot12 = plot(subplot12, chrom,n_uns_modes');
% % Create xlabel
% xlabel({'Chromaticity'},'FontSize',16);
% % Create ylabel
% ylabel({'#'},'FontSize',16);
% % Create title
% title({'Number of Unstable Coupled Bunch Modes'},'FontSize',16);
% 
% % Create textbox
% 
% if strcmp(plane,'v')
%     string = [{'Plane: Vertical'}, {' '}]; 
% else
%     string = [{'Plane: Horizontal'}, {' '}]; 
% end
% string = [string , {'Impedances used:'}];
% for i=imped
%     string = [string, {sprintf('- %s', budget{i}.name)}];
% end
% string = [string , {' '}, {'Stage of the Machine:'}, {ringdata.stage}];
% 
% annotation(figure1,'textbox',...
%     [(hs*3/2+h) vs h*3/5 v],...
%     'String',string,...
%     'FontSize',16,...
%     'FitBoxToText','off');
% 
% if save
%     saveas(figure1,['transverse_cbi_' plane '_' ringdata.stage '.fig']);
% end
% 
% 

%% Calc coherent tune-shifts

fprintf('Calculation of %s Coupled Bunch Instability and Tune Shifts\n', label);
fprintf('%-20s: %-20.4g\n','Betatron Tune', nut);
fprintf('%-20s: %-20.4g\n\n','Damping Time [ms]', tau*1e3);
fprintf('%-7s # of unstable Modes Modes\n','Chrom');
fprintf('%-7s ',' ');fprintf('%-5d ', azi_modes);fprintf('\n');
fprintf('%s',repmat('-',1,40)); fprintf('\n');

imtune_shift = [];
retune_shift = [];
imidx = [];
reidx = [];
n_uns_modes = [];
deltaw = zeros(length(azi_modes),nb,length(chrom));
for i=1:length(chrom)
    g_rate = [];
    ig_rate = [];
    tu_shift = [];
    itu_shift = [];
    n_uns_mode = [];
    fprintf('%-5.3g:  ',chrom(i));
    for m=1:length(azi_modes)
        deltaw(m,:,i) = lnls_calc_transverse_cbi(w,Zt, sigma, nb, w0, nus,...
              nut, chrom(i), eta, azi_modes(m), E, I_tot);
        [gmax igmax] = max(imag(deltaw(m,:,i)));
        if chrom(i)==0 && azi_modes(m)==0;
            deltaw0 = deltaw(m,:,i);
        end
        g_rate = [g_rate, gmax];
        ig_rate = [ig_rate, igmax];
        [tumax itumax] = max(abs(real(deltaw(m,:,i))));
        tu_shift = [tu_shift, real(deltaw(m,itumax,i))];
        itu_shift = [itu_shift, itumax];
        n = imag(deltaw(m,:,i)) > 1/tau ;
        n = sum(n);
        n_uns_mode = [n_uns_mode, n];
        fprintf('%-5d ',n);
    end
    fprintf('\n');
    imtune_shift = [imtune_shift; g_rate];
    imidx        = [imidx; ig_rate];
    retune_shift = [retune_shift; tu_shift];
    reidx        = [reidx; itu_shift];
    n_uns_modes  = [n_uns_modes; n_uns_mode];
end
fprintf('\n\n');
rel_tuneshift = retune_shift'/w0/nus;
rel_tuneshift0  = real(deltaw0)/w0/nus;
rel_growth    = imtune_shift';
rel_growth0   = imag(deltaw0);


%% Plot Results

f1 = figure('Position',[1 1 850 528]);

axes1 = axes('Parent',f1,'Position',[0.12 0.103 0.850 0.415], 'FontSize',16);
box(axes1,'on'); hold(axes1,'all'); grid(axes1,'on');

plot(chrom,rel_tuneshift,'Parent',axes1,'LineWidth',2);
xlabel({'Norm. Chromaticity'},'FontSize',16);
ylabel({'Re(\Omega - \nu_{\beta})/\nu_s'},'FontSize',16);
xlim([min(chrom), max(chrom)]);


axes2 = axes('Parent',f1,'Position',[0.12 0.518 0.850 0.415],...
             'FontSize',16, 'XTickLabel',{''});
box(axes2,'on');hold(axes2,'all');grid(axes2,'on');
plot2 = plot(axes2, chrom,rel_growth,'LineWidth',2);
plot(axes2, chrom,repmat(1/tau,size(chrom)),'LineWidth',2,...
     'LineStyle','--','Color','k','DisplayName','1/\tau_d');
for i=1:length(azi_modes)
    leg = sprintf('m = %d',azi_modes(i));
    set(plot2(i),'DisplayName',leg);
end
legend(axes2,'show','Location','Best');
ylabel({'1/\tau [Hz]'},'FontSize',16);
xlim([min(chrom), max(chrom)]);


f2 = figure('Position',[1 1 850 528]);

axes3 = axes('Parent',f2,'Position',[0.12 0.103 0.850 0.415], 'FontSize',16);
box(axes3,'on'); hold(axes3,'all'); grid(axes3,'on');

plot((1:nb)-1,rel_tuneshift0,'Parent',axes3,'LineWidth',2);
xlabel({'Coupled Bunch Mode'},'FontSize',16);
ylabel({'Re(\Omega - \nu_{\beta})/\nu_s'},'FontSize',16);
xlim([0 (nb-1)]);


axes4 = axes('Parent',f2,'Position',[0.12 0.518 0.850 0.415],...
             'FontSize',16, 'XTickLabel',{''});
box(axes4,'on');hold(axes4,'all');grid(axes4,'on');
plot(axes4, (1:nb)-1,rel_growth0,'LineWidth',2,'DisplayName','Chrom=0, m=0');
plot(axes4, (1:nb)-1,repmat(1/tau,1,nb),'LineWidth',2,...
     'LineStyle','--','Color','k','DisplayName','1/\tau_d');
legend(axes4,'show','Location','Best');
ylabel({'1/\tau_{damp} [Hz]'},'FontSize',16);
xlim([0 (nb-1)]);


if save
    saveas(f1,['cbi_',plane,'_chrom_', ringdata.stage '.fig']);
    saveas(f2,['cbi_',plane,'_mode_',  ringdata.stage '.fig']);
end



