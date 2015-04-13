function varargout = lnls_characterize_tmci(ringdata, budget, params, save)

plane = params.plane;
n_azi = params.n_azi;
n_rad = params.n_rad;
I     = params.I;
chrom = params.chrom;

w = ringdata.w;
%% Load parameters
if strcmp(plane,'v')
    betas = 'betay'; Zts    = 'Zv';
    tau  = ringdata.tauy;
    nut = ringdata.nuty;
    label = 'Vertical';
else
    betas = 'betax'; Zts    = 'Zh';
    tau  = ringdata.taux;
    nut = ringdata.nutx;
    label = 'Horizontal';
end

Zt = zeros(size(w));
for i=1:length(budget)
    Zt   = Zt + budget{i}.quantity*budget{i}.(Zts)*budget{i}.(betas);
end

w0    = ringdata.w0;      % revolution angular frequency [Hz]
nus   = ringdata.nus;     % synchrotron tune
eta   = ringdata.eta;     % momentum compaction factor
E     = ringdata.E;       % energy [GeV];

%% Calcula Transeverse mode Coupling
fprintf('Calculation of %s Mode Coupling Instability\n', label);
fprintf('%-20s: %-20.4g\n','Betatron Tune', nut);
fprintf('%-20s: %-20.4g\n','Chromaticity', chrom);
fprintf('%-20s: %-20d\n','Azimuthal Modes', n_azi);
fprintf('%-20s: %-20d\n','Radial Modes', n_rad);
fprintf('%-20s: %-20.4g\n\n','Damping Time [ms]', tau*1e3);
fprintf('I [mA]: '); fprintf('%-5.3g ', I*1e3);fprintf('\n');
fprintf('Stable? ');

params.E     = E;   
params.w0    = w0;   
params.nus   = nus;
params.nut   = nut;
params.eta   = eta;   

delta =lnls_calc_transverse_mode_couplingopt(w, Zt, params);

first = true;
for i = 1:length(I)
    if any(imag(delta(:,i))*nus*w0*tau > 1)
        fprintf('%-6s','N');
        if first
            varargout{1} = I(i);
        end
        first = false;
    else
        fprintf('%-6s','Y');
    end
end
fprintf('\n\n\n');

varargout{2} = delta;
%% Plot Results

% [real_delta ind] = sort(real(delta));
% for j = 1:size(delta,2)
%     imag_delta(:,j) = imag(delta(ind(:,j),j));
% end
[real_delta, ~] = sort(real(delta));
[imag_delta, ~] = sort(imag(delta));


f1 = figure('Position',[1 1 850 528]);

axes1 = axes('Parent',f1,'Position',[0.12 0.103 0.850 0.415], 'FontSize',16);
box(axes1,'on'); hold(axes1,'all'); grid(axes1,'on');
plot(I*1e3, real_delta,'Parent',axes1,'LineWidth',2,'Color','b');
xlabel({'Current per bunch [mA]'},'FontSize',16)
ylabel({'Re(\Omega - \omega_\beta)/\omega_s'},'FontSize',16);
xlim([min(I), max(I)]*1e3);


axes2 = axes('Parent',f1,'Position',[0.12 0.518 0.850 0.415],...
             'FontSize',16, 'XTickLabel',{''});
box(axes2,'on');hold(axes2,'all');grid(axes2,'on');
plot(axes2, I*1e3, imag_delta*nus*w0*tau,'LineWidth',2,'Color','b');
ylabel({'\tau_{damp}/\tau_g'},'FontSize',16);
xlim([min(I), max(I)]*1e3);


if save
    saveas(f1,['tmci_' plane '_' ringdata.stage '.fig']);
end
