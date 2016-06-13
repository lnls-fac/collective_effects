function varargout = lnls_characterize_lmci(ringdata, budget, params, save)
%function varargout = lnls_characterize_lmci(ringdata, budget, params, save)

I     = params.I;
n_azi = params.n_azi; 
n_rad = params.n_rad;

w = ringdata.w;

%% Load parameters

Zls    = 'Zl';
tau  = ringdata.taue;

imped  = 1:length(budget);
Zl = zeros(size(w));
for i=imped
    quant = budget{i}.quantity;
    Zl   = Zl + quant*budget{i}.(Zls);
end

stage    = ringdata.stage;
w0 = ringdata.w0;       % revolution angular frequency [Hz]
nus = ringdata.nus;     % synchrotron tune
eta = ringdata.eta;     % momentum compaction factor
E = ringdata.E;         % energy [GeV];

%% Calcula Longitudinal mode Coupling
fprintf('Calculation of Longitudinal Mode Coupling Instability\n');
fprintf('%-20s: %-20.4g\n','Shynchrotron Tune', nus);
fprintf('%-20s: %-20d\n','Azimuthal Modes', n_azi);
fprintf('%-20s: %-20d\n','Radial Modes', n_rad);
fprintf('%-20s: %-20.4g\n\n','Damping Time [ms]', tau*1e3);
fprintf('I [mA]: '); fprintf('%-5.3g ', I*1e3);fprintf('\n');
fprintf('Stable? ');

params.E     = E;   
params.w0    = w0;   
params.nus   = nus;
params.eta   = eta;   

delta = lnls_calc_longitudinal_mode_couplingopt(w, Zl, params);

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
fprintf('\n\n');

varargout{2} = delta;
%% Plot Results

[real_delta, ~] = sort(real(delta));
[imag_delta, ~] = sort(imag(delta));




f1 = figure('Position',[1 1 850 528]);

axes1 = axes('Parent',f1,'Position',[0.12 0.103 0.850 0.415], 'FontSize',16);
box(axes1,'on'); hold(axes1,'all'); grid(axes1,'on');
plot(I*1e3, real_delta,'Parent',axes1,'LineWidth',2,'Color','b');
xlabel({'Current per bunch [mA]'},'FontSize',16)
ylabel({'Re(\Omega/\omega_s)'},'FontSize',16);
xlim([min(I), max(I)]*1e3);


axes2 = axes('Parent',f1,'Position',[0.12 0.518 0.850 0.415],...
             'FontSize',16, 'XTickLabel',{''});
box(axes2,'on');hold(axes2,'all');grid(axes2,'on');
plot(axes2, I*1e3, imag_delta*nus*w0*tau,'LineWidth',2,'Color','b');
ylabel({'\tau_{damp}/\tau_g'},'FontSize',16);
xlim([min(I), max(I)]*1e3);

if save
    saveas(f1,['lmci_' ringdata.stage '.fig']);
end
