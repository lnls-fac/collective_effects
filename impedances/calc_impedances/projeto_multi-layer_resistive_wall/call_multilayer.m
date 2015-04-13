clear all

mu0 = 4*pi*1e-7;
c   = 299792458;
ep0 = 1/c^2/mu0;

w = logspace(11,14,10000);
% w = [-fliplr(w) w];

epb     = [1 1];
mub     = [1 1];
ange    = [0 0];
angm    = [0 0];
sigmadc = [0 3.344e7];
tau     = [0 0]*27e-15;
b       = [5]*1e-2;
L       = 1;

for j = 1: length(epb)
        epr(j,:) = epb(j)*(1-1i.*sign(w).*tan(ange(j))) + sigmadc(j)./(1+1i*w*tau(j))./(1i*w*ep0);
        mur(j,:) = mub(j)*(1-1i.*sign(w).*tan(angm(j)));
end

[Zl Zv Zh] = lnls_calc_impedance_multilayer_round_pipe(w, epr, mur, b, L, 3, false, 0, 0);

budget{1}.name = 'Mounet';
budget{1}.w = w;
budget{1}.Zh = Zh;
budget{1}.Zv = Zv;
budget{1}.Zl = Zl;
budget{1}.quantity = 1;


figure1 = figure('XVisual',...
    '0x21 (TrueColor, depth 24, RGB mask 0xff0000 0xff00 0x00ff)');

% Create axes
axes1 = subplot(2,1,2,'Parent',figure1,'XScale','log','XMinorTick','on',...
    'FontSize',16);
box(axes1,'on');
hold(axes1,'all');

% Create multiple lines using matrix input to semilogx
semilogx1 = semilogx(w,[real(Zh);imag(Zh)],'Parent',axes1);
set(semilogx1(1),'DisplayName','Real');
set(semilogx1(2),'DisplayName','Imaginária');

% Create xlabel
xlabel('\omega [rad/s]','FontSize',16);

% Create ylabel
ylabel('Z_t/L [Ohm/m^2]','FontSize',16);

% Create legend
legend(axes1,'show');


% Create axes
axes2 = subplot(2,1,1,'Parent',figure1,'XScale','log','XMinorTick','on',...
    'FontSize',16);
box(axes2,'on');
hold(axes2,'all');

% Create multiple lines using matrix input to semilogx
semilogx2 = semilogx(w,[real(Zl);imag(Zl)],'Parent',axes2);
set(semilogx2(1),'DisplayName','Real');
set(semilogx2(2),'DisplayName','Imaginária');

% Create xlabel
xlabel('\omega [rad/s]','FontSize',16);

% Create ylabel
ylabel('Z_l/L [Ohm/m]','FontSize',16);

% Create legend
legend(axes2,'show');

return;


%% Comparison with Chao's formulae
cond = sigmadc(2);
thick = 2e-3;
perm = 1;
[Zt Zl w_val] = lnls_calc_impedance_resistive_wall(thick, ...
        cond, perm, b, L, w);

budget{2}.name = 'Chao';
budget{2}.Zh = Zt;
budget{2}.Zv = Zt;
budget{2}.Zl = Zl;
budget{2}.quantity = 1;


plot_impedances(budget,'log','linear',false);