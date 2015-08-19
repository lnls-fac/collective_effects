function  plot_results(globdata)
% function  plot_results(globdata)
% Function for plotting results requested by the user

% Data info
dsrc  = globdata.simpar.datasource;
m     = globdata.simpar.m;
waxis = globdata.simpar.whichaxis;
if strcmp(waxis,'y')
    wplane = 'Vertical';
elseif strcmp(waxis,'x')
    wplane = 'Horizontal';
end
taxis = globdata.simpar.whichaxis;

% Wakepotential
w     = globdata.results.W;
s     = globdata.results.s;
sigs  = globdata.simpar.sigma;

% Impedance
if m==0
    rez = globdata.results.ReZlong;
    imz = globdata.results.ImZlong;
elseif m>0
    rez = globdata.results.ReZt/1000;
    imz = globdata.results.ImZt/1000;
end
ImZoN = globdata.results.ImZoN;
f     = globdata.results.freq;
naxis = globdata.results.naxis;

% Loss / Kick Factor
sigi   = globdata.results.sigmak;
kZi    = globdata.results.klossZ;
kW     = globdata.results.klossW;
kickZi = globdata.results.kickZ;
kickW  = globdata.results.kickW;


%% Tick Position # 0: Plot wakepotential

%% Short Range
figure('Units','characters','Position',[15 37 80 30]);

%========= Plot bunch shape =========
sbunch = linspace(-5*sigs,5*sigs,1000); % 5 sigma
ss = sbunch.^2;
norm = max(w);                          % Normalization factor
bunchshape = norm*exp(-ss./(2*sigs^2));

plot(sbunch*1000,bunchshape,'b','LineWidth',2);

hold on
%========= Plot wakepotential =========
plot(s*1000,w,'r','LineWidth',2);
grid on;

%========= Longitudinal Options =========
if m==0
    title (['Longitudinal Wakepotential (' dsrc ')'],'fontsize',13)
    xlabel('s (mm)','fontsize',13);
    ylabel('W (V)','fontsize',13);
    %========= Transverse Options =========
elseif m==1
    title ([wplane ' Dipole Wakepotential (' dsrc ')'],'fontsize',13)
    xlabel('s (mm)','fontsize',13);
    ylabel(['W_D_' waxis '(V/m)'],'fontsize',13);
elseif m==2
    title ([wplane ' Quadrupole Wakepotential (' dsrc ')'],'fontsize',13)
    xlabel('s (mm)','fontsize',13);
    ylabel(['W_Q_' waxis '(V/m)'],'fontsize',13);
end

set(gca,'FontSize',12)
xlim( [s(1)*1000 7000*sigs])
ylim ([-inf norm*1.5])
legend ('Bunch Shape','Wakepotential')

%% Long Range

figure('Units','characters','Position',[15 10 150 20]);

plot(s,w,'r','LineWidth',2);
grid on;

%========= Longitudinal Options =========
if m==0
    title (['Longitudinal Wakepotential (' dsrc ')'],'fontsize',13)
    xlabel('s (m)','fontsize',13);
    ylabel('W (V)','fontsize',13);
    %========= Transverse Options =========
elseif m==1
    title ([wplane ' Dipole Wakepotential (' dsrc ')'],'fontsize',13)
    xlabel('s (m)','fontsize',13);
    ylabel(['W_D_' waxis '(V/m)'],'fontsize',13);
elseif m==2
    title ([wplane ' Quadrupole Wakepotential (' dsrc ')'],'fontsize',13)
    xlabel('s (m)','fontsize',13);
    ylabel(['W_Q_' waxis '(V/m)'],'fontsize',13);
end

%% Tick Position # 3: Plot Impedance
figure('Units','characters','Position',[100 37 140 25]);

plot(f/1e9,rez,'r','LineWidth',2);
hold on
plot(f/1e9,imz,'b--','LineWidth',2);

xlabel('Frequency, GHz','fontsize',13)
set(gca,'FontSize',12);

%========= Longitudinal Options =========
if m==0
    title(['Longitudinal Impedance (' dsrc ')'],'fontsize',13)
    ylabel('Z_{||}, \Omega','fontsize',13)
    %========= Transverse Options =========
elseif m==1
    title ([wplane ' Dipole Impedance (' dsrc ')'],'fontsize',13)
    ylabel(['Z_D_' waxis ', k\Omega/m'],'fontsize',13);
elseif m==2
    title ([wplane ' Quadrupole Impedance (' dsrc ')'],'fontsize',13)
    ylabel(['Z_Q_' waxis ', k\Omega/m'],'fontsize',13);
end
grid
legend ('Re', 'Im')
xlim([0 f(end)]/1e9)
%             ylim([-0.01 1.1*max(rezlong)])
%  ylim([-0.01 40])


%       ========= Plot ImZ/n =========
%         if m==0
%             figure('Units','characters','Position',[245 43 100 30]);
%             plot(naxis,ImZoN,'b','LineWidth',2);
%
%             grid minor;
%             title(['Im(Z_|_|) / n (' dsrc ')'],'fontsize',13)
%             xlabel('Normalized Frequency ( \omega/\omega_0 )','fontsize',13)
%             ylabel('Im(Z) / n, \Omega','fontsize',13)
%             set(gca,'FontSize',12);
%             xlim([0 naxis(end)])
%     %         ylim([-0.1 1.1*max(rezlong)])
%         end


%% Tick Position # 6: Plot Loss/Kick Factor vs. Sigma
if m==0
    
    %========= Plot k_lossZ x sigma =========
    figure('Units','characters','Position',[170 10 120 20]);
    plot(sigi * 1e3, kZi * 1e3, 'o','MarkerSize',2)
    hold on
    plot(sigs * 1e3, kW * 1e3, '*','MarkerSize',5,'LineWidth',2,'color',[1 0 0])
    xlabel('\fontsize{12} \sigma, mm')
    ylabel('\fontsize{12} \kappa_{loss}, mV/pC')
    legend('\kappa_{loss}^Z','\kappa_{loss}^W')
    set(gca,'FontSize',12);
    grid on
    text(sigs * 1.1e3, kW * 1e3, ['\kappa_{loss}^W= ' num2str(kW*1e3,2) ' mV/pC'],'fontsize',12,'fontweight','bold');
    
elseif m > 0
    if m==1
        subind = 'D';
    else
        subind = 'Q';
    end
    
    %========= Plot kickZ x sigma =========
    figure('Units','characters','Position',[170 13 120 20]);
    plot(sigi * 1e3, kickZi, 'o','MarkerSize',2)
    hold on
    plot(sigs * 1e3, kickW, '*','MarkerSize',5,'LineWidth',2,'color',[1 0 0])
    xlabel('\fontsize{12} \sigma, mm')
    ylabel(['\fontsize{12} \kappa_{' subind waxis '}, V/pC/m'])
    legend(['\kappa_{' subind waxis '}^Z'],['\kappa_{' subind waxis '}^W'])
    set(gca,'FontSize',12);
    grid on
    text(sigs * 1.1e3, kickW, ['\kappa_{' subind waxis '}^W= ' num2str(kickW,2) ' V/pC/m'],'fontsize',12,'fontweight','bold');
end

%% Tick Position # 18: Plot CBM Analysis
%
% % CBM analysis
% whalf    = globdata.results.peakinfo.wmodel;
% Zbhalf   = globdata.results.peakinfo.ReZmodel;
% Rshs     = globdata.results.peakinfo.Rshunt;
% GRs      = globdata.results.GRs;
% GR_HOM   = globdata.results.GR_HOM;
% ifast    = globdata.results.ifast;
% ReZsampl = globdata.results.ReZsampl;
% fsampl   = globdata.results.fsampl;
% omegar   = globdata.results.peakinfo.omegar;
% iCBM     = 0:1:globdata.ringpar.h-1;
% p        = -100:1:100;
% units    = globdata.simpar.units;

%% ========= Plot Impedance Spectra =========
%
%     scrsz = get(0,'ScreenSize');
%     figure('Position',[scrsz(4)/2 50 scrsz(3)/2 scrsz(4)-150]);
%     subplot(3,1,1);
%     plot(whalf/(2e9*pi),Zbhalf,'LineWidth',2,'Color','blue');
%     hold on;
%
%% ========= plot growth rate for each mode =========
%     subplot(3,1,2);
%     stem(iCBM,GRs);
%     grid;
%     xlabel('Coupled Bunch Mode','fontsize',13);
%     ylabel('Growth Rate (1/s)','fontsize',13);
%     % eixox = 1:9:M;
%     % set(gca,'XTick',eixox);
%
%% =========  fastest growth time/rate =========
%     g_rate = GRs(ifast)*units;
%     g_time = 1/GRs(ifast)/units;
%
%     disp(['GR (' num2str(units) ' units) = ' num2str(g_rate) ' 1/s']);
%     disp(['GT (' num2str(units) ' units) = ' num2str(g_time*1000) ' ms']);
%
%
%% ========= plots sampled impedance =========
%     subplot(3,1,1);
%     stem(fsampl(1:floor(length(p)/2)),ReZsampl(1:floor(length(p)/2)),'Color','black','LineWidth',2,'Marker','none','MarkerSize',5);
%     stem(fsampl(floor(length(p)/2)+1:end),ReZsampl(floor(length(p)/2)+1:end),'Color','red','LineWidth',2,'Marker','none','MarkerSize',5);
%     grid;
%     xlim([0 1.2*max(omegar)/(2e9*pi)]);
%     legend('Re[Z_|_|(\omega)]',['Negative CBM - ', num2str(ifast-1)],['Positive CBM - ', num2str(ifast-1)]);
%     xlabel('Frequency (GHz)','fontsize',13);
%     ylabel('R_|_| (\Omega)','fontsize',13);
%     title('Longitudinal CBM Threshold Analysis','fontsize',14);
%
%
%% ========= plots accumulated growth rate across HOMs =========
%     m = 1:1:length(Rshs)+1;
%     subplot(3,1,3);
%     stairs(m,[GR_HOM(ifast,:) GR_HOM(ifast,end)],'b','LineWidth',2);
%     grid;
%     legend(['Fastest Mode ', num2str(ifast-1)],'Location','NorthEast');
%     xlabel('Higher Order Mode','fontsize',13);
%     ylabel('Acc. Growth Rate (1/s)','fontsize',13);
%     xlim([1 m(end)])
