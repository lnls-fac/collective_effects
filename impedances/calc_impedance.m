%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                    %
%                This script extracts impedance from Wakepotential data              %
%                                                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c = 299792458;

%% Extracts Needed Variables

path = '/home/ABTLUS/fernando.sa/MATLAB/Impedancias/Simulacoes/bellows_booster/';

[~, wakes] = hdrload([path 'wake.dat']);
m = 0;
sigs = 3e-3;
wake = wakes(:,2)'*1e12; % rescale W to [V/C]
saxis = wakes(:,1)'*1e-2;
f0 = c/496.8;

figure1 = figure('Units','normalized','Position',[0.109, 0.331,0.507,0.431]);
axes1 = axes('Parent',figure1,'YGrid','on','XGrid','on','FontSize',16);
xlim(axes1,[-0.02 1]);
box(axes1,'on');
hold(axes1,'all');
plot(saxis,wake/1e12);
xlabel('s [m]','FontSize',16);
ylabel('Wake Potential [V/pC]','FontSize',16);

sigt = sigs/c;

%% Impedance calculation based on Wakepotential

% initializes vectors

freq = zeros(1,length(saxis));
time = freq;

% builds time and frequency scales

% translates scale to zero origin
wakeorigin = saxis-min(saxis);

% frequency scale (Hz)
Np = length(saxis);
ds = max(wakeorigin)/(Np-1);
dt = ds/c;

% Scale shift/zero shift

shift = (saxis(1))/c;

% calculates FFT and frequency scale
NFFT = Np;%2^nextpow2(Np);
%         if strcmp(wake_type,'T')
%             ffttt = -1i*fft(-wake(:,2),NFFT)';
%         else
ffttt = fft(-wake',NFFT)';
%         end
p = 1:1:NFFT;
omegalist = 2*pi*(p-1)*c/(NFFT*ds);
VHat = (dt*exp(1i*omegalist*shift).*ffttt);
freq = omegalist/(2*pi);

% calculates bunch spectrum (current roll-off)

Jwlist = exp(-(omegalist*sigt/sqrt(2)).^2);

Z = VHat./Jwlist;

%% Exports calculated FFT - there is a 10 MHz difference between the peaks, which corresponds to 1/2 df
data = [freq(1:round(NFFT/2)+1)'/1e9 real(Z(1:round(NFFT/2)+1))' imag(Z(1:round(NFFT/2)+1))']';

fmax = c/(pi*sigs*1e9); %defines the frequency range due to the bunch length
ifmax = round(fmax*1e9/freq(2))+1;
lines2save = ifmax;

if m==0
    
    freq = data(1,1:lines2save);                           % assignment to structure
    ReZlong = data(2,1:lines2save);                        % assignment to structure
    ImZlong = -data(3,1:lines2save);                       % assignment to structure
    naxis = freq(2:round(lines2save/20)).*1e9./f0;    % assignment to structure
    ImZoN = ImZlong(2:round(lines2save/20))./naxis; % assignment to structure
    
    R = 17.2e-3;
    L = 1.156e-3;
    t = 6.8e-3;
    comprimento = 23.8e-3/L;
    [Zl, Zv, Zh, Zl2] = lnls_calc_impedance_bellow(2*pi*freq*1e9,R,L,t, comprimento);
    
   
    figure1 = figure('Units','normalized','Position', [0.073,0.36,0.693,0.478]);
        
    axes1 = axes('Parent',figure1,...
        'Position',[0.0827067669172932 0.11 0.381952323991798 0.815],...
        'FontSize',16);
    box(axes1,'on');
    hold(axes1,'all');
    plot(freq,[-ImZlong; imag(Zl); imag(Zl2)],'Parent',axes1);
    xlabel('frequency [GHz]','FontSize',16);
    ylabel('Im(Zl) [\Omega]','FontSize',16);
    
    axes2 = axes('Parent',figure1,...
        'Position',[0.570340909090909 0.11 0.380035030758715 0.815],...
        'FontSize',16);
    box(axes2,'on');
    hold(axes2,'all');
    plot1 = plot(freq,[ReZlong; real(Zl); real(Zl2)],'Parent',axes2);
    set(plot1(1),'DisplayName','Simulation');
    set(plot1(2),'DisplayName','Zotter');
    set(plot1(3),'DisplayName','Kheifets');
    xlabel('frequency [GHz]','FontSize',16);
    ylabel('Re(Zl) [\Omega]','FontSize',16);
    legend(axes2,'show');
    
    
    fp = fopen([path 'ReZlongECHO.txt'], 'w');
    fprintf(fp,'%10.6f %10.6f\n', [freq; ReZlong]);
    fclose(fp);
    fp = fopen([path 'ImZlongECHO.txt'], 'w');
    fprintf(fp,'%10.6f %10.6f\n', [freq; ImZlong]);
    fclose(fp);
    
elseif m>0
    
    freq = data(1,1:lines2save);                           % assignment to structure
    ReZt = data(3,1:lines2save);                          % assignment to structure
    ImZt = data(2,1:lines2save);                          % assignment to structure
    
    plot(freq,[ReZt, ImZt]);
end



