function budget = booster_impedance_budget(w, select)



mu0 = 4*pi*1e-7;
c   = 299792458;
ep0 = 1/c^2/mu0;
E = 3;

i=1;
%% Resistive wall form cilindrical vaccum chamber;
if (any(strcmp(select,'resistive_wall_sc')) || strcmp(select,'all') || strcmp(select,'ring'))
    budget{i}.name = 'Resistive Wall SS';
    budget{i}.type = 'rw';
    budget{i}.quantity = 1;
    budget{i}.betax = 10;
    budget{i}.betay = 12;
    epb     = [1 1];
    mub     = [1 1];
    ange    = [0 0];
    angm    = [0 0];
    sigmadc = [0 1/7.4e-7];
    tau     = [0 1]*16e-15;
    b       = 18.000*1e-3;
    L       = 440;
    budget{i}.mub = mub;
    budget{i}.ange = ange;
    budget{i}.angm = angm;
    budget{i}.tau   = tau;
    budget{i}.sigmadc = sigmadc;
    budget{i}.epb = epb;
    budget{i}.b = b;
    budget{i}.L = L;
    
    for j = 1: length(epb)
        epr(j,:) = epb(j)*(1-1i.*sign(w).*tan(ange(j))) + sigmadc(j)./(1+1i*w*tau(j))./(1i*w*ep0);
        mur(j,:) = mub(j)*(1-1i.*sign(w).*tan(angm(j)));
    end
    [Zl Zv Zh] = lnls_calc_impedance_multilayer_round_pipe(w, epr, mur, b, L, E, false, 0, 0);
    budget{i}.Zv = Zv;
    budget{i}.Zh = Zh;
    budget{i}.Zl = Zl;
    budget{i}.escala = 'log';
    i=i+1;
end

if (any(strcmp(select,'resistive_wall_dip')) || strcmp(select,'all') || strcmp(select,'ring'))
    budget{i}.name = 'Resistive Wall Bends';
    budget{i}.type = 'rw';
    budget{i}.quantity = 1;
    budget{i}.betax = 1.5;
    budget{i}.betay = 22;
    epb     = [1 1];
    mub     = [1 1];
    ange    = [0 0];
    angm    = [0 0];
    sigmadc = [0 1/7.4e-7];
    tau     = [0 1]*16e-15;
    b       = 11.700*1e-3;
    L       = 58;
    budget{i}.mub = mub;
    budget{i}.ange = ange;
    budget{i}.angm = angm;
    budget{i}.tau   = tau;
    budget{i}.sigmadc = sigmadc;
    budget{i}.epb = epb;
    budget{i}.b = b;
    budget{i}.L = L;
    for j = 1: length(epb)
        epr(j,:) = epb(j)*(1-1i.*sign(w).*tan(ange(j))) + sigmadc(j)./(1+1i*w*tau(j))./(1i*w*ep0);
        mur(j,:) = mub(j)*(1-1i.*sign(w).*tan(angm(j)));
    end
    [Zl Zv Zh] = lnls_calc_impedance_multilayer_round_pipe(w, epr, mur, b, L,E, false, 0, 0);
    budget{i}.Zv = Zv;
    budget{i}.Zh = Zh;
    budget{i}.Zl = Zl;
    budget{i}.escala = 'log';
    i=i+1;
end

%% Ferrite Kickers for injection and extraction
if (any(strcmp(select,'kicker')) || strcmp(select,'all') || strcmp(select,'ring') )
    budget{i}.name = 'Ferrite Kickers';
    budget{i}.type = 'misto';
    budget{i}.quantity = 1;
    budget{i}.betax = 23;
    budget{i}.betay = 4;
% Valores que peguei com o F??bio.
    a = 63/2*1e-3;
    b = (10+7.5)*1e-3;
    d = b + 20e-3;
    L = 0.5 + 0.2;  % extracao e injecao
    % Ferrite CMD5005
    epl = 12;
    rho = 1e6;
    epr = epl - 1i./w/ep0/rho*0;
    mui = 1600;
    ws = 20e6;
    mur = 1 + mui./(1+1i*w/ws);
    Zg = (1/50 + 1i*w*30e-12).^-1;
    % mean case
    [Zlw Zhw Zvw] = lnls_calc_ferrite_kicker_impedance(w,a,b,d,L,epr,mur,Zg,'pior');
    [Zlb Zhb Zvb] = lnls_calc_ferrite_kicker_impedance(w,a,b,d,L,epr,mur,Zg,'melhor');
    theta = atan(b/a);
    Zl = (theta*Zlb + (pi/2-theta)*Zlw)/(pi/2);
    Zv = (theta*Zvb + (pi/2-theta)*Zvw)/(pi/2);
    Zh = (theta*Zhb + (pi/2-theta)*Zhw)/(pi/2);
        
    budget{i}.L = L;
    budget{i}.Zv = Zv;
    budget{i}.Zh = Zh;
    budget{i}.Zl = Zl;
    budget{i}.escala = 'log';
    i=i+1;
end

%% Bellows;
if (any(strcmp(select,'bellows')) || strcmp(select,'all') || strcmp(select,'ring'))
    budget{i}.name = 'Bellows';
    budget{i}.type = 'geo';
    budget{i}.quantity = 100;
    budget{i}.betax = 5;
    budget{i}.betay = 16;
    R = 17.2e-3;
    L = 1.156e-3;
    t = 6.8e-3;
    comprimento = 23e-3/L;
    
    [~, Zv, Zh, ~] = lnls_calc_impedance_bellow(w,R,L,t, comprimento);
    [Zlor wor] = lnls_load_impedance(...
        '/home/fernando/MATLAB/Impedancias/Simulacoes/bellows_booster/', ...
        'Longitudinal','ECHO');
    Zl = interp1(wor,Zlor,w,'linear',0);
    ind = real(Zl)<0;
    Zl(ind) = 1i*imag(Zl(ind));
   
    budget{i}.Zv = Zv;
    budget{i}.Zh = Zh;
    budget{i}.Zl = Zl;
    budget{i}.escala = 'linear';
    i=i+1;
end

%% BPMs;
if (any(strcmp(select,'bpm')) || strcmp(select,'all') || strcmp(select,'ring'))
    budget{i}.name = 'BPM-3-BNcernew';
    budget{i}.type = 'geo';
    budget{i}.quantity = 50;
    budget{i}.betax = 10;
    budget{i}.betay = 10;
    % Fiz uma mescla entre impedancia fittada por ressonadores e a impedancia
    % interpolada a partir da original. Usei a fittada para reproduzir o
    % broad-band e a interpolada para os picos narrow-band.
    Zlf = lnls_load_fitted_impedance(w, ...
        '/home/fernando/MATLAB/Impedancias/Simulacoes/BPM/004-BPM_buttons/3/BNcernew', ...
        'Longitudinal','GdfidL');
    [Zlor wor] = lnls_load_impedance(...
        '/home/fernando/MATLAB/Impedancias/Simulacoes/BPM/004-BPM_buttons/3/BNcernew', ...
        'Longitudinal','GdfidL');
    Zlint = interp1(wor,Zlor,w,'linear',0);
    Zl = Zlint;
    ind = (Zlint ==0);
    Zl(ind) = Zlf(ind);
    % Zl = Zlf;
    
    % Zxf = lnls_load_fitted_impedance(w, ...
    %     '/home/fernando/MATLAB/Impedancias/Simulacoes/Masks/_results/Ridge/w10/SOFT_HARD/h2', ...
    %     'Horizontal','GdfidL');
    % [Zxor wor] = lnls_load_impedance(...
    %     '/home/fernando/MATLAB/Impedancias/Simulacoes/Masks/_results/Ridge/w10/SOFT_HARD/h2', ...
    %     'Horizontal','GdfidL');
    % Zxint = interp1(wor,Zxor,w,'linear',0);
    % ind = (Zxint ==0) | abs(w)<1e10;
    % Zxint(ind) = Zxf(ind);
    
    budget{i}.Zv = Zl*0;
    budget{i}.Zh = Zl*0;
    budget{i}.Zl = Zl;
    budget{i}.escala = 'linear';
    i=i+1;
end

%% Broad-band impedance;

if (any(strcmp(select,'broad_band')) || strcmp(select,'all') || strcmp(select,'ring'))
    budget{i}.name = 'Broad Band';
    budget{i}.type = 'geo';
    budget{i}.quantity = 1;
    budget{i}.betax = 9;
    budget{i}.betay = 12;                    
    budget{i}.Rsx =  1.43*1e6*518.25/754.0/budget{i}.betax;
    budget{i}.wrx =  22*1e9*2*pi;
    budget{i}.Qx =   1;                        
    budget{i}.Rsy =  1.43*1e6*518.25/754.0/budget{i}.betay;
    budget{i}.wry =  22*1e9*2*pi;              
    budget{i}.Qy =   1;                        
    budget{i}.Rsl = 3.6*518.25/354.0*1e3;
    budget{i}.wrl = 20*1e9*2*pi;
    budget{i}.Ql =  1;
    Zv = lnls_calc_impedance_transverse_resonator(budget{i}.Rsy, budget{i}.Qy, budget{i}.wry, w);
    Zh = lnls_calc_impedance_transverse_resonator(budget{i}.Rsx, budget{i}.Qx, budget{i}.wrx, w);
    Zl = lnls_calc_impedance_longitudinal_resonator(budget{i}.Rsl, budget{i}.Ql, budget{i}.wrl, w);
    
    budget{i}.Zv = Zv;
    budget{i}.Zh = Zh;
    budget{i}.Zl = Zl;
    budget{i}.escala = 'linear';
    i=i+1;
end


%% Petra-5-cell cavity;

if (any(strcmp(select,'petra_cav')) || strcmp(select,'all') || strcmp(select,'ring'))
    budget{i}.name = 'Petra Cavity';
    budget{i}.type = 'geo';
    budget{i}.quantity = 1;
    budget{i}.betax = 9;
    budget{i}.betay = 12; 
    
    data = importdata('/home/fernando/MATLAB/Impedancias/Simulacoes/RF_cav/_results/Petra5c/200_first_homs_petra.dat');
    data = data.data;

    budget{i}.Rsl =  data(:,3);
    budget{i}.wrl =  (data(:,2)*1e9 + 000*1e3)*2*pi;
    budget{i}.Ql =   data(:,4);                        
    Zl = lnls_calc_impedance_longitudinal_resonator(budget{i}.Rsl, budget{i}.Ql, budget{i}.wrl, w);
    
    budget{i}.Zv = Zl*0;
    budget{i}.Zh = Zl*0;
    budget{i}.Zl = Zl;
    budget{i}.escala = 'linear';
    i=i+1;
end