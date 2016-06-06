mu0 = 4*pi*1e-7;
c   = 299792458;
ep0 = 1/c^2/mu0;
E = 3;
% Broad Band Parameters
Zovern = 0.2;
radius = 12e-3;
fr  = 2.4* c/(radius*2*pi); 
Rs = Zovern*fr*(518.396/c);
Q = 1;
wr = fr*2*pi;
Ql = sqrt(Q.^2 - 1/4);
wrl = wr .* Ql ./ Q;

% Resistive Wall parameters:
% Bane and Sands Formula with Krinsky approximation.
radius = 12e-3;
sigma = 5.9e7;
Z0 = c*4*pi*1e-7;
a = 3/sqrt(2*pi)/4;
p = 2.7;
L = 480;
s0 = (2*radius^2/Z0/sigma)^(1/3);
W0 = c*Z0/4/pi * 2*s0*L/radius^4;

% Mounet Formula:
epb     = [1 1 1];
mub     = [1 1 1];
ange    = [0 0 0];
angm    = [0 0 0];
sigmadc = [0 sigma 1];
relax_time = [0 1 0]*27e-15;
b       = [radius 13.000*1e-3];


% I tried several methods to take the fourier transform of the impedance:
% First I tried to do the inverse fft with several frequency intervals

%% Inverse FFT Test

% w =  [    -w,  fliplr(w)];
% Zl = [ conj(Zl), fliplr(Zl)];
% Zv = [-conj(Zv), fliplr(Zv)];
% Zh = [-conj(Zh), fliplr(Zh)];
% [~,inds] = sort(w);
% Zv(inds) = conv(Zv(inds),fil,'same');

% [wake, tau] = inverse_fft(w, -1i/pi*Zv);

% w2_max = (100000+1/2)*w_max;
% w = logspace(1,14,4000);
% epr = zeros(length(epb),length(w));
% mur = zeros(length(epb),length(w));
% for j = 1: length(epb)
%     epr(j,:) = epb(j)*(1-1i.*sign(w).*tan(ange(j))) + sigmadc(j)./(1+1i*w*relax_time(j))./(1i*w*ep0);
%     mur(j,:) = mub(j)*(1-1i.*sign(w).*tan(angm(j)));
% end
% [Zl, Zv, Zh] = lnls_calc_impedance_multilayer_round_pipe(w, epr, mur, b, L, E);
% 
% figure; loglog(w,[imag(Zv);real(Zv)]);
% return
% w_max = [4e14,1e13,4e10];%,1e8,1e6];
% n_samples = 1000000;
% 
% tau  = zeros(length(w_max),n_samples/2);
% wake = zeros(length(w_max),n_samples/2);
% 
% for i = 1:length(w_max);
%     
%     w_min = w_max(i)/(n_samples/2-1);
%     wi = linspace(w_min/2, w_max(i), n_samples/2);
%     
%     epr = zeros(length(epb),length(wi));
%     mur = zeros(length(epb),length(wi));
%     for j = 1: length(epb)
%         epr(j,:) = epb(j)*(1-1i.*sign(wi).*tan(ange(j))) + sigmadc(j)./(1+1i*wi*relax_time(j))./(1i*wi*ep0);
%         mur(j,:) = mub(j)*(1-1i.*sign(wi).*tan(angm(j)));
%     end
%     [Zli, Zvi, Zhi] = lnls_calc_impedance_multilayer_round_pipe(wi, epr, mur, b, L, E);
%     
%     wi  =  [    -wi,  fliplr(wi)];
%     Zli = [ conj(Zli), fliplr(Zli)];
%     Zvi = [-conj(Zvi), fliplr(Zvi)];
%     Zhi = [-conj(Zhi), fliplr(Zhi)];
%     
% %     [~,inds] = sort(wi);
% %     Zvi(inds) = conv(Zvi(inds),fil,'same');
%     
%     [wake(i,:), tau(i,:)] = inverse_fft(wi, -1i/pi*Zvi);
% end
% 
% figure; axes; 
% 
% for i = 1:size(tau,1)
%     loglog(-tau(i,:),-conv(wake(i,:),fil,'same'));
%     hold all;
% end


% %% Broad Band Test
% nfreq = 300001;
% w_max = 4e12;
% w = linspace(0,w_max,nfreq); 
% w =  [    -w,  fliplr(w)];
% 
% 
% 
% Zl = Rs./(1+1i*Q*(wr./w - w/wr));
% Zv = wr*Rs/radius./(w + 1i*Q*(wr - w.^2/wr));
% 
% [wake, tau] = inverse_fft(w, -1i/pi*Zv);
% 
% wakeST = wr*Rs/radius/Ql*sin(wrl*tau).*exp(wr*tau/(2*Q));
% figure; plot(tau,[imag(wake);real(wake);wakeST]);
% figure; plot(w,[real(Zv);imag(Zv)]);

%% Quadrature Integration Test:

% But the method wich have showed the best results was the quadrature
% integral of the fourier transform (at least for the short range wake);
tau  = -[0, logspace(-14,-9,300)];

wake = fourier_transform_by_quadgk(@multilayer_impedance,tau,1e-5,16);
% Lets compare it with Bane's formula:
wakeST = W0*(8/3*exp(tau*c/s0).*sin(-tau*c*sqrt(3)/s0-pi/6) ...
    + 1/sqrt(2*pi)*1./(a^p + (-tau*c/(4*s0)).^(p/2)).^(1/p));

figure; loglog(-tau,[abs(wake);abs(wakeST)]);


% Zv = (@(x) wr*Rs/radius./(x + 1i*Q*(wr - x.^2/wr)));
% wake = fourier_transform_by_quadgk(Zv,tau,0,16);
% 
% wakeST = wr*Rs/radius/Ql*sin(wrl*tau).*exp(wr*tau/(2*Q));
% figure; plot(tau,[wake;wakeST]);
