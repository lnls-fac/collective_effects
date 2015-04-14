n1 = 8000;
n2 = 8000;

mu0 = 4*pi*1e-7;
c   = 299792458;
ep0 = 1/c^2/mu0;
E = 3;

epb     = [1 1 1 1];
mub     = [1 1 1 1];
ange    = [0 0 0 0];
angm    = [0 0 0 0];
sigmadc = [0 4e6 5.9e7 1]; % Not sure about NEG resistivity
tau     = [0 0 1 0]*27e-15;
b       = [12.000 12.001 13.000]*1e-3;
L       = 480;
Zovern = 0.2;
radius = 12e-3;
fr  = 2.4* c/(radius*2*pi)/10; 
Rs = Zovern*fr*(518.396/c);
Q = 1;
wr = fr*2*pi;
Ql = sqrt(Q.^2 - 1/4);
wrl = wr .* Ql ./ Q;

% w1  = linspace(1e4,1e7,n1);
% w1  = [-fliplr(w1) w1];
% 
% % epr = zeros(length(epb),length(w1));
% % mur = zeros(length(epb),length(w1));
% % for j = 1: length(epb)
% %     epr(j,:) = epb(j)*(1-1i.*sign(w1).*tan(ange(j))) + sigmadc(j)./(1+1i*w1*tau(j))./(1i*w1*ep0);
% %     mur(j,:) = mub(j)*(1-1i.*sign(w1).*tan(angm(j)));
% % end
% % [Zl1, Zv1, Zh1] = lnls_calc_impedance_multilayer_round_pipe(w1, epr, mur, b, L, E);
% 
% Zl1 = Rs./(1+1i*Q*(wr./w1 - w1/wr));
% Zv1 = (c./w1) .* Rs./(1+1i*Q*(wr./w1 - w1/wr));
% 
% 
% tau1 = -((1:n1)-1)*2*pi/w1(end);
% wake1 = -1i*ifft(Zv1);
% % figure; plot(tau1,real(wake1(1:end/2))');
% 
% 
% 
% % w2 = linspace(1e7,1e10,n2); 
% % 
% % epr = zeros(length(epb),length(w2));
% % mur = zeros(length(epb),length(w2));
% % for j = 1: length(epb)
% %     epr(j,:) = epb(j)*(1-1i.*sign(w2).*tan(ange(j))) + sigmadc(j)./(1+1i*w2*tau(j))./(1i*w2*ep0);
% %     mur(j,:) = mub(j)*(1-1i.*sign(w2).*tan(angm(j)));
% % end
% % [Zl2, Zv2, Zh2] = lnls_calc_impedance_multilayer_round_pipe(w2, epr, mur, b, L, E);
% 
% % Zl2 = Rs./(1+1i*Q*(wr./w2 - w2/wr));
% % Zv2 = (c./w2) .* Rs./(1+1i*Q*(wr./w2 - w2/wr));
% 
% % tau2 = -((1:n2)-1)*2*pi/w2(end);
% % wake2 = -1i*ifft(Zv2);
% % % figure; plot(tau2,real(wake2(1:end/2))');
% % figure; plot(1:length(wake3),[imag(wake3);real(wake3)]);

n3 = 200001;
w_max = 1e14;
w3 = linspace(0,w_max,n3); 
w3 = [-w3(2:end), fliplr(w3(2:end))]; %;

epr = zeros(length(epb),length(w3));
mur = zeros(length(epb),length(w3));
for j = 1: length(epb)
    epr(j,:) = epb(j)*(1-1i.*sign(w3).*tan(ange(j))) + sigmadc(j)./(1+1i*w3*tau(j))./(1i*w3*ep0);
    mur(j,:) = mub(j)*(1-1i.*sign(w3).*tan(angm(j)));
end
[Zl3, Zv3, Zh3] = lnls_calc_impedance_multilayer_round_pipe(w3, epr, mur, b, L, E);
% Zv3 = conv(Zv3,ones(1,100)/100,'same');


% Zl3 = Rs./(1+1i*Q*(wr./w3 - w3/wr));
% Zv3 = wr*Rs/radius./(w3 + 1i*Q*(wr - w3.^2/wr));

tau3 = -(0:n3-1)*pi/w_max;
wake3 = ifft(-1i/pi*Zv3,'symmetric')*w_max;
wake3 = wake3(1:n3);


% wakeST = wr*Rs/radius/Ql*sin(wrl*tau3).*exp(wr*tau3/(2*Q));
radius = 12e-3;
sigma = 5.9e7;
Z0 = c*4*pi*1e-7;
a = 3/sqrt(2*pi)/4;
p = 2.7;
L = 480;
% tau = -(0:1000)*1e-12;
s0 = (2*radius^2/Z0/sigma)^(1/3);
W0 = c*Z0/4/pi * 2*s0*L/radius^4;
wakeST = beta_imp*W0*(8/3*exp(tau3*c/s0).*sin(-tau3*c*sqrt(3)/s0-pi/6) ...
    + 1/sqrt(2*pi)*1./(a^p + (-tau3*c/(4*s0)).^(p/2)).^(1/p));


figure; plot(tau3,[imag(wake3);real(wake3);wakeST]);
figure; plot(w3,[real(Zv3);imag(Zv3)]);

% tau = (-1000:0)*1e-12;
% FT = exp(-1i * w' * tau);
% 
% Wake = Zv * FT
% figure; plot(tau,real(Wake));