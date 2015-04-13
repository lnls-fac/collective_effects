function [ave_bun, rms_bun] = tracking_test()
el_ch = 1.602e-19;
c = 299792458;

% ring data
param.nus = 3.93e-3;
nus = param.nus;
param.nut0 = 14.647;
param.chrom = 0.0;
param.T0  = 518.39/299792458;
param.E   = 3e9;
param.alpha = 1.7e-4;
param.beta = 11;
param.nb = 2;
nb = param.nb;
param.np = [50 1]; %repmat(1,1,nb);
np = param.np;
param.Ne = [2.5 0.01]*1e-3/el_ch*param.T0; %repmat(2.5e-3/el_ch*param.T0,1,nb);
param.nr = 10000;
nr = param.nr;
param.ntrack = 5000;


part = zeros(4,sum(np));
ave_bun = zeros(4,nb,nr);
rms_bun = zeros(4,nb,nr);

% % definicao da impedancia
% Zovern = 0.2;
% fr  = 2.4* 299792458/12e-3/2/pi; % 2.4 c/b/2/pi;
% beta_imp = 11;
% Rs = Zovern*fr/0.578e6/12e-3;
% Q = 1;
% wr = fr*2*pi;
% Ql = sqrt(Q^2 - 1/4); wrl = wr*Ql/Q;
% tau = linspace(1e-8,1,10000)/c;
% wakefun = beta_imp*wr*Rs/Ql*sin(wrl*tau).*exp(-wr*tau/2/Q);
% 
% b = 5;
% sigma = 6e7;
% Z0 = c*4*pi*1e-7;
% a = 3/sqrt(2*pi)/4;
% p = 2.7;
% s0 = (2*b^2/Z0/sigma)^(1/3);
% L = 480e9;
% W0 = c*Z0/4/pi * 2*s0*L/b^4;
% % tau = linspace(1e-8,1,10000)/c;
% wakeres= W0*(8/3*exp(-tau*c/s0).*sin(tau*c*sqrt(3)/s0-pi/6) ...
%     + 1/sqrt(2*pi)*1./(a^p + (tau*c/(4*s0)).^p).^(1/p));
% figure; plot(tau*c, [wakeres; wakefun]);

RandStream.setGlobalStream(RandStream('mt19937ar','seed', 190488));
%definicao dos pacotes
for ii = 1:nb
    ini = 1 + sum(np(1:ii-1));
    fim = sum(np(1:ii));
    vec = ini:fim;
    part(1,vec) = randn(1,np(ii))*7e-4;
    part(2,vec) = randn(1,np(ii))*3.8e-3/299792458;
    part(3,vec) = randn(1,np(ii))*1e-9;
    part(4,vec) = randn(1,np(ii))*1e-9;
    if np(ii) >= 1
        ave_bun(:,ii,1) = mean(part(:,vec),2);
        rms_bun(:,ii,1) = rms(part(:,vec),2);
    end
end
part_old = part;



for ii=2:nr;
    part(1,:) = part_old(1,:) + (2*pi*nus).^2/param.alpha/param.T0.*part_old(2,:);
    part(2,:) = part_old(2,:) - param.alpha*param.T0.*part(1,:);
    nut  =  param.nut0*(1 + param.chrom*part_old(1,:));
    part(3,:) = part_old(3,:).*cos(2*pi*nut) + param.beta*part_old(4,:).*sin(2*pi*nut);
    part(4,:) = -1/param.beta*part_old(3,:).*sin(2*pi*nut) + part_old(4,:).*cos(2*pi*nut);
    if mod(ii,100)==0
        fprintf('%d\n',ii);
    end
    for jj=1:param.nb
        ini = 1 + sum(param.np(1:jj-1));
        fim = sum(param.np(1:jj));
        vec = ini:fim;
        if param.np(jj) >= 1
            ave_bun(:,jj,ii) = mean(part(:,vec),2);
            rms_bun(:,jj,ii) = rms(part(:,vec),2);
        end
    end
    kick = kick_impedance(part([2 3],:), ave_bun([2 3],:,:), ii, param, 0, 0);
    part(4,:) = part(4,:) + kick;
    part_old = part;
end

figure; plot(squeeze(ave_bun(3,1,:)));



function kick = kick_impedance(part, ave_bun, turn, param, wake, tau)
el_ch = 1.602e-19;
c = 299792458;
Z0 = c*4*pi*1e-7;

kick = zeros(size(part(1,:)));

Y = el_ch*param.Ne./param.E;

% definicao da impedancia: está mais rápido que a interpolacao
Zovern = 0.2;
fr  = 2.4* 299792458/12e-3/2/pi; % 2.4 c/b/2/pi;
beta_imp = 11;
Rs = Zovern*fr/0.578e6/12e-3;
Q = 1;
wr = fr*2*pi;
Ql = sqrt(Q^2 - 1/4); wrl = wr*Ql/Q;

b = 12e-3;
sigma = 5.9e7;
a = 3/sqrt(2*pi)/4;
p = 2.7;
s0 = (2*b^2/Z0/sigma)^(1/3);
L = 4800;
W0 = c*Z0/4/pi * 2*s0*L/b^4;

for ii=1:param.nb
    ini = 1 + sum(param.np(1:ii-1));
    fim = sum(param.np(1:ii));
    vec = ini:fim;
    if param.np(ii) > 1
        difft = bsxfun(@minus,part(1,vec)',part(1,vec));
%         kik = beta_imp*wr*Rs/Ql*sin(wrl*difft).*exp(-wr*difft/2/Q);
        kik = beta_imp*W0*(8/3*exp(-difft*c/s0).*sin(difft*c*sqrt(3)/s0-pi/6) ...
            + 1/sqrt(2*pi)*1./(a^p + (difft*c/(4*s0)).^p).^(1/p));
        ind = difft <= 0;
        kik(ind) = 0;
        desl = part(2,vec); desl = desl(:)';
        kick(1,vec) = Y(ii).*(desl*kik)/param.beta/param.np(ii);
    end
    
    if turn <= param.ntrack
        track = turn;
    else
        track = param.ntrack;
    end
    time_turns = param.T0*(0:track-1);
    time_bun = param.T0/param.nb*((-ii+1):(param.nb-ii));
    difft = bsxfun(@plus,time_turns',time_bun);
    difft = difft + squeeze(ave_bun(1,:,turn:-1:(turn-track+1)))' - ave_bun(1,ii,turn);
    kik = beta_imp*W0*( ... %8/3*exp(-difft*c/s0).*sin(difft*c*sqrt(3)/s0-pi/6) ...
        + 1/sqrt(2*pi)*1./(a^p + (difft*c/(4*s0)).^p).^(1/p));
    ind = difft <= 0;
    kik(ind) = 0;
    kick(1,vec) =  kick(1,vec) + sum(Y*(squeeze(ave_bun(2,:,turn:-1:(turn-track+1))).*kik')/param.beta);
end




% for ii=1:param.nb
%     difft = bsxfun(@minus,squeeze(part(1,ii,:))',squeeze(part(1,ii,:)));
%     kik = interp1(tau,wake,difft,'linear',0);
%     desl = squeeze(part(2,ii,:)); desl = desl(:)';
%     kick(1,ii,:) = Y(ii).*(desl*kik)/param.beta;
% end
