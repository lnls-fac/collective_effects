function [espread, ave_tau,std_tau, tune] = longitudinal_threshold_experiment(ring, bunch, wake, I_b, plota)

espread   = zeros(1,length(I_b));
distr     = zeros(length(I_b),length(bunch.tau));
std_tau   = zeros(1,length(I_b));
ave_tau   = zeros(1,length(I_b));
tune      = zeros(1,length(I_b));
% if numel(bunch.espread) > 1
%     espr = interp1(bunch.I, bunch.espread, bunch.I_b, 'linear', bunch.espread(end));
% else
%     espr = bunch.espread;
% end
for i=1:length(I_b)
    bunch.I_b = I_b(i);
    [~, espread(i), ~, distr(i,:), tune(i)] = generate_longitudinal_bunch(bunch, ring, wake);%, espr);
%     espr = espread(i);
    distr(i,:) = distr(i,:);
    ave_tau(i) = trapz(bunch.tau,bunch.tau.*distr(i,:));
    std_tau(i) = sqrt(trapz(bunch.tau,bunch.tau.^2.*distr(i,:))-ave_tau(i)^2);
    fprintf('.');
end
fprintf('\n');
ind = abs(bunch.tau) < 3*max(std_tau);

if plota
    figure; plot(I_b*1e3,ave_tau*299792458*1e3);
    figure; plot(I_b*1e3,std_tau*299792458*1e3);
    figure; plot(I_b*1e3,tune*1e3);
    figure; plot(I_b*1e3, espread*100);
    
    [I,T] = meshgrid(I_b*1e3,bunch.tau*299792458*1e3);
    figure;  surface(I(ind,:),T(ind,:),distr(:,ind)','LineStyle','none');
end