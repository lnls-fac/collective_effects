function [rmsx, avex, gtime] = transverse_threshold_experiment(ring, bunch, wake, I_b, plota)

fft_ave = zeros(length(I_b),ring.nturns);
rmsx    = zeros(length(I_b),ring.nturns);
avex    = zeros(length(I_b),ring.nturns);
gtime   = zeros(1,length(I_b));
n = ring.nturns;

for i=1:length(I_b)
    bunch.I_b = I_b(i);
    [ave_bun,rms_bun, ~, ~, ~] = single_bunch_tracking(ring, bunch, wake);
    fft_ave(i,:) = 2*abs(fft(ave_bun(1,:)));
    rmsx(i,:) = rms_bun(1,:);
    avex(i,:) = ave_bun(1,:);
    growth = lnls_polyfit((n/2:n)*ring.rev_time,log(rmsx(i,(n/2:n))),[0,1]);
    gtime(i) = growth(2);
    fprintf('%d : %5.3f mA,  %9.5g Hz\n',i,I_b(i)*1e3, gtime(i));
end

if plota
    nom_tune = rem(ring.tune,1);
    tune = ((0:n/2)/n - nom_tune)*1e3;
    ind = tune > -10 & tune < 10;
    [I,T] = meshgrid(I_b*1e3,tune);
    
    
    pfft = fft_ave(:,1:n/2+1)';
    pfft = pfft./repmat(max(pfft),n/2+1,1);
    figure;  surface(I(ind,:),T(ind,:),pfft(ind,:),'LineStyle','none');
    xlim([min(I_b),max(I_b)]*1e3);ylim([min(tune(ind)),max(tune(ind))]);
    
    figure; plot(I_b*1e3,gtime);
end
