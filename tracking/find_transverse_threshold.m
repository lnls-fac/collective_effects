function [I0, tune] = find_transverse_threshold(ring,bunch,wake, tau, lim_inf, lim_sup)

np = 1000;%[500,5000];
Tol = 0.03e-3;

for i=1:length(np)
%     bunch.num_part = np(i);
    
    I = [lim_inf,lim_sup];
    [~, avex, gtime] = transverse_threshold_experiment(ring, bunch, wake, I, false);
    gtime_inf = gtime(1);
    gtime_sup = gtime(2);
    lim_inf_old = lim_inf/1.5;
    lim_sup_old = lim_sup*1.5;
    gtime_sup_old = 2*1/tau;
    gtime_inf_old = 1/tau/2;
    while (lim_sup-lim_inf) > Tol
        
        if (gtime_inf-1/tau) > 0
            lim_sup_old = lim_sup;
            gtime_sup_old = gtime_sup;
            lim_sup = lim_inf;
            gtime_sup = gtime_inf;
            lim_inf = lim_inf_old;
            gtime_inf = gtime_inf_old;
%             [~, avex, gtime_inf] = transverse_threshold_experiment(ring, bunch, wake, lim_inf, false);
        elseif (gtime_sup-1/tau) > 0
            lim_sup_old = lim_sup;
            gtime_sup_old = gtime_sup;
            lim_sup = (lim_sup+lim_inf)/2;
            [~, avex, gtime_sup] = transverse_threshold_experiment(ring, bunch, wake, lim_sup, false);
        else
            lim_inf_old = lim_inf;
            gtime_inf_old = gtime_inf;
            lim_inf = lim_sup;
            gtime_inf = gtime_sup;
            lim_sup = lim_sup_old;
            gtime_sup = gtime_sup_old;
%             [~, avex, gtime_sup] = transverse_threshold_experiment(ring, bunch, wake, lim_sup, false);
        end

    end
    I0 = (lim_sup+lim_inf) /2;
    lim_sup = I0 + 4*(length(np)-i)*Tol;
    lim_inf = I0 - 4*(length(np)-i)*Tol;
end

n = ring.nturns;
tune = (0:n/2)/n;
fft_ave = 2*abs(fft(avex));
[~,idx] = max(fft_ave(1:n/2+1));
tune = tune(idx);