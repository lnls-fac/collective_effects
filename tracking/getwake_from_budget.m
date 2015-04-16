function [wakel,wakeh,wakev] = getwake_from_budget(w,budget,tau)


Zl = zeros(size(w));
Zv = zeros(size(w));
Zh = zeros(size(w));
for i=1:length(budget)
    Zl = Zl + budget{i}.quantity * budget{i}.Zl;
    Zv = Zv + budget{i}.quantity * budget{i}.Zv * budget{i}.betay;
    Zh = Zh + budget{i}.quantity * budget{i}.Zh * budget{i}.betax;
end
% total_budget{1}.name = 'All Impedances';
% total_budget{1}.Zl = Zl;
% total_budget{1}.Zv = Zv;
% total_budget{1}.Zh = Zh;
% total_budget{1}.quantity = 1;

% plot_impedances(w,total_budget,false,'log','log',false);

wakel = zeros(size(tau));
wakev = zeros(size(tau));
wakeh = zeros(size(tau));
%% This lines were an attempt to fit the wake with only resonators
% Rs = zeros(1,length(budget));
% Q  = zeros(1,length(budget));
% wr = zeros(1,length(budget));
for i=1:length(budget)
    if ~strcmpi(budget{i}.name,'Broad Band')
        Zh = budget{i}.quantity * budget{i}.Zh * budget{i}.betax;
        
        %% This lines were an attempt to fit the wake with only resonators
%         radius = budget{i}.b(1);
%         sigma = budget{i}.sigmadc(2);
%         Z0 = c*4*pi*1e-7;
%         L = budget{i}.quantity * budget{i}.L * budget{i}.betax;
%         s0 = (2*radius^2/Z0/sigma)^(1/3);
%         W0 = c*Z0/4/pi * 2*s0*L/radius^4;
%         Qx = sqrt(3)/2;
%         wrx = c*sqrt(3)/s0;
%         Rsx = W0*8/3*Qx/wrx;
        
        wakeh_fun = (@(x)interp1(w,Zh,x,'linear',0));
        % wakel_fun = (@(x)wr*Rs/radius./(x + 1i*Q*(wr - x.^2/wr)));
        wake = fourier_transform_by_quadgk(wakeh_fun,tau,'trans',1e-6,20);
        
%         [Rsx,wrx,Qx, wake_fit] = lnls_ajusta_wakes(tau,wake,'trans',Rsx,wrx,Qx);
        
        wakeh = wakeh + wake;
    else
        Rsx = budget{i}.Rsx * budget{i}.quantity * budget{i}.betax;
        wrx = budget{i}.wrx;
        Qx = budget{i}.Qx;
        
        wakeh = wakeh +  lnls_calc_wake_transverse_resonator(Rsx,Qx,wrx,tau);
    end
%     Rs(i) = Rsx; Q(i) = Qx; wr(i) = wrx;
end
% [Rsx,wrx,Qx, wake_fit] = lnls_ajusta_wakes(tau,wakeh,'trans',Rs,wr,Q);

for i=1:length(budget)
    if ~strcmpi(budget{i}.name,'Broad Band')
        Zv = budget{i}.quantity * budget{i}.Zv * budget{i}.betay;
        
        wakev_fun = (@(x)interp1(w,Zv,x,'linear',0));
        % wakel_fun = (@(x)wr*Rs/radius./(x + 1i*Q*(wr - x.^2/wr)));
        wakev = wakev + fourier_transform_by_quadgk(wakev_fun,tau,'trans',1e-6,20);
    else
        Rsy = budget{i}.Rsy;
        wry = budget{i}.wry;
        Qy = budget{i}.Qy;
        
        wakev = wakev + budget{i}.quantity * budget{i}.betay * lnls_calc_wake_transverse_resonator(Rsy,Qy,wry,tau);
    end
end
for i=1:length(budget)
    if ~strcmpi(budget{i}.name,'Broad Band')
        Zl = budget{i}.quantity * budget{i}.Zl;
        
        wakel_fun = (@(x)interp1(w,Zl,x,'linear',0));
        % wakel_fun = (@(x)Rs./(1+1i*Q*(wr./x - x/wr)));
        wakel = wakel + fourier_transform_by_quadgk(wakel_fun,tau,'long',1e-6,20);
    else
        Rsl = budget{i}.Rsl;
        wrl = budget{i}.wrl;
        Ql = budget{i}.Ql;
        
        wakel = wakel + budget{i}.quantity *lnls_calc_wake_longitudinal_resonator(Rsl,Ql,wrl,tau);
    end
end