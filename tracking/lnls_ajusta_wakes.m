function [Rs,wr,Q, wake] = lnls_ajusta_wakes(tau,wake,plane,Rs,wr,Q,coef_ang,coef_lin)

x0 = [Rs/1e14,wr/1e13,Q];

if strcmpi(plane,'long')
    fun = @funl;
else
    fun = @funt;
end
options = optimset('MaxIter',3000);
x = lsqcurvefit(fun,x0,tau,wake,[],[],options);

wake = fun(x,tau);

len = length(x);
Rs = x(1:len/3)*1e14;
wr = x((len/3+1):(2*len/3))*1e13;
Q =  x((2*len/3+1):len);


function wake_test = funt(x,tau)

len = length(x);
Rs = x(1:len/3)*1e14;
wr = x((len/3+1):(2*len/3))*1e13;
Q =  x((2*len/3+1):len);

wake_test = lnls_calc_wake_transverse_resonator(Rs,Q,wr,tau);


function wake_test = funl(x,tau)

len = length(x);
Rs = x(1:len/3)*1e14;
wr = x((len/3+1):(2*len/3))*1e13;
Q =  x((2*len/3+1):len);

wake_test = lnls_calc_wake_longitudinal_resonator(Rs,Q,wr,tau);