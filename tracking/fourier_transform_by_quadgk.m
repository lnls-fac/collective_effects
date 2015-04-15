function wake = fourier_transform_by_quadgk(impedance_func,tau, multi_min, multi_max)

wmax =  abs(multi_max*2*pi./tau);
wmin =  abs(multi_min*2*pi./tau);

ind0 = find(tau == 0);
wmax(ind0) = 1e16;
wmin(ind0) = 1e9;

wake = zeros(1,length(tau));
for i = 1:length(tau)
    fun_wake = (@(x)-1i/2/pi*impedance_func(x) .* exp(1i*x*tau(i)));
    [wake(i), err] = quadgk(fun_wake,wmin(i),wmax(i),'AbsTol',0,'RelTol',1e-4);%,'MaxIntervalCount',2000);
    if err/wake(i) >1e-4, error('integral did not converge.');end
    fprintf('.');
    if ~mod(i,80), fprintf('\n');end
end
fprintf('\n');

wake = 2*real(wake);
wake(ind0) = 0;