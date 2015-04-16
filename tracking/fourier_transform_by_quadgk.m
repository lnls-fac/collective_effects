function wake = fourier_transform_by_quadgk(impedance_func,tau,plane, multi_min, multi_max)

wmax =  abs(multi_max*2*pi./tau);
wmin =  abs(multi_min*2*pi./tau);

ind0 = find(tau == 0);
wmax(ind0) = wmax(ind0+1);
wmin(ind0) = wmin(ind0+1);

wake = zeros(1,length(tau));
RelTol = 1e-4;
for i = 1:length(tau)
    if strcmpi(plane,'long')
        fun_wake = (@(x)1/2/pi*(real(impedance_func(x)) .* cos(x*tau(i)) - ...
                                imag(impedance_func(x)) .* sin(x*tau(i))));
    else
        fun_wake = (@(x)1/2/pi*(real(impedance_func(x)) .* sin(x*tau(i)) + ...
                                imag(impedance_func(x)) .* cos(x*tau(i))));
    end
    [wake(i), err] = quadgk(fun_wake,wmin(i),wmax(i),'AbsTol',0,'RelTol',RelTol);%,'MaxIntervalCount',2000);
    if err/wake(i) >RelTol, disp(tau(i)*1e12); error('integral did not converge.');end
    fprintf('.');
    if ~mod(i,80), fprintf('\n');end
end
fprintf('\n');

wake = 2*wake;
if ~strcmpi(plane,'long')
    wake(ind0) = 0;
end