function [wake, tau] = inverse_fft(w,Z)
w_max = max(w);
n = ceil(length(w)/2);

tau = -(0:n-1)*pi/w_max;
wake = ifft(Z,'symmetric')*w_max;
wake = wake(1:n);