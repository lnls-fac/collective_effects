

freq = 0.11;
phase = 3/4*pi;

len = 8;
samp = 1:len;

fil = cos(2*pi*freq*samp + phase).*sin(2*pi*samp/(2*len))./(samp*pi);

n = (0:199)/400;
x = cos(2*pi*n'*samp);

figure; plot(n,20*log(abs(fil*x')))

% phi = (0:199)/200 * 2 * pi;
% y = cos(bsxfun(@plus,2*pi*0.12*samp',phi));
% 
% figure; plot(y(end,:), fil*y);