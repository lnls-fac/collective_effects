function deltaw = lnls_calc_longitudinal_cbi(w, Zl, sigma, nb, w0, nus, eta, E, I_tot, m)
% Calcula a impedancia longitudinal efetiva dos nb modos de oscilacao,
% considerando um feixe gaussiano, para o modo azimutal m e radial k=0;
% E calcula as instabilidades de Coupled_bunch a partir dela.

c = 299792458;

%% Calculate Effective Impedance
pmin = ceil((w(1)-m*nus*w0)/(w0*nb)); % arredonda em direcao a +infinito
pmax = ceil((w(end)-(nb-1+m*nus)*w0)/(w0*nb))-1; % arredonda em direcao a -infinito

p = pmin:pmax;
wp = w0*(bsxfun(@plus, p*nb, (0:nb-1)') + m*nus);

h = (wp*sigma/c).^(2*abs(m)).*exp(-(wp*sigma/c).^2);
interpol_Z = interp1(w,Zl(:),wp);
Zl_eff = diag((interpol_Z./wp)*h').';


deltaw = 1i/(2*pi*2^m*factorial(m - 1)) * I_tot*eta/(E*1e9*nus*(sigma/c)^2) * Zl_eff;