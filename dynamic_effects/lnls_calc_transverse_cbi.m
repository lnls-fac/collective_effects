function deltaw = lnls_calc_transverse_cbi(w,Z, sigma, nb, w0, nus, nut, chrom, eta, m, E, I_tot)
% Calcula a impedância transversal efetiva dos nb modos de oscilacao,
% considerando um feixe gaussiano, para o modo azimutal m e radial k=0;
% E calcula as instabilidades de Coupled_bunch a partir dela.
%
% deltaw = lnls_calc_transverse_cbi(w,Z, sigma, nb, w0, nus, nut, chrom, eta, m, E, I_tot)

c = 299792458;

%% Calculate Effective Impedance
nucro = chrom/eta;
pmin = ceil((w(1)-(m*nus + nut)*w0)/(w0*nb)); % arredonda em direcao a +infinito
pmax = ceil((w(end)-(nb-1 + nut + m*nus)*w0)/(w0*nb))-1; % arredonda em direcao a -infinito

p = pmin:pmax;
wp = w0*(bsxfun(@plus, p*nb, (0:nb-1)') + nut + m*nus);

wpcro = wp - nucro*w0;
h = (wpcro*sigma/c).^(2*abs(m)).*exp(-(wpcro*sigma/c).^2);
interpol_Z = interp1(w,Z(:),wp);
Zt_eff = diag(interpol_Z*h').';

%% Calculate Coupled_bunch Instabilities

deltaw = -1i*I_tot/(4*pi*2^abs(m)*factorial(abs(m))*E*1e9) *w0* (Zt_eff);
