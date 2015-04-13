function delta = lnls_calc_transverse_mode_couplingopt(w, Z, params)

n_rad = params.n_rad;
n_azi = params.n_azi;
sigma = params.sigma;
I_b   = params.I;
E     = params.E;
w0    = params.w0;
nus   = params.nus;
nut   = params.nut;
chrom = params.chrom;
eta   = params.eta;
nb    = params.nb;
mu    = params.mu;

c = 299792458;

nucro = nut/eta*chrom;
pmin = ceil((w(1)-(n_azi*nus + mu + nut)*w0)/(w0*nb)); % arredonda em direcao a +infinito
pmax = ceil((w(end)-(nb-1 + mu + nut + n_azi*nus)*w0)/(w0*nb))-1; % arredonda em direcao a -infinito

p = pmin:pmax;
wp = w0*(p*nb + mu + nut + 0*nus);
wpcro = wp - nucro*w0;
interpol_Z = interp1(w,Z,wp);



A = zeros(1 + 2*n_azi, 1+n_rad, 1 + 2*n_azi, 1+n_rad);
M = zeros(1 + 2*n_azi, 1+n_rad, 1 + 2*n_azi, 1+n_rad);
delta = zeros((1 + 2*n_azi)*(1+n_rad),length(I_b));

if length(sigma)~=1
    for ii=1:length(I_b)
        K = I_b(ii)*nb*w0/(4*pi)/(nus*w0)/(E*1e9);
        for k = 0:n_rad
            for m = 0:n_azi
                Imk = 1/sqrt(factorial(abs(m)+k)*factorial(k))* ...
                    (wpcro*sigma(ii)/c/sqrt(2)).^(abs(m)+2*k).*exp(-(wpcro*sigma(ii)/c/sqrt(2)).^2);
                A(n_azi+1 + m, 1+k, n_azi+1 + m, 1+k) = m;
                A(n_azi+1 - m, 1+k, n_azi+1 - m, 1+k) = -m;
                for n = m:n_azi
                    Inl = 1/sqrt(factorial(abs(n)+k)*factorial(k)) .* ...
                        (wpcro*sigma(ii)/c/sqrt(2)).^(abs(n)+2*k).*exp(-(wpcro*sigma(ii)/c/sqrt(2)).^2);
                    M(n_azi+1 + m, 1+k, n_azi+1 + n, 1+k) = -1i*(1i)^(m-n)*K*sum(interpol_Z.*Imk.*Inl);
                    M(n_azi+1 - m, 1+k, n_azi+1 + n, 1+k) = M(n_azi+1 + m, 1+k, n_azi+1 + n, 1+k);
                    M(n_azi+1 + m, 1+k, n_azi+1 - n, 1+k) = M(n_azi+1 + m, 1+k, n_azi+1 + n, 1+k);
                    M(n_azi+1 - m, 1+k, n_azi+1 - n, 1+k) = M(n_azi+1 + m, 1+k, n_azi+1 + n, 1+k);
                    M(n_azi+1 + n, 1+k, n_azi+1 + m, 1+k) = (-1)^(m-n)*M(n_azi+1 + m, 1+k, n_azi+1 + n, 1+k);
                    M(n_azi+1 + n, 1+k, n_azi+1 - m, 1+k) = M(n_azi+1 + n, 1+k, n_azi+1 + m, 1+k);
                    M(n_azi+1 - n, 1+k, n_azi+1 + m, 1+k) = M(n_azi+1 + n, 1+k, n_azi+1 + m, 1+k);
                    M(n_azi+1 - n, 1+k, n_azi+1 - m, 1+k) = M(n_azi+1 + n, 1+k, n_azi+1 + m, 1+k);
                end
            end
            for l = (k+1):n_rad
                for m = 0:n_azi
                    Imk = 1/sqrt(factorial(abs(m)+k)*factorial(k))* ...
                        (wpcro*sigma(ii)/c/sqrt(2)).^(abs(m)+2*k).*exp(-(wpcro*sigma(ii)/c/sqrt(2)).^2);
                    for n = 0:n_azi
                        Inl = 1/sqrt(factorial(abs(n)+l)*factorial(l)) .* ...
                            (wpcro*sigma(ii)/c/sqrt(2)).^(abs(n)+2*l).*exp(-(wpcro*sigma(ii)/c/sqrt(2)).^2);
                        M(n_azi+1 + m, 1+k, n_azi+1 + n, 1+l) = -1i*(1i)^(m-n)*K*sum(interpol_Z.*Imk.*Inl);
                        M(n_azi+1 - m, 1+k, n_azi+1 + n, 1+l) = M(n_azi+1 + m, 1+k, n_azi+1 + n, 1+l);
                        M(n_azi+1 + m, 1+k, n_azi+1 - n, 1+l) = M(n_azi+1 + m, 1+k, n_azi+1 + n, 1+l);
                        M(n_azi+1 - m, 1+k, n_azi+1 - n, 1+l) = M(n_azi+1 + m, 1+k, n_azi+1 + n, 1+l);
                        M(n_azi+1 + n, 1+l, n_azi+1 + m, 1+k) = (-1)^(m-n)*M(n_azi+1 + m, 1+k, n_azi+1 + n, 1+l);
                        M(n_azi+1 + n, 1+l, n_azi+1 - m, 1+k) = M(n_azi+1 + n, 1+l, n_azi+1 + m, 1+k);
                        M(n_azi+1 - n, 1+l, n_azi+1 + m, 1+k) = M(n_azi+1 + n, 1+l, n_azi+1 + m, 1+k);
                        M(n_azi+1 - n, 1+l, n_azi+1 - m, 1+k) = M(n_azi+1 + n, 1+l, n_azi+1 + m, 1+k);
                    end
                end
            end
        end
        
        B = A + M;
        B = reshape(B,(1 + 2*n_azi)*(1+n_rad), (1 + 2*n_azi)*(1+n_rad));
        delta(:,ii) = eig(B);
    end
    
else
    for k = 0:n_rad
        for m = 0:n_azi
            Imk = 1/sqrt(factorial(abs(m)+k)*factorial(k))* ...
                (wpcro*sigma/c/sqrt(2)).^(abs(m)+2*k).*exp(-(wpcro*sigma/c/sqrt(2)).^2);
            A(n_azi+1 + m, 1+k, n_azi+1 + m, 1+k) = m;
            A(n_azi+1 - m, 1+k, n_azi+1 - m, 1+k) = -m;
            for n = m:n_azi
                Inl = 1/sqrt(factorial(abs(n)+k)*factorial(k)) .* ...
                    (wpcro*sigma/c/sqrt(2)).^(abs(n)+2*k).*exp(-(wpcro*sigma/c/sqrt(2)).^2);
                M(n_azi+1 + m, 1+k, n_azi+1 + n, 1+k) = -1i*(1i)^(m-n)*sum(interpol_Z.*Imk.*Inl);
                M(n_azi+1 - m, 1+k, n_azi+1 + n, 1+k) = M(n_azi+1 + m, 1+k, n_azi+1 + n, 1+k);
                M(n_azi+1 + m, 1+k, n_azi+1 - n, 1+k) = M(n_azi+1 + m, 1+k, n_azi+1 + n, 1+k);
                M(n_azi+1 - m, 1+k, n_azi+1 - n, 1+k) = M(n_azi+1 + m, 1+k, n_azi+1 + n, 1+k);
                M(n_azi+1 + n, 1+k, n_azi+1 + m, 1+k) = (-1)^(m-n)*M(n_azi+1 + m, 1+k, n_azi+1 + n, 1+k);
                M(n_azi+1 + n, 1+k, n_azi+1 - m, 1+k) = M(n_azi+1 + n, 1+k, n_azi+1 + m, 1+k);
                M(n_azi+1 - n, 1+k, n_azi+1 + m, 1+k) = M(n_azi+1 + n, 1+k, n_azi+1 + m, 1+k);
                M(n_azi+1 - n, 1+k, n_azi+1 - m, 1+k) = M(n_azi+1 + n, 1+k, n_azi+1 + m, 1+k);
            end
        end
        for l = (k+1):n_rad
            for m = 0:n_azi
                Imk = 1/sqrt(factorial(abs(m)+k)*factorial(k))* ...
                    (wpcro*sigma/c/sqrt(2)).^(abs(m)+2*k).*exp(-(wpcro*sigma/c/sqrt(2)).^2);
                for n = 0:n_azi
                    Inl = 1/sqrt(factorial(abs(n)+l)*factorial(l)) .* ...
                        (wpcro*sigma/c/sqrt(2)).^(abs(n)+2*l).*exp(-(wpcro*sigma/c/sqrt(2)).^2);
                    M(n_azi+1 + m, 1+k, n_azi+1 + n, 1+l) = -1i*(1i)^(m-n)*sum(interpol_Z.*Imk.*Inl);
                    M(n_azi+1 - m, 1+k, n_azi+1 + n, 1+l) = M(n_azi+1 + m, 1+k, n_azi+1 + n, 1+l);
                    M(n_azi+1 + m, 1+k, n_azi+1 - n, 1+l) = M(n_azi+1 + m, 1+k, n_azi+1 + n, 1+l);
                    M(n_azi+1 - m, 1+k, n_azi+1 - n, 1+l) = M(n_azi+1 + m, 1+k, n_azi+1 + n, 1+l);
                    M(n_azi+1 + n, 1+l, n_azi+1 + m, 1+k) = (-1)^(m-n)*M(n_azi+1 + m, 1+k, n_azi+1 + n, 1+l);
                    M(n_azi+1 + n, 1+l, n_azi+1 - m, 1+k) = M(n_azi+1 + n, 1+l, n_azi+1 + m, 1+k);
                    M(n_azi+1 - n, 1+l, n_azi+1 + m, 1+k) = M(n_azi+1 + n, 1+l, n_azi+1 + m, 1+k);
                    M(n_azi+1 - n, 1+l, n_azi+1 - m, 1+k) = M(n_azi+1 + n, 1+l, n_azi+1 + m, 1+k);
                end
            end
        end
    end
    for ii=1:length(I_b)
        K = I_b(ii)*nb*w0/(4*pi)/(nus*w0)/(E*1e9);
        B = A + K*M;
        B = reshape(B,(1 + 2*n_azi)*(1+n_rad), (1 + 2*n_azi)*(1+n_rad));
        delta(:,ii) = eig(B);
    end
end