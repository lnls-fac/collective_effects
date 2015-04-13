function delta = lnls_calc_longitudinal_mode_couplingopt(w,Z, params)


n_rad = params.n_rad;
n_azi = params.n_azi;
sigma = params.sigma;
I_b   = params.I;
E     = params.E;
w0    = params.w0;
nus   = params.nus;
eta   = params.eta;
nb    = params.nb;
mu    = params.mu;

c = 299792458;

pmin = ceil((w(1)-(mu + n_azi*nus)*w0)/(w0*nb)); % arredonda em direcao a +infinito
pmax = ceil((w(end)-(mu + n_azi*nus)*w0)/(w0*nb))-1; % arredonda em direcao a -infinito

p = pmin:pmax;
wp = w0*(p*nb + mu + 1*nus);
interpol_Z = interp1(w,Z,wp);


A = zeros(1 + 2*n_azi + n_rad*(2*n_azi+1), 1 + 2*n_azi + n_rad*(2*n_azi+1));
M = zeros(1 + 2*n_azi + n_rad*(2*n_azi+1), 1 + 2*n_azi + n_rad*(2*n_azi+1));
delta = zeros(1 + 2*n_azi + n_rad*(2*n_azi+1),length(I_b));


if length(sigma)~=1
    for ii=1:length(I_b)
        K = I_b(ii)*nb*w0*eta/(2*pi)/(nus*w0)^2/(E*1e9)/(sigma(ii)/c)^2;
        for k = 0:n_rad
            for m = 0:n_azi
                Imk = 1/sqrt(factorial(abs(m)+k)*factorial(k))* ...
                    (wp*sigma(ii)/c/sqrt(2)).^(abs(m)+2*k).*exp(-(wp*sigma(ii)/c/sqrt(2)).^2);
                A(n_azi+1 + m + k*(2*n_azi+1), n_azi+1 + m + k*(2*n_azi+1)) = m;
                A(n_azi+1 - m + k*(2*n_azi+1), n_azi+1 - m + k*(2*n_azi+1)) = -m;
                for n = m:n_azi
                    Inl = 1/sqrt(factorial(abs(n)+k)*factorial(k)) .* ...
                        (wp*sigma(ii)/c/sqrt(2)).^(abs(n)+2*k).*exp(-(wp*sigma(ii)/c/sqrt(2)).^2);
                    Mmknl = 1i*(1i)^(m-n)*K*sum((interpol_Z./wp).*Imk.*Inl);
                    M(n_azi+1 + m + k*(2*n_azi+1), n_azi+1 + n + k*(2*n_azi+1)) = m*Mmknl;
                    M(n_azi+1 - m + k*(2*n_azi+1), n_azi+1 + n + k*(2*n_azi+1)) = - m*Mmknl;
                    M(n_azi+1 + m + k*(2*n_azi+1), n_azi+1 - n + k*(2*n_azi+1)) = m*Mmknl;
                    M(n_azi+1 - m + k*(2*n_azi+1), n_azi+1 - n + k*(2*n_azi+1)) = - m*Mmknl;
                    M(n_azi+1 + n + k*(2*n_azi+1), n_azi+1 + m + k*(2*n_azi+1)) = (-1)^(m-n)*n*Mmknl;
                    M(n_azi+1 + n + k*(2*n_azi+1), n_azi+1 - m + k*(2*n_azi+1)) = (-1)^(m-n)*n*Mmknl;
                    M(n_azi+1 - n + k*(2*n_azi+1), n_azi+1 + m + k*(2*n_azi+1)) = - (-1)^(m-n)*n*Mmknl;
                    M(n_azi+1 - n + k*(2*n_azi+1), n_azi+1 - m + k*(2*n_azi+1)) = - (-1)^(m-n)*n*Mmknl;
                end
            end
            for l = (k+1):n_rad
                for m = 0:n_azi
                    Imk = 1/sqrt(factorial(abs(m)+k)*factorial(k))* ...
                        (wp*sigma(ii)/c/sqrt(2)).^(abs(m)+2*k).*exp(-(wp*sigma(ii)/c/sqrt(2)).^2);
                    for n = 0:n_azi
                        Inl = 1/sqrt(factorial(abs(n)+l)*factorial(l)) .* ...
                            (wp*sigma(ii)/c/sqrt(2)).^(abs(n)+2*l).*exp(-(wp*sigma(ii)/c/sqrt(2)).^2);
                        Mmknl = 1i*(1i)^(m-n)*K*sum((interpol_Z./wp).*Imk.*Inl);
                        M(n_azi+1 + m + k*(2*n_azi+1), n_azi+1 + n + l*(2*n_azi+1)) = m*Mmknl;
                        M(n_azi+1 - m + k*(2*n_azi+1), n_azi+1 + n + l*(2*n_azi+1)) = - m*Mmknl;
                        M(n_azi+1 + m + k*(2*n_azi+1), n_azi+1 - n + l*(2*n_azi+1)) = m*Mmknl;
                        M(n_azi+1 - m + k*(2*n_azi+1), n_azi+1 - n + l*(2*n_azi+1)) = - m*Mmknl;
                        M(n_azi+1 + n + l*(2*n_azi+1), n_azi+1 + m + k*(2*n_azi+1)) = (-1)^(m-n)*n*Mmknl;
                        M(n_azi+1 + n + l*(2*n_azi+1), n_azi+1 - m + k*(2*n_azi+1)) = (-1)^(m-n)*n*Mmknl;
                        M(n_azi+1 - n + l*(2*n_azi+1), n_azi+1 + m + k*(2*n_azi+1)) = - (-1)^(m-n)*n*Mmknl;
                        M(n_azi+1 - n + l*(2*n_azi+1), n_azi+1 - m + k*(2*n_azi+1)) = - (-1)^(m-n)*n*Mmknl;
                    end
                end
            end
        end
        B = A + M;
        delta(:,ii) = eig(B);
    end
else
    for k = 0:n_rad
        for m = 0:n_azi
            Imk = 1/sqrt(factorial(abs(m)+k)*factorial(k))* ...
                (wp*sigma/c/sqrt(2)).^(abs(m)+2*k).*exp(-(wp*sigma/c/sqrt(2)).^2);
            A(n_azi+1 + m + k*(2*n_azi+1), n_azi+1 + m + k*(2*n_azi+1)) = m;
            A(n_azi+1 - m + k*(2*n_azi+1), n_azi+1 - m + k*(2*n_azi+1)) = -m;
            for n = m:n_azi
                Inl = 1/sqrt(factorial(abs(n)+k)*factorial(k)) .* ...
                    (wp*sigma/c/sqrt(2)).^(abs(n)+2*k).*exp(-(wp*sigma/c/sqrt(2)).^2);
                Mmknl = 1i*(1i)^(m-n)*sum((interpol_Z./wp).*Imk.*Inl);
                M(n_azi+1 + m + k*(2*n_azi+1), n_azi+1 + n + k*(2*n_azi+1)) = m*Mmknl;
                M(n_azi+1 - m + k*(2*n_azi+1), n_azi+1 + n + k*(2*n_azi+1)) = - m*Mmknl;
                M(n_azi+1 + m + k*(2*n_azi+1), n_azi+1 - n + k*(2*n_azi+1)) = m*Mmknl;
                M(n_azi+1 - m + k*(2*n_azi+1), n_azi+1 - n + k*(2*n_azi+1)) = - m*Mmknl;
                M(n_azi+1 + n + k*(2*n_azi+1), n_azi+1 + m + k*(2*n_azi+1)) = (-1)^(m-n)*n*Mmknl;
                M(n_azi+1 + n + k*(2*n_azi+1), n_azi+1 - m + k*(2*n_azi+1)) = (-1)^(m-n)*n*Mmknl;
                M(n_azi+1 - n + k*(2*n_azi+1), n_azi+1 + m + k*(2*n_azi+1)) = - (-1)^(m-n)*n*Mmknl;
                M(n_azi+1 - n + k*(2*n_azi+1), n_azi+1 - m + k*(2*n_azi+1)) = - (-1)^(m-n)*n*Mmknl;
            end
        end
        for l = (k+1):n_rad
            for m = 0:n_azi
                Imk = 1/sqrt(factorial(abs(m)+k)*factorial(k))* ...
                    (wp*sigma/c/sqrt(2)).^(abs(m)+2*k).*exp(-(wp*sigma/c/sqrt(2)).^2);
                for n = 0:n_azi
                    Inl = 1/sqrt(factorial(abs(n)+l)*factorial(l)) .* ...
                        (wp*sigma/c/sqrt(2)).^(abs(n)+2*l).*exp(-(wp*sigma/c/sqrt(2)).^2);
                    Mmknl = 1i*(1i)^(m-n)*sum((interpol_Z./wp).*Imk.*Inl);
                    M(n_azi+1 + m + k*(2*n_azi+1), n_azi+1 + n + l*(2*n_azi+1)) = m*Mmknl;
                    M(n_azi+1 - m + k*(2*n_azi+1), n_azi+1 + n + l*(2*n_azi+1)) = - m*Mmknl;
                    M(n_azi+1 + m + k*(2*n_azi+1), n_azi+1 - n + l*(2*n_azi+1)) = m*Mmknl;
                    M(n_azi+1 - m + k*(2*n_azi+1), n_azi+1 - n + l*(2*n_azi+1)) = - m*Mmknl;
                    M(n_azi+1 + n + l*(2*n_azi+1), n_azi+1 + m + k*(2*n_azi+1)) = (-1)^(m-n)*n*Mmknl;
                    M(n_azi+1 + n + l*(2*n_azi+1), n_azi+1 - m + k*(2*n_azi+1)) = (-1)^(m-n)*n*Mmknl;
                    M(n_azi+1 - n + l*(2*n_azi+1), n_azi+1 + m + k*(2*n_azi+1)) = - (-1)^(m-n)*n*Mmknl;
                    M(n_azi+1 - n + l*(2*n_azi+1), n_azi+1 - m + k*(2*n_azi+1)) = - (-1)^(m-n)*n*Mmknl;
                end
            end
        end
    end
    for ii=1:length(I_b)
        K = I_b(ii)*nb*w0*eta/(2*pi)/(nus*w0)^2/(E*1e9)/(sigma/c)^2;
        B = A + K*M;
        delta(:,ii) = eig(B);
    end
end
