function [Zl Zv Zh] = lnls_symb_calc_impedance_multilayer_round_pipe(w, epb, mur, ange, angm, sigmadc, tau, b, L)

old = digits;
digits(35);
c = 299792458;
mu0 = 4*pi*1e-7;
ep0 = 1/c^2/mu0;
gam = 3e9/511e3;
bet = sqrt(1-1/gam^2);

for i = 1: length(epb)
    ep1(i,:) = epb(i)*(1-1i.*sign(w).*tan(ange(i))) + sigmadc(i)./(1+1i*w*tau(i))./(1i*w*ep0);
    mu1(i,:) = mur(i)*(1-1i.*sign(w).*tan(angm(i)));
    nu(i,:) = abs(w/c).*sqrt(1 - bet^2*ep1(i,:).*mu1(i,:));
    nu2(i,:) = abs(w/c).^2.*(1 - bet^2*ep1(i,:).*mu1(i,:));
end

Zl = squeeze(1i*L*w/(2*pi*ep0*(bet*c)^2*gam^2) .*alphaTM(0, ep1, mu1, bet, nu, nu2, b));
Zv = squeeze(1i*L*w.^2/(4*pi*ep0*c^2*(bet*c)*gam^4) .*alphaTM(1, ep1, mu1, bet, nu, nu2, b));
Zh = Zv;
digits(old);



function alphaTM = alphaTM(m, ep1, mu1, bet, nu, nu2, b)

j=1;
for i = 1:length(b) % lembrando que length(b) = número de camadas - 1
    
    x = nu(i+1,:)*b(i);
    y = nu(i,:)*b(i);
    
    if i<length(b)
        z = nu(i+1,:)*b(i+1);
        
        if ~any(real(z)<0)
            
            ind = (real(z)<50);
            
            A = besseli(m,z(ind));
            B = besselk(m,z(ind));
            C = besseli(m,x(ind));
            E = besselk(m,x(ind));
            
            simbol{j}  = ['D' num2str(i) '1'];
            numero{j}  = 1; j = j+1;
            
            simbol{j}       = ['D' num2str(i) '2'];
            numero{j}(~ind) =  - exp(-2*(z(~ind)-x(~ind)));
            numero{j}(ind)  = -B.*C./(A.*E); j = j+1;
            
            A = sym(['D' num2str(i) '1']);
            Bi = sym(['D' num2str(i) '2']);
            D = diag([A Bi A Bi]);
            clear A Bi;
            
            
        end
    end
    
    numero{j} = -nu2(i+1,:)*b(i)./ep1(i+1,:).*( ep1(i+1,:)./nu(i+1,:).*(-besselk(m-1,x,1)./besselk(m,x,1) - m./x) ...
                                                - ep1(i,:)./nu(i,:).*(besseli(m-1,y,1)./besseli(m,y,1) - m./y));
    simbol{j}  = ['P' num2str(i) '1_1'];j = j+1;
                                            
    numero{j} = -nu2(i+1,:)*b(i)./ep1(i+1,:).*( ep1(i+1,:)./nu(i+1,:).*(-besselk(m-1,x,1)./besselk(m,x,1) - m./x) ...
                                                - ep1(i,:)./nu(i,:).*(-besselk(m-1,y,1)./besselk(m,y,1) - m./y));
    simbol{j}  = ['P' num2str(i) '1_2']; j = j+1;   
    
    numero{j} = -nu2(i+1,:)*b(i)./ep1(i+1,:).*( ep1(i+1,:)./nu(i+1,:).*(besseli(m-1,x,1)./besseli(m,x,1) - m./x) ...
                                                - ep1(i,:)./nu(i,:).*(besseli(m-1,y,1)./besseli(m,y,1) - m./y));
    simbol{j}  =  ['P' num2str(i) '2_1']; j = j+1;
    
    numero{j} = -nu2(i+1,:)*b(i)./ep1(i+1,:).*( ep1(i+1,:)./nu(i+1,:).*(besseli(m-1,x,1)./besseli(m,x,1) - m./x) ...
                                                - ep1(i,:)./nu(i,:).*(-besselk(m-1,y,1)./besselk(m,y,1) - m./y));
    simbol{j}  = ['P' num2str(i) '2_2']; j = j+1;
    
    numero{j} = (nu2(i+1,:)./nu2(i,:) - 1).*m./(bet*ep1(i+1,:));  
    simbol{j}  = ['Q' num2str(i)]; j = j+1;
    
    numero{j} = -nu2(i+1,:)*b(i)./mu1(i+1,:).*( mu1(i+1,:)./nu(i+1,:).*(-besselk(m-1,x,1)./besselk(m,x,1) - m./x) ...
                                                - mu1(i,:)./nu(i,:).*(besseli(m-1,y,1)./besseli(m,y,1) - m./y));
    simbol{j}  = ['R' num2str(i) '1_1'];j = j+1; 
    
    numero{j} = -nu2(i+1,:)*b(i)./mu1(i+1,:).*( mu1(i+1,:)./nu(i+1,:).*(-besselk(m-1,x,1)./besselk(m,x,1) - m./x) ...
                                                - mu1(i,:)./nu(i,:).*(-besselk(m-1,y,1)./besselk(m,y,1) - m./y));
    simbol{j}  = ['R' num2str(i) '1_2'];j = j+1;
    
    numero{j} = -nu2(i+1,:)*b(i)./mu1(i+1,:).*( mu1(i+1,:)./nu(i+1,:).*(besseli(m-1,x,1)./besseli(m,x,1) - m./x) ...
                                                - mu1(i,:)./nu(i,:).*(besseli(m-1,y,1)./besseli(m,y,1) - m./y));
    simbol{j}  = ['R' num2str(i) '2_1'];j = j+1;
    
    numero{j} = -nu2(i+1,:)*b(i)./mu1(i+1,:).*( mu1(i+1,:)./nu(i+1,:).*(besseli(m-1,x,1)./besseli(m,x,1) - m./x) ...
                                                - mu1(i,:)./nu(i,:).*(-besselk(m-1,y,1)./besselk(m,y,1) - m./y));
    simbol{j}  = ['R' num2str(i) '2_2'];j = j+1;
    
    numero{j} = (nu2(i+1,:)./nu2(i,:).^2 - 1).*m./(bet*mu1(i+1,:));  
    simbol{j}  = ['S' num2str(i)];j = j+1;
    
    P = sym(['P' num2str(i)],[2,2]);
    Q = sym(['P' num2str(i)],[2,2]); Q(:,:) = ['Q' num2str(i)];
    S = sym(['P' num2str(i)],[2,2]); S(:,:) = ['S' num2str(i)];
    R = sym(['R' num2str(i)],[2,2]);
    
    Mt = [P Q; S R];
    
    if length(b) == 1
        M = Mt;
    else
        if (i ==1)
            M = D*Mt;
        elseif i < length(b)
            M = D*(Mt*M);
        else
            M = Mt*M;
        end
    end
end

Bj = squeeze((M(1,2).*M(3,3) - M(3,2).*M(1,3)) ./ (M(1,1).*M(3,3) - M(1,3).*M(3,1)));
Bj = expand(Bj);
Bt = double(subs(Bj,simbol,numero));

alphaTM = besselk(m,nu(1,:)*b(1))./besseli(m,nu(1,:)*b(1)).* Bt;

