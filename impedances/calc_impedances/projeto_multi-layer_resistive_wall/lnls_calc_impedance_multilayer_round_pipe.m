function [Zl Zv Zh] = lnls_calc_impedance_multilayer_round_pipe(w, epr, mur, b, L,E)%, filt, n, lim_w)


c = 299792458;
mu0 = 4*pi*1e-7;
ep0 = 1/c^2/mu0;
gam = E*1e9/511e3;
bet = sqrt(1-1/gam^2);


nu = (ones(size(epr,1),1)*abs(w/c)).*sqrt(1 - bet^2*epr.*mur);
nu2 = (ones(size(epr,1),1)*abs(w/c)).^2.*(1 - bet^2*epr.*mur);


Zl = -squeeze(1i*L*w/(2*pi*ep0*(bet*c)^2*gam^2) .*alphaTM(0, epr, mur, bet, nu, nu2, b)); % minus signal to compatibilize with Chao's convention
Zv = -squeeze(1i*L*w.^2/(4*pi*ep0*c^2*(bet*c)*gam^4) .*alphaTM(1, epr, mur, bet, nu, nu2, b));% minus signal to compatibilize with Chao's convention

% The code cant handle w = 0;
ind0 = find(w == 0);
if ~isempty(ind0)
    indht0 = find(w > 0);
    [~,indht0s] = sort(w(indht0));
    Zl(ind0) = real(Zl(indht0(indht0s(1))));
    Zv(ind0) = 0 + 1i*imag(Zv(indht0(indht0s(1))));
end

Zh = Zv;



function alphaTM = alphaTM(m, epr, mur, bet, nu, nu2, b)

M = Mtil(m, epr, mur, bet, nu, nu2, b);

B = squeeze((M(1,2,:).*M(3,3,:) - M(3,2,:).*M(1,3,:)) ./ (M(1,1,:).*M(3,3,:) - M(1,3,:).*M(3,1,:)));
alphaTM = besselk(m,nu(1,:)*b(1))./besseli(m,nu(1,:)*b(1)).* B';




function M = Mtil(m, epr, mur, bet, nu, nu2, b)



for i = 1:length(b) % lembrando que length(b) = número de camadas - 1
   
    x = nu(i+1,:)*b(i);
    y = nu(i,:)*b(i);
    
    if i<length(b)
        
        D = zeros(4,4,size(nu,2));
        
        z = nu(i+1,:)*b(i+1);
        
        if ~any(real(z)<0)
            
            ind = (real(z)<50);
            
            A = besseli(m,z(ind));
            B = besselk(m,z(ind));
            C = besseli(m,x(ind));
            E = besselk(m,x(ind));
            
            D(1,1,:) =  1;
            D(2,2,ind) = -B.*C./(A.*E); D(2,2,~ind) =  - exp(-2*(z(~ind)-x(~ind)));
            D(3,3,:) =  1;
            D(4,4,ind) = -B.*C./(A.*E); D(4,4,~ind) =  - exp(-2*(z(~ind)-x(~ind)));
        end
    end
    
    Mt(1,1,:) = -nu2(i+1,:)*b(i)./epr(i+1,:).*( epr(i+1,:)./nu(i+1,:).*(-besselk(m-1,x,1)./besselk(m,x,1) - m./x) ...
                                                - epr(i,:)./nu(i,:).*(besseli(m-1,y,1)./besseli(m,y,1) - m./y));
            
    Mt(1,2,:) = -nu2(i+1,:)*b(i)./epr(i+1,:).*( epr(i+1,:)./nu(i+1,:).*(-besselk(m-1,x,1)./besselk(m,x,1) - m./x) ...
                                                - epr(i,:)./nu(i,:).*(-besselk(m-1,y,1)./besselk(m,y,1) - m./y));
            
    Mt(2,1,:) = -nu2(i+1,:)*b(i)./epr(i+1,:).*( epr(i+1,:)./nu(i+1,:).*(besseli(m-1,x,1)./besseli(m,x,1) - m./x) ...
                                                - epr(i,:)./nu(i,:).*(besseli(m-1,y,1)./besseli(m,y,1) - m./y));
            
    Mt(2,2,:) = -nu2(i+1,:)*b(i)./epr(i+1,:).*( epr(i+1,:)./nu(i+1,:).*(besseli(m-1,x,1)./besseli(m,x,1) - m./x) ...
                                                - epr(i,:)./nu(i,:).*(-besselk(m-1,y,1)./besselk(m,y,1) - m./y));
    
    Mt(1,3,:) = (nu2(i+1,:)./nu2(i,:) - 1).*m./(bet*epr(i+1,:));  
    Mt(1,4,:) = Mt(1,3,:);
    Mt(2,3,:) = Mt(1,3,:); 
    Mt(2,4,:) = Mt(1,3,:); 
    
    Mt(3,3,:) = -nu2(i+1,:)*b(i)./mur(i+1,:).*( mur(i+1,:)./nu(i+1,:).*(-besselk(m-1,x,1)./besselk(m,x,1) - m./x) ...
                                                - mur(i,:)./nu(i,:).*(besseli(m-1,y,1)./besseli(m,y,1) - m./y));
            
    Mt(3,4,:) = -nu2(i+1,:)*b(i)./mur(i+1,:).*( mur(i+1,:)./nu(i+1,:).*(-besselk(m-1,x,1)./besselk(m,x,1) - m./x) ...
                                                - mur(i,:)./nu(i,:).*(-besselk(m-1,y,1)./besselk(m,y,1) - m./y));
            
    Mt(4,3,:) = -nu2(i+1,:)*b(i)./mur(i+1,:).*( mur(i+1,:)./nu(i+1,:).*(besseli(m-1,x,1)./besseli(m,x,1) - m./x) ...
                                                - mur(i,:)./nu(i,:).*(besseli(m-1,y,1)./besseli(m,y,1) - m./y));
            
    Mt(4,4,:) = -nu2(i+1,:)*b(i)./mur(i+1,:).*( mur(i+1,:)./nu(i+1,:).*(besseli(m-1,x,1)./besseli(m,x,1) - m./x) ...
                                                - mur(i,:)./nu(i,:).*(-besselk(m-1,y,1)./besselk(m,y,1) - m./y));
    
    Mt(3,1,:) = (nu2(i+1,:)./nu2(i,:) - 1).*m./(bet*mur(i+1,:));  
    Mt(3,2,:) = Mt(3,1,:);
    Mt(4,1,:) = Mt(3,1,:); 
    Mt(4,2,:) = Mt(3,1,:); 
    
    if length(b) == 1
        M = Mt;
    else
        if (i ==1) 
            M = produto(D,Mt);
        elseif i < length(b)
            M = produto(D,produto(Mt,M));
        else
            M = produto(Mt,M);
        end
    end
    
end
    
function C = produto(A,B)
C = zeros(size(A));
for i =1:size(A,1)
    for j = 1:size(B,2)
        for k = 1:size(A,2)
            C(i,j,:) = C(i,j,:) + A(i,k,:).*B(k,j,:);
        end
    end
end
