function [x_old, tau_old] = multi_bunch_tracking(ring, beam, wake)


RandStream.setGlobalStream(RandStream('mt19937ar','seed', 190488));

betx  = ring.beta;
etax  = ring.eta;
etxp  = ring.etap;
tune  = ring.tune;
chrom = ring.dtunedp;
tu_sh = ring.dtunedj;

fillPat_time = -(beam.fillPat/ring.har_num) * ring.rev_time;

x   = zeros(size(fillPat_time));
xp  = zeros(size(fillPat_time));
en  = zeros(size(fillPat_time));
tau = 1e-12*randn(size(fillPat_time)); %zeros(size(fillPat_time));
x_old   = zeros(100,length(fillPat_time));
tau_old = zeros(100,length(fillPat_time));

pl_tau = true;
pl_x   = false;

if pl_tau, figure; ax_tau = axes; plot(ax_tau,fillPat_time,tau,'b.', 'MarkerSize',1); end
if pl_x,  figure; ax_x  = axes; plot(ax_x ,fillPat_time,x, 'b.', 'MarkerSize',1); end
drawnow;
for ii=1:ring.nturns;
    % First do single particle tracking:
    % define one phase advance per particle
    phi = 2*pi*(tune + chrom*en + tu_sh*((x-etax*en).^2/betx + (xp-etxp*en).^2*betx));
    
    [x, xp] = transverse_tracking(x,xp,en,phi,betx,etax,etxp);
    % The longitudinal time evolution equations are not in the differential
    % form. Thus, any longitudinal potential can be taken into account.
    % I don't have to normalize the potential by the charge, because the
    % energy is already in electronVolts
    [en, tau] = longitudinal_tracking(en,tau, beam.tau,beam.potential,ring);
    
    [kickt, kickx, wake] = kick_from_wake(x/betx, tau, fillPat_time, wake, ii, ring, beam.I_b, beam.num_part);
    xp = xp + kickx;
    en = en + kickt;
    
    x_old   = circshift(x_old,  [1,0]);
    tau_old = circshift(tau_old,[1,0]);
    x_old(1,:) = x;
    tau_old(1,:) = tau;
    
    if mod(ii,1000)==0
        fprintf('%d\n',ii);
    end
    if mod(ii,10)==0
        if pl_x
            cury = get(ax_x,'YLim');
            nexty = [min([x,  cury(1)]), max([x,  cury(2)])];
            plot(ax_x,fillPat_time,x,'b.', 'MarkerSize',1);
            ylim(ax_x,nexty);
        end
        if pl_tau
            cury = get(ax_tau,'YLim');
            nexty = [min([tau, cury(1)]), max([tau, cury(2)])];
            plot(ax_tau,fillPat_time,tau,'b.', 'MarkerSize',1);
            ylim(ax_tau,nexty);
        end
        drawnow;
    end
end


function [x_new, xp_new] = transverse_tracking(x,xp,en,phix,betx,etax,etxp)

x_new  =         (x-etax*en).*cos(phix) + betx*(xp-etxp*en).*sin(phix) + etax*en;
xp_new = -1/betx*(x-etax*en).*sin(phix) +      (xp-etxp*en).*cos(phix) + etxp*en;


function [en,tau] = longitudinal_tracking(en,tau,ref_tau,potential, ring)

for i = 1:length(en)
    en(i) = en(i) + interp1q(ref_tau', potential(i,:)', tau(i))/ring.E;
end
tau = tau - ring.rev_time.*(en*ring.mom_comp);


function [kickt, kickx, wake] = kick_from_wake(x, tau, fillPat_time, wake, volta, ring, I_b, num_part)

kickx  = zeros(size(x));
kickt = zeros(size(tau));
% Remember that tau is the position ahead of the synchronous particle.
% Thus, a positive tau passes through the impedance before a negative tau.

tau = tau + fillPat_time - volta/ring.rev_time;

if wake.dipo.track && isfield(wake.dipo,'general')
    difft = bsxfun(@minus,tau,tau');
    kik = interp1(wake.dipo.general.tau,wake.dipo.general.wake,difft,'linear',0);
    kickx = kickx - (x*kik);
end
if wake.dipo.track && isfield(wake.dipo,'resonator')
    Memory = wake.dipo.resonator.Memory;
    wr = wake.dipo.resonator.wr;
    Rs = wake.dipo.resonator.Rs;
    Q  = wake.dipo.resonator.Q;
    betax = wake.dipo.resonator.beta;
    Ql = sqrt(Q.^2 - 1/4);
    wrl = wr .* Ql ./ Q;
    
    W_pot = zeros(1,length(x)+1);
    for ii=1:length(wr)
        W_pot(1) = 0;
        W_pot(2:end) = cumsum((ring.rev_time/ring.E) * (I_b / num_part).* ...
            exp(-tau*(1i*wrl(ii)+wr(ii)/(2*Q(ii))))*betax(ii).*x);
        kik = (W_pot(1:end-1) + Memory(ii)) .* ...
            exp( tau*(1i*wrl(ii)+wr(ii)/(2*Q(ii))));
        
        % Memory of the cavity
        Memory(ii) = Memory(ii) + W_pot(end);
        
        kickx = kickx - wr(ii)*Rs(ii)/Ql(ii) * imag(kik);
        %         % This is the longitudinal kick given by the transverse wake:
        %         % Chao's book, p. 212
        %         kickt(I) = kickt(I) - betax(ii)*x(I).* wr(ii)*Rs(ii)/Q(ii).*(real(kik) + 1/(2*Ql(ii))*imag(kik));
    end
    wake.dipo.resonator.Memory = Memory;
end


if  wake.quad.track && isfield(wake.quad,'general')
    if ~exist('difft','var'), difft = bsxfun(@minus,tau,tau'); end
    kik = interp1(wake.quad.general.tau,wake.quad.general.wake,difft,'linear',0);
    kickx = kickx - (sum(kik).*x);
end
if wake.quad.track && isfield(wake.quad,'resonator')
    Memory = wake.quad.resonator.Memory;
    wr = wake.quad.resonator.wr;
    Rs = wake.quad.resonator.Rs;
    Q  = wake.quad.resonator.Q;
    betax = wake.quad.resonator.beta;
    
    Ql = sqrt(Q.^2 - 1/4);
    wrl = wr .* Ql ./ Q;
    
    W_pot = zeros(1,length(x)+1);
    for ii=1:length(wr)
        W_pot(1) = 0;
        W_pot(2:end) = cumsum((ring.rev_time/ring.E) * (I_b / num_part).* ...
            exp(-tau*(1i*wrl(ii)+wr(ii)/(2*Q(ii)))));
        kik = (W_pot(1:end-1) + Memory(ii)) .* exp( tau*(1i*wrl(ii)+wr(ii)/(2*Q(ii))));
        
        % Memory of the cavity
        Memory(ii) = Memory(ii) + W_pot(end);
        
        kickx = kickx - x .* betax(ii)* wr(ii)*Rs(ii)/Ql(ii) .* imag(kik);
        
        %         % This is the longitudinal kick given by the transverse wake:
        %         % Not in Chao's book, but easily derived
        %         kickt(I) = kickt(I) - (betax(ii)*x(I)).^2.* wr(ii)*Rs(ii)/Q(ii).*(real(kik) + 1/(2*Ql)*imag(kik));
    end
    wake.quad.resonator.Memory = Memory;
end


if wake.long.track && isfield(wake.long,'general')
    if ~exist('difft','var'), difft = bsxfun(@minus,tau,tau'); end
    kik = interp1(wake.long.general.tau,wake.long.general.wake,difft,'linear',0);
    kickt = kickt - sum(kik);
end
if wake.long.track && isfield(wake.long,'resonator')
    Memory = wake.long.resonator.Memory;
    wr = wake.long.resonator.wr;
    Rs = wake.long.resonator.Rs;
    Q  = wake.long.resonator.Q;
    
    Ql = sqrt(Q.^2 - 1/4);
    wrl = wr .* Ql ./ Q;
    
    W_pot = zeros(1,length(tau)+1);
    for ii=1:length(wr)
        W_pot(1) = 0;
        W_pot(2:end) = cumsum( (ring.rev_time/ring.E) * (I_b / num_part).* ...
            exp(-tau*(1i*wrl(ii)+wr(ii)/(2*Q(ii)))));
        kik = (W_pot(1:end-1) + Memory(ii)) .* exp(tau*(1i*wrl(ii)+wr(ii)/(2*Q(ii))));
        
        % Memory of the cavity
        Memory(ii) = Memory(ii) + W_pot(end);
        
        kickt = kickt - wr(ii)*Rs(ii)/Q(ii) * (1/2 + real(kik) + 1/(2*Ql)*imag(kik));
    end
    wake.long.resonator.Memory = Memory;
end
