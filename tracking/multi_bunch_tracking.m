function [ave_bun, rms_bun] = multi_bunch_tracking(ring, beam, wake)

ave_bun = zeros(4,ring.nturns);
rms_bun = zeros(4,ring.nturns);

RandStream.setGlobalStream(RandStream('mt19937ar','seed', 190488));
%definicao dos pacotes
%generate the longitudinal phase space;

betx  = ring.beta;
alpx  = ring.alpha;
etax  = ring.eta;
etxp  = ring.etap;
tune  = ring.tune;
chrom = ring.dtunedp;
tu_sh = ring.dtunedj;

npoints = ring.har_num^2*wake.nturns;
time_struct = -((1:npoints) - 1) * ring.rev_time/ring.har_num;


figure;
for ii=1:ring.nturns;
    phi  =  2*pi*(tune + chrom*en + tu_sh*((x-etax*en).^2/betx + ...
                 ((xp-etxp*en).^2 + alpx/betx*(x-etax*en).^2)*betx));
    
    [x, xp] = transverse_tracking(x,xp,en,phi,betx,alpx,etax,etxp);
    en      = en  - interp1(beam.tau, beam.potential, tau)/ring.E;
    tau     = tau - ring.rev_time.*(en*ring.mom_comp);
    
    kickx = kickx_from_wake(x, tau, wake, ii, ring.E, beam.I_b, ring.rev_time);
    kickt = kickt_from_wake(tau, wake, ii, ring.E, beam.I_b, ring.rev_time);
    xp = xp + kickx/betx;
    en = en + kickt;
    ave_bun(:,ii) = mean([x;xp;en;tau],2);
    rms_bun(:,ii) =  rms([x;xp;en;tau],2);
    if mod(ii,100)==0
        fprintf('%d\n',ii);
    end
end

figure; plot(squeeze(ave_bun(1,:)));


function [x_new, xp_new] = transverse_tracking(x,xp,en,phix,betx,alpx,etax,etxp)

x_new  =         (x-etax*en).*cos(phix) + betx*(xp-etxp*en).*sin(phix) + etax*en;
xp_new = -1/betx*(x-etax*en).*sin(phix) +      (xp-etxp*en).*cos(phix) + etxp*en;


function kick = kickx_from_wake(x, tau, wake, volta, E, I_b, T0)

kick = zeros(size(x));
np   = length(tau);

difft = bsxfun(@minus,tau',tau);
if wake.sing.dipo.sim
    kik = interp1(wake.sing.dipo.tau,wake.sing.dipo.wake,difft,'linear',0);
    kick = kick + T0*I_b/E/np*(x*kik);
end
if wake.sing.quad.sim
    kik = interp1(wake.sing.quad.tau,wake.sing.quad.wake,difft,'linear',0);
    kick = kick + T0*I_b/E/np*(sum(kik).*x);
end
if wake.mult.trans.sim
    kick = kick + wake.mult.trans.wake(volta);
end

function kick = kickt_from_wake(tau, wake, volta, E, I_b, T0)

kick = zeros(size(tau));
np   = length(tau);

if wake.sing.long.sim
    difft = bsxfun(@minus,tau',tau);
    kik = interp1(wake.sing.long.tau,wake.sing.long.wake,difft,'linear',0);
    kick = kick + T0*I_b/E/np*sum(kik);
end
if wake.mult.long.sim
    kick = kick + wake.mult.long.wake(volta);
end