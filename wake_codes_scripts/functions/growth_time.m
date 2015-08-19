
% script for calculating growth rates/times based on impedance data (Fc, Rs and Q)
% references: 
% 1) ..\HC_SEM\Alexei_Meetings\Alexei\2012-07-10_recvd_LKickerLNLS
% (T_Wake, T_Coupl_bunch).ppt
% 2) ..\HC_SEM\Scripts\Alexei\CBThresholdzLN.nb

function  alt = growth_time(globdata)  
    %% Loading and calculating parameters and results
    alt = globdata; 
    delta = 0;         %  center frequency shift (MHz)

    % basic selections

    omegar = alt.results.peakinfo.omegar + delta*2e6*pi;
    Qs = alt.results.peakinfo.Q;
    Rshs = alt.results.peakinfo.Rshunt;
    RovQ = Rshs./Qs;

    max_HOM = length(Qs);                 %  max HOM to be taken in the calculations
    
    sigs = alt.ringpar.sigma;
    Cecam = alt.ringpar.circumf;
    R = Cecam/(2*pi);
    bettaref = alt.ringpar.beta;
    cvel = 299792458;
    E0 = alt.ringpar.E0;
    M = alt.ringpar.h;
    Iav = alt.ringpar.Iavg;
    alpha = alt.ringpar.alpha;
    omega0 = alt.ringpar.omega0;

    % bunch length (s)
    sigt = sigs*cvel*bettaref;
    % electron rest energy (eV)
    Ee = 0.510998910e6;
    % relativistic gamma
    gam = E0/Ee;
    % average bunch current
    Ib = Iav/M;
    % Revolution Period-Frequency
    T0 = 2*pi/omega0;
    % Horiz. tune
    nux = alt.ringpar.nux;
    omegab = nux*omega0;
    % Sincrotron tune
    omegas = alt.ringpar.omegas;
    nus = omegas/omega0;

    % equation coefficient - ref: ..\HC_SEM\Documents\Formulae\growth_rate.docx
    Wcoef = Ib*M*omega0^2*alpha/(4*pi*E0*omegas);

    % fundamental charge
    e= 1.6e-19;
    % average bunch charge
    Ne = Iav*T0/(M*e);
    % momentum spread
    sigp = nus*sigs/(alpha*R);

    % x scale for impedance and y initialization
    wfull = linspace(-1.2*max(omegar),1.2*max(omegar),1e6); 
    Zbfull = zeros(length(wfull),1);

%% Modeling Impedance Spectrum According to Modes Info

    for m=1:1:length(Qs)  % all HOMs
        for i=1:1:length(wfull)   % whole spectrum

            Zbfull(i) = Zbfull(i) + Rshs(m) / (1 + j*Qs(m)*(wfull(i)/omegar(m) - omegar(m)/wfull(i)));    

        end
    end
    
    whalf = wfull(floor(length(wfull)/2):end);
    Zbhalf = Zbfull(floor(length(wfull)/2):end);
    
    
%% Calculates growth time/rate for modes up to 100th rf harmonic

    p = -100:1:100;

    iCBM = 0:1:M-1; 
    GRs = zeros(M,1);
    GR_HOM = zeros(M,max_HOM);

    for is = 1:1:length(iCBM) % all CBMs
        for m = 1:1:max_HOM % all HOMs
            for i=1:1:length(p) % all RF harmonics
                ReZlocal = 0;

                ReZlocal = real( Rshs(m) / (1 + j*Qs(m)*(((p(i)*M)+iCBM(is)+nus)*omega0/omegar(m) - omegar(m)/(((p(i)*M)+iCBM(is)+nus)*omega0))) );
                calc = Wcoef*(((p(i)*M) + iCBM(is))^2)*exp(-(((p(i)*M) + iCBM(is))^2)*((sigs/R)^2))*ReZlocal/(((p(i)*M)+iCBM(is)+nus));

                GRs(is) = GRs(is) + calc;
            end
            GR_HOM(is,m) = GRs(is);
        end
    end

    % gets the peak growth rate coordinate
    [aux,ifast] = max(GRs);

%% Samples impedance at the CBM frequencies

    Zsampl = zeros(length(p),1);

    for m=1:1:length(Qs);  % all HOMs
        for i=1:1:length(p)  % whole spectrum     
            Zsampl(i) = Zsampl(i) + Rshs(m) / (1 + j*Qs(m)*(((p(i)*M)+(ifast-1)+nus)*omega0/omegar(m) - omegar(m)/(((p(i)*M)+(ifast-1)+nus)*omega0)));   
        end
    end

    for i=1:1:length(p)    
        fsampl(i) = abs((((p(i)*M)+(ifast-1)+nus)*omega0)/(2e9*pi));     
    end

%% Results assignement to structure

    alt.results.peakinfo.wmodel = whalf;
    alt.results.peakinfo.ReZmodel = real(Zbhalf);
    
    alt.results.GRs = GRs;
    alt.results.GR_HOM = GR_HOM;
    alt.results.ifast = ifast;
    alt.results.ReZsampl = real(Zsampl);
    alt.results.fsampl = fsampl;
    

end
