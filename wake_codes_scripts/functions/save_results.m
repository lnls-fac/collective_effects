function  save_results(globdata)
% Function for exporting results requested by the user on txt file
%  according to the nomenclature below:
%
%   1 - Wakepotential: W<wtype><codename>.txt
%   2 - Real Part of Impedance: Re<wtype><codename>.txt
%   3 - Imaginary Part of Impedance: Re<wtype><codename>.txt
%   4 - Loss Parameters: Loss info_<codename>.txt
%   5 - Loss Factor per sigma: Kloss<codename>.txt
%
%
% where
%   <wtype> = long
%             xdip or xquad
%             ydip or yquad
%
%   <codename> = GdfidL
%                CST
%                ACE3P
%                ECHO
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filesout = globdata.simpar.targetdir;
dsrc     = globdata.simpar.datasource;
m        = globdata.simpar.m;

if m==0
    wtype = 'long';
elseif m==1
    wtype = [globdata.simpar.whichaxis 'dip'];
elseif m==2
    wtype = [globdata.simpar.whichaxis 'quad'];
end

wake = globdata.results.W;
spos = globdata.results.s;

if m==0
    rez = globdata.results.ReZlong;
    imz = globdata.results.ImZlong;
elseif m>0
    rez = globdata.results.ReZt/1000;
    imz = globdata.results.ImZt/1000;
end

ImZoN = globdata.results.ImZoN;
f     = globdata.results.freq;
naxis = globdata.results.naxis;
sigi  = globdata.results.sigmak;

kZi = globdata.results.klossZ;
kW  = globdata.results.klossW;

kickZi = globdata.results.kickZ;
kickW  = globdata.results.kickW;

Ploss = globdata.results.Ploss;
T0    = 2*pi/globdata.ringpar.omega0;

%% Tick Position # 2: Export wakepotential
fid = fopen(fullfile(filesout, ['W', wtype, dsrc, '.txt']), 'w');
fprintf(fid, '%30.16g %30.16g\n', [spos; wake]);
fclose(fid);

%% Tick Position # 5: Export Impedance
fid = fopen(fullfile(filesout,['ReZ', wtype, dsrc, '.txt']), 'w');
fprintf(fid,'%30.16g    %30.16g\n',[f; rez]);
fclose(fid);

fid = fopen(fullfile(filesout,['ImZ', wtype, dsrc, '.txt']), 'w');
fprintf(fid,'%30.16g    %30.16g\n',[f; imz]);
fclose(fid);
if m==0
    fid = fopen(fullfile(filesout,['ImZoN', wtype, dsrc, '.txt']), 'w');
    fprintf(fid,'%30.16g    %30.16g\n',[naxis; ImZoN]);
    fclose(fid);
end

%% Tick Position # 8: Export Loss Factor vs. Sigma and Loss Info
if m==0
    fid = fopen(fullfile(filesout,['Loss info_', dsrc, '.txt']), 'w');
    fprintf(fid,'Loss factor Z = %10.6f mV/pC  \nLoss factor W = %10.6f mV/pC  \nPower Loss = %10.5f W \nfor I = %9.4f mA ; h =%5.0f ; T0 = %8.4f ns ', kZi(1)*1e3, kW*1e3, Ploss, globdata.ringpar.Iavg*1e3, globdata.ringpar.h, T0*1e9);
    fclose(fid);
    
    fid = fopen(fullfile(filesout,['Kloss', dsrc, '.txt']), 'w');
    fprintf(fid,'%12.8f    %12.8f\n', [sigi./1e-3 ; abs(kZi)]);
    fclose(fid);
elseif m>0
    fid = fopen(fullfile(filesout,['Kick info_', wtype, dsrc, '.txt']), 'w');
    fprintf(fid, 'Kick Z = %10.6f V/pC/m  \nKick W = %10.6f V/pC/m', kickZi(1), kickW);
    fclose(fid);
    
    fid = fopen(fullfile(filesout,['K', wtype, dsrc, '.txt']), 'w');
    fprintf(fid, '%12.8f    %12.8f\n', [sigi./1e-3 ; abs(kickZi)]);
    fclose(fid);
end

save(fullfile(filesout,'globdata.mat'),'globdata');
