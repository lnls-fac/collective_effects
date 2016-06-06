function  globdata = load_wake(globdata)
% Function to read wakepotential files from any software and:
% - Rescale s axis to m (for CST and ECHO)
% - ACE3P (exclusive): Performs displacement of -Nsigmas*sigma over s axis
% - Exports new file with standard nomenclature: W<waketype><codename>.txt
% - Plot Wakepotential Result

nsigmas = 5; % Only used for ACE3P displacement, standard value!

dsrc   = globdata.simpar.datasource;
m      = globdata.simpar.m;
sym    = globdata.simpar.sym;
whaxis = globdata.simpar.whichaxis;
wdir   = globdata.simpar.wakepath;
tardir = globdata.simpar.targetdir;
%% Each Software uses a different nomenclature. Setting complete path for loading:

if strcmp(dsrc,'ACE3P')
    wakepath = fullfile(wdir,'wakefield.out');
    headerL = 3;
    declar = '%f%f%f';
elseif strcmp(dsrc,'GdfidL')
    declar = '%f%f';
    headerL = 11;
    
    if m==0
        wakepath = fullfile(wdir,'Results-Wq_AT_XY.0001');
    elseif m==1
        if sym
            shiftpath = fullfile(wdir,'Results-Wq_AT_XY.0001');
            wakepath1 = fullfile(wdir,['Results-W' upper(whaxis) '_AT_XY.0001']);
            wakepath2 = fullfile(wdir,['Results-W' upper(whaxis) '_AT_XY.0002']);
        else
            shiftpath = fullfile(wdir,'dxdpl','Results-Wq_AT_XY.0001');
            wakepath1 = fullfile(wdir,['d', whaxis, 'dpl'],['Results-W', upper(whaxis), '_AT_XY.0001']);
            wakepath2 = fullfile(wdir,['d', whaxis, 'dpl'],['Results-W', upper(whaxis), '_AT_XY.0002']);
            wakepath3 = fullfile(wdir,['d', whaxis, 'dmi'],['Results-W', upper(whaxis), '_AT_XY.0001']);
            wakepath4 = fullfile(wdir,['d', whaxis, 'dmi'],['Results-W', upper(whaxis), '_AT_XY.0002']);
        end
    elseif m==2
        wakepath1 = fullfile(wdir,['Results-W', upper(whaxis), '_AT_XY.0001']);
        if ~sym
            wakepath2 = fullfile(wdir,['Results-W', upper(whaxis), '_AT_XY.0002']);
        end
    end
elseif strcmp(dsrc,'CST')
    wakepath = fullfile(wdir,'wake.txt');
    headerL = 2;
    declar = '%f%f';
elseif strcmp(dsrc,'ECHO')
    if m==0
        wakepath = fullfile(wdir,'wake.dat');
        headerL = 0;
        declar = '%f%f';
    elseif m > 0
        if sym
            wakepath = fullfile(wdir,'wakeT.dat');
            headerL = 0;
            declar = '%f%f';
        else
            wrt = 'symetry error: wrong set';
        end
    end
end

%% Read specified file(s)

if (m==0 || ~strcmp(dsrc,'GdfidL'))
    fid = fopen(wakepath,'rt');
    loadres = textscan(fid, declar, 'HeaderLines', headerL);
    fclose(fid);
    copyfile(wakepath, tardir);
    
    spos = transpose(loadres{1});
    wake = transpose(loadres{2});
elseif m==1
    fid1 = fopen(wakepath1,'rt');
    loadres1 = textscan(fid1, declar, 'HeaderLines', headerL);
    fclose(fid1);
    copyfile(wakepath1, tardir);
    
    fid2 = fopen(wakepath2,'rt');
    loadres2 = textscan(fid2, declar, 'HeaderLines', headerL);
    fclose(fid2);
    copyfile(wakepath2, tardir);
    
    spos = transpose(loadres1{1});
    wabs = transpose(loadres1{2}+loadres2{2})/2;
    
    if ~sym
        fid3 = fopen(wakepath3,'rt');
        loadres3 = textscan(fid3, declar, 'HeaderLines', headerL);
        fclose(fid3);
        copyfile(wakepath3, tardir);
        
        fid4 = fopen(wakepath4,'rt');
        loadres4 = textscan(fid4, declar, 'HeaderLines', headerL);
        fclose(fid4);
        copyfile(wakepath4, tardir);
        
        wabs2 = transpose(loadres3{2}+loadres4{2})/2;
        wabs = (wabs - wabs2)/2;
    end
    
    %obtaining shift value
    fid5 = fopen(shiftpath,'rt');
    for i=1:3
        loadres5 = fgetl(fid5);
    end
    fclose(fid5);
    copyfile(shiftpath, tardir);
    
    coord = textscan(loadres5,' %% subtitle= "W_l (x,y)= ( %f, %f ) [m]"');
    
    if strcmp(whaxis,'x')
        shift = coord{1};
    elseif strcmp(whaxis,'y')
        shift = coord{2};
    end
    
    wake = wabs/shift;
elseif m==2
    fid1 = fopen(wakepath1,'rt');
    loadres1 = textscan(fid1, declar, 'HeaderLines', headerL);
    fclose(fid1);
    copyfile(wakepath1, tardir);
    
    spos = transpose(loadres1{1});
    wabs = transpose(loadres1{2});
    
    if ~sym
        fid2 = fopen(wakepath2,'rt');
        loadres2 = textscan(fid2, declar, 'HeaderLines', headerL);
        fclose(fid2);
        copyfile(wakepath2, tardir);
        
        w2 = transpose(loadres2{2});
        wabs = (wabs - w2)/2;
    end
    
    %obtaining offset value
    fid3 = fopen(wakepath1,'rt');
    for i=1:3
        loadres5 = fgetl(fid3);
    end
    fclose(fid3);
    copyfile(wakepath1, tardir);
    
    if strcmp(whaxis,'x')
        coord = textscan(loadres5,' %% subtitle= "integral d/dx W(z) dz, (x,y)=( %f, %f )"');
        shift = coord{1};
    elseif strcmp(whaxis,'y')
        coord = textscan(loadres5,' %% subtitle= "integral d/dy W(z) dz, (x,y)=( %f, %f )"');
        shift = coord{2};
    end
    
    wake = -wabs/shift;
end

%% Adjust s-axis (rescale or shift)

if strcmp(dsrc,'ACE3P')
    spos = spos - nsigmas*globdata.simpar.sigma; % Performs displacement over s axis
elseif strcmp(dsrc,'CST')
    spos = spos/1000; % Rescaling mm to m
    if m>0
        wake = -wake;
    end
elseif strcmp(dsrc,'ECHO')
    spos = spos/100; % Rescaling cm to m
end

%% Assign to Structure
globdata.results.W = wake;
globdata.results.s = spos;