%% Pipeline for Calculations from Wakepotential Results

addpath(genpath(fileparts(mfilename('fullpath'))));

globdata = sirius_structures;

newdir = 'linear_taper_2';

globdata.simpar.wakepath = '/home/fernando/Windows/ECHO_Release2D_v_2_0/Release_2_0/chamBC';
rootdir = globdata.simpar.wakepath;

m = 1;
sigma = 0.5e-3;

mkdir(rootdir,newdir);
globdata.simpar.targetdir = fullfile(rootdir, newdir);
globdata.simpar.datasource = 'ECHO';

globdata.ringpar.sigma = sigma;
globdata.simpar.m = m;
globdata.simpar.sym = 1;
globdata.simpar.whichaxis = 'y';
globdata.simpar.units = 1;

% Load wakepotential result from referred software, rescale and save txt-file on default format
globdata = load_wake(globdata);

% Calculates Impedance Spectrum from Wakepotential Results
globdata = calc_impedance(globdata);

% Calculates Loss Factor 
if m == 0
    globdata = calc_loss_factor(globdata);
elseif m > 0
    globdata = calc_kick_factor(globdata);
end

% Plot Results 
plot_results(globdata);

% Export Results 
save_results(globdata);


