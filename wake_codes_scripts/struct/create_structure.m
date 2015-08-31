%   SIRIUS_PARAMETERS
%
%   Type "edit sirius_parameters" to edit parameters.

function wfpstruct = create_structure() 

        globdata = struct( ...
            'ringpar' , [], ...                   % Storage Ring Parameters
            'guipar', [], ...                     % GUI parameters
            'simpar', [], ...                     % Simulation parameters
            'results', [] ...                    % Script Results
        );

        globdata.ringpar = struct(...
            'circumf', 518.246, ...               % Ring circumference [m]
            'omega0', 3.63466682e6, ...           % Revolution frequency [rad/s]
            'sigma', 2.65e-3, ...                % Bunch Length [m]
            'Iavg', 500e-3, ...                   % Beam current [A]
            'h', 864, ...                         % Harmonic number 
            'beta', 0.9999999855, ...             % Beam speed (fraction of c)
            'E0', 3e9, ...                        % Beam Energy [eV]
            'Ee', 0.510998910e6, ...              % Electron Rest Energy [eV]
            'alpha', 1.74e-4, ...                 % Momentum Compaction Factor
            'nux', 16.941, ...                    % Horizontal Tune
            'omegas', 2*pi*2.61e3 ...             % Synchrotron Frequency [rad/s]
        );

        globdata.guipar = struct( ...
            'procoptions', [], ...                % Processing Options (uint64)
            'message', [] ...                     % GUI Message box 
        );


        globdata.simpar = struct( ...
            'wakepath' , [], ...                  % Path where the wake file (from any of the softwares) is
            'targetdir' , [], ...                 % Path where all results will be saved
            'datasource', [], ...                 % CST, ACE3P, ECHO or GdfidL
            'waketype', [],  ...                  % long, y or x
            'units', [] ...                       % # of components in the ring
        );


        globdata.results = struct( ...
            's' , [], ...                         % axis: distance from following to drive bunch [m]
            'Wlong' , [], ...                     % Longitudinal Wakepotential [V/pC]
            'Wy', [], ...                         % Vertical Dipole Wakepotential [V/pC/m]
            'Wx', [], ...                         % Vertical Dipole Wakepotential [V/pC/m]
            'freq', [], ...                       % axis: frequency obtained from FFT [GHz]
            'naxis', [], ...                      % axis: omega/omega0
            'interfreq', [], ...                  % interpolated frequency
            'ReZlong' , [], ...                   % Real Part of Longitudinal Impedance [Ohm]
            'interReZlong', [], ...               % interpolated impedance
            'ImZlong' , [], ...                   % Imaginary Part of Longitudinal Impedance [Ohm]
            'ImZoN' , [], ...                     % Imaginary Part of Longitudinal Impedance over n [Ohm]
            'ReZy' , [], ...                      % Real Part of Vertical Dipole Impedance [KOhm/m]
            'interReZy', [], ...                  % interpolated impedance
            'ImZy' , [], ...                      % Imaginary Part of Vertical Dipole Impedance [KOhm/m]
            'ReZx' , [], ...                      % Real Part of Horizontal Dipole Impedance [KOhm/m]
            'interReZx', [], ...                  % interpolated impedance
            'ImZx' , [], ...                      % Imaginary Part of Horizontal Dipole Impedance [KOhm/m]
            'peakinfo', [], ...'                  % Omegar [rad/s], Rshunt [Ohm] and Q from ReZ 
            'klossW' , [], ...                    % Single-bunch Loss Factor Calculated from Wlong [mV/pC]
            'klossZ' , [], ...                    % Single-bunch Loss Factor Calculated from ReZlong [mV/pC]
            'kyW' , [], ...                       % Vertical Kick Factor Calculated from Wy [V/pC/m]
            'kyZ' , [], ...                       % Vertical Kick Factor Calculated from ImZy [V/pC/m]
            'kxW' , [], ...                       % Horizontal Kick Factor Calculated from Wx [V/pC/m]
            'kxZ' , [], ...                       % Horizontal Kick Factor Calculated from ImZx [V/pC/m]
            'sigmak' , [], ...                    % axis: bunch length for kloss|kick integration [mm]
            'Ploss', [], ...                      % Power loss from single-bunch loss factor [W]
            'Plossvec', [], ...                   % Power loss vector from single-bunch loss factor, for different sigmas [W]
            'klossWM' , [], ...                   % Multi-bunch Loss Factor Calculated from Wlong [mV/pC]
            'klossZM' , [], ...                   % Multi-bunch Loss Factor Calculated from ReZlong [mV/pC]
            'PlossM', [], ...                      % Power loss from multi-bunch loss factor [W]
            'ifast', [], ...                      % # of fastest CBM
            'GRs', [], ...                        % Growth Rate value for each CBM
            'GR_HOM', [], ...                     % Growth Rate value for each CBM accumulated through each HOM
            'ReZsampl', [], ...                   % Impedance Spectrum Sampled by fastest CBM
            'fsampl', [] ...                     % Frequency axis for sampled impedance
        );

        globdata.results.peakinfo = struct( ...
            'omegar', [], ...                     % Modes Center Frequencies
            'BW', [], ...                         % Modes Bandwidths (half height)
            'Rshunt', [], ...                     % Modes Peak Values
            'Q', [], ...                          % Modes Quality Factors
            'ReZmodel', [], ...                   % Impedance Spectrum Modeled from Peaks Info
            'wmodel', [] ...                      % Frequency axis referred to ReZmodel
        );
    
           
        % outputs created structre
        wfpstruct = globdata;
        
        % load structure in the workspace
        assignin('base','wfpstruct', globdata);
    
return;