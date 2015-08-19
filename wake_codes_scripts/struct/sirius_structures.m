function globdata = sirius_structures()
%Loads sirius parameters
%
%   Type "edit sirius_parameters" to edit parameters.

globdata = struct( ...
    'ringpar' , [], ...                   % Storage Ring Parameters
    'guipar', [], ...                     % GUI parameters
    'simpar', [], ...                     % Simulation parameters
    'results', [] ...                     % Script Results
    );

globdata.ringpar = struct(...
    'circumf', 518.396, ...               % Ring circumference [m]
    'omega0', 3.63361516e6, ...           % Revolution frequency [rad/s]
    'sigma', 2.5e-3, ...                      % Nominal Bunch Length[m]
    'sigmamax', 15e-3, ...                   % Maximum Bunch Length of interest[m]
    'Iavg', 500e-3, ...                   % Beam current [A]
    'h', 864, ...                         % Harmonic number
    'beta', 0.9999999855, ...             % Beam speed (fraction of c)
    'E0', 3e9, ...                        % Beam Energy [eV]
    'Ee', 0.510998910e6, ...              % Electron Rest Energy [eV]
    'alpha', 1.74e-4, ...                 % Momentum Compaction Factor
    'nux', 48.0977, ...                   % Horizontal Tune
    'omegas', 2*pi*299792458/518.396*0.00393 ...             % Synchrotron Frequency [rad/s]
    );

globdata.simpar = struct( ...
    'wakepath' , [], ...                  % Path where the wake file (from any of the softwares) is
    'targetdir' , [], ...                 % Path where all results will be saved
    'datasource', [], ...                 % CST, ACE3P, ECHO or GdfidL
    'm', [],  ...                         % 0=long, 1=dipole trans, 2=quadrupole trans
    'offset', [], ...                     % offset for transverse analysis
    'sym', [],  ...                       % mirror symmetry for transverse?
    'whichaxis', [], ...                  % y or x
    'units', [], ...                       % # of components in the ring
    'bunlen', [] ...                     % Bunch Length Used in simulation[m]
    );

globdata.results = struct( ...
    's' , [], ...                         % axis: distance from following to drive bunch [m]
    'W' , [], ...                         % Longitudinal or Transverse Wakepotential [V/pC or V/pC/m]
    'freq', [], ...                       % axis: frequency obtained from FFT [GHz]
    'naxis', [], ...                      % axis: omega/omega0
    'interfreq', [], ...                  % interpolated frequency
    'ReZlong' , [], ...                   % Real Part of Longitudinal Impedance [Ohm]
    'interReZlong', [], ...               % interpolated impedance
    'ImZlong' , [], ...                   % Imaginary Part of Longitudinal Impedance [Ohm]
    'ImZoN' , [], ...                     % Imaginary Part of Longitudinal Impedance over n [Ohm]
    'interReZt', [], ...                  % interpolated impedance
    'ImZt' , [], ...                      % Imaginary Part of Vertical Dipole Impedance [KOhm/m]
    'ReZt' , [], ...                      % Real Part of Horizontal Dipole Impedance [KOhm/m]
    'interImZt', [], ...                  % interpolated impedance
    'peakinfo', [], ...'                  % Omegar [rad/s], Rshunt [Ohm] and Q from ReZ
    'klossW' , [], ...                    % Single-bunch Loss Factor Calculated from Wlong [mV/pC]
    'klossZ' , [], ...                    % Single-bunch Loss Factor Calculated from ReZlong [mV/pC]
    'kickW' , [], ...                     % Vertical Kick Factor Calculated from Wy [V/pC/m]
    'kickZ' , [], ...                     % Vertical Kick Factor Calculated from ImZy [V/pC/m]
    'sigmak' , [], ...                    % axis: bunch length for kloss|kick integration [mm]
    'Ploss', [], ...                      % Power loss from single-bunch loss factor [W]
    'Plossvec', [], ...                   % Power loss vector from single-bunch loss factor, for different sigmas [W]
    'klossWM' , [], ...                   % Multi-bunch Loss Factor Calculated from Wlong [mV/pC]
    'klossZM' , [], ...                   % Multi-bunch Loss Factor Calculated from ReZlong [mV/pC]
    'PlossM', [], ...                     % Power loss from multi-bunch loss factor [W]
    'ifast', [], ...                      % # of fastest CBM
    'GRs', [], ...                        % Growth Rate value for each CBM
    'GR_HOM', [], ...                     % Growth Rate value for each CBM accumulated through each HOM
    'ReZsampl', [], ...                   % Impedance Spectrum Sampled by fastest CBM
    'fsampl', [] ...                      % Frequency axis for sampled impedance
    );

globdata.results.peakinfo = struct( ...
    'omegar', [], ...                     % Modes Center Frequencies
    'BW', [], ...                         % Modes Bandwidths (half height)
    'Rshunt', [], ...                     % Modes Peak Values
    'Q', [], ...                          % Modes Quality Factors
    'ReZmodel', [], ...                   % Impedance Spectrum Modeled from Peaks Info
    'wmodel', [] ...                      % Frequency axis referred to ReZmodel
    );
