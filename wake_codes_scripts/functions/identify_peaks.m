% This function was created to integrate the wake field processor pipeline
% and support data input for the coupled-bunch mode analysis
%
% It extracts peak information from the from the real part of the impedance
% spectra. The ReZ vs freq is interpolated and the findpeaks function is
% executed over it, generating information as center frequencis - omegar, shunt
% impedances and quality factors

function outstruct = identify_peaks(wfpstruct)

    addpath '\\CENTAURUS\Repositorio\Grupos\DIG\Projetos\Ativos\HC_SEM\_SEM\Scripts\_pipeline';

    % checks data type and perform factor 10 interpolation (needed for proper peak search and BW detection 
    if strcmp(wfpstruct.simpar.waketype,'long') && ~isempty(wfpstruct.results.ReZlong)
        wfpstruct.results.interReZlong = interp(wfpstruct.results.ReZlong,5);
        wake = wfpstruct.results.interReZlong;
    elseif strcmp(wfpstruct.simpar.waketype,'x') && ~isempty(wfpstruct.results.ReZx)
        wfpstruct.results.interReZx = interp(wfpstruct.results.ReZx,10);
        wake = wfpstruct.results.interReZx;
    elseif strcmp(wfpstruct.simpar.waketype,'y') && ~isempty(wfpstruct.results.ReZy)
        wfpstruct.results.interReZy = interp(wfpstruct.results.ReZy,10);
        wake = wfpstruct.results.interReZy;
    else
        wfpstruct.guipar.message = 'Invalid input data for peak identification. Canceled!';
        outstruct = wfpstruct;
        return
    end
    
    % perform factor 10 interpolation (needed for proper peak search and BW detection 
    wfpstruct.results.interfreq = interp(wfpstruct.results.freq,1);
      
   
    %  P=findpeaks(x,y,SlopeThreshold,AmpThreshold,SmoothWidth,FitWidth,smoothtype). Returns a list (in matrix P) containing the peak number and the
    %  estimated position, height, and width of each peak. Source: http://terpconnect.umd.edu/~toh/spectrum/PeakFindingandMeasurement.htm
    peakinfo = findpeaks(wfpstruct.results.interfreq,wake,0,2.0,5,5);
    
    % loads information of interest in the structure
    wfpstruct.results.peakinfo.omegar = 2*pi*1e9*peakinfo(:,2);
    wfpstruct.results.peakinfo.Rshunt = peakinfo(:,3);
    wfpstruct.results.peakinfo.Q = int32(peakinfo(:,2)./peakinfo(:,4));
    wfpstruct.results.peakinfo.BW = 2*pi*1e9*peakinfo(:,4);
    
    
    % informs the hard work is done
    wfpstruct.guipar.message = 'Peak identification sucessfully completed!';
    
    outstruct = wfpstruct;
    
    plot(wfpstruct.results.freq,wfpstruct.results.ReZlong,'r');
    hold on;
    plot(wfpstruct.results.interfreq,wfpstruct.results.interReZlong,'b');
    
    Fc = wfpstruct.results.peakinfo.omegar/(2*pi*1e9)
    Rsh = wfpstruct.results.peakinfo.Rshunt
    BW = wfpstruct.results.peakinfo.BW/(2*pi*1e6)

return