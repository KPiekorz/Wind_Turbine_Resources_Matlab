function missingdates = fcnmissingdates(time)
% fcnmissingdates.m
% 
% This function searches the time vector or gaps or sampling intervals that
% appear to inidicate that samples are missing.  
% 
% Usage: missingdates = fcnmissingdates(time)
% 
% Inputs:
% time = serial date vector of the recorded samples
% 
% Outputs:
% missingdates = matrix of the intervals where missing samples is detected

% Copyright 2009 - 2011 MathWorks, Inc.

% Find the sampling interval (using the median to resist outliers)
difft = diff(time);
dt = median(difft);      
% Find dates with a sampling interval larger then dt and store
I = find(difft > dt);
if isempty(I)
    disp('No log times are missing.')
else
    disp('Log times are missing. They are between (inclusive):')
    disp('      Start                 End')
    missingdates = zeros(length(I),2);
    for ii = 1:length(I)
        missingdates(ii,1) = time(I(ii)) + dt;
        missingdates(ii,2) = time(I(ii)+1) - dt;
        txt = [datestr(missingdates(ii,1),'yyyy-mm-dd HH:MM:SS') ...
                '  ' ...
                datestr(missingdates(ii,2),'yyyy-mm-dd HH:MM:SS')];
        disp(txt)
    end        
    disp(' ')
end