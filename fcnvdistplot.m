function fcnvdistplot(wresults,vnames)
% fcnvdistplot.m plot 
% 
% This function plots the velocity distributions.  
% 
% Usage: fcnvdistplot(wresults,vnames)
% 
% Inputs:
% wresutls = data structure with vdist data
% varnames = variable name for vdist
% 
% Outputs: None

% Copyright 2009 - 2011 MathWorks, Inc.

bar(wresults.vdist.vbins,wresults.vdist.(vnames))
axis tight
axis 'auto y'
xlabel('Wind velocity (m/s)')
ylabel('Fraction of Time')
title(['Wind Speed Distribution for ' vnames])

% Subplots are too small for most monitors. Use seperate figures for better
% images. 
% for ii = 1:length(vnames)
%     subplot(length(vnames),1,ii)
%         bar(wresults.vdist.vbins,wresults.vdist.(vnames{ii}))
%         axis tight
%         axis 'auto y'
%         xlabel('Wind velocity (m/s)')
%         ylabel('Fraction of Time')
%         title(['Wind Speed Distribution for ' vnames{ii}])
% end

% [EOF]