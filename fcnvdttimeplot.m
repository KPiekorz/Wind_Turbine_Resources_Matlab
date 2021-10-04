function fcnvdttimeplot(wind)
% fcnvdttimeplot.m plot 
% 
% This function plots the velocity, direction, and temperature time-series
% data. 
% 
% Usage: fcnvdttimeplot(wind)
% 
% Inputs:
% wind = dataset array with all data
% 
% Outputs: None

% Copyright 2009 - 2011 MathWorks, Inc.

% Find begining of all months in date range
m = getMonth(wind.t);
y = getYear(wind.t);

dateGroups = unique([y m], 'rows');
dateGroups = sortrows(dateGroups,[1 2]);
dateGroups = dateGroups(1:2:length(dateGroups),:);
xt = datenum([dateGroups ones(size(dateGroups,1),1) ... 
    zeros(size(dateGroups,1),1) ones(size(dateGroups,1),1) ... 
    zeros(size(dateGroups,1),1)]);

% Create figure
sp1 = subplot(3,1,1);
    plot(wind.t,wind.v49Avg1,wind.t,wind.v49Avg2, ...
         wind.t,wind.v38Avg1,wind.t,wind.v38Avg2,wind.t,wind.v20Avg)
    legend('v49a','v49b','v38a','v38b','v20')
    ylabel('v_{avg} (m/s)')
    set(gca,'XTick',xt)
    set(gca,'XTickLabel',cell(1,length(xt)))
sp2 = subplot(3,1,2);
    plot(wind.t,wind.d49Avg,wind.t,wind.d38Avg,wind.t,wind.d20Avg)
    legend('d49','d38','d20')
    ylabel('d_{avg} (deg)')
    set(gca,'XTick',xt)
    set(gca,'XTickLabel',cell(1,length(xt)))
sp3 = subplot(3,1,3);
    plot(wind.t,wind.T3Avg)
    legend('T3')
    ylabel('T_{avg} (\circC)')
    set(gca,'XTick',xt)
    datetick('x','mmm-yy','keeplimits','keepticks')
linkaxes([sp1, sp2, sp3], 'x');

% [EOF]