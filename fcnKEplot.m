function fcnKEplot(wind,ivh,wresults)
% fcnKEplot.m plot 
% 
% This function plots the velocity, and kinetic energy at the hub height
% time-series data. 
% 
% Usage: fcnKEplot(wind,ivh,wresults)
% 
% Inputs:
% wind = dataset array with all data
% ivh = index of velocity at hub height
% wresults = analysis output structure variable
% 
% Outputs: None

% Copyright 2009 - 2011 MathWorks, Inc.

%fcnKEfluxplot(wind.t, wind.vhub, wresults.overall.velocity(length(ivh)),...
%              wind.phub, wresults.overall.phub)
%fcnKEfluxplot(X1, Y1, Y1avg, Y3, Y3avg)

% Variable reassignments
t = wind.t;
vhub = wind.vhub;
vhubbar = wresults.overall.velocity(length(ivh));
phub = wind.phub;
phubbar = wresults.overall.phub;

clear wind wresults ivh

% Find begining of all months in date range
m = getMonth(t);
y = getYear(t);

dateGroups = unique([y m], 'rows');
dateGroups = sortrows(dateGroups,[1 2]);
dateGroups = dateGroups(1:2:length(dateGroups),:);
xt = datenum([dateGroups ones(size(dateGroups,1),1) ... 
    zeros(size(dateGroups,1),1) ones(size(dateGroups,1),1) ... 
    zeros(size(dateGroups,1),1)]);

% Create figure
sp1 = subplot(2,1,1);
    box(sp1,'on');
    hold(sp1,'all');
    plot(t,vhub,'Parent',sp1,'DisplayName','Data');
    plot([t(1) t(end)],[vhubbar vhubbar],'Parent',sp1, ... 
        'LineWidth',2,'Color',[1 0 0],'DisplayName','Mean');
    ylabel('v_{avg} (m/s)')
    set(gca,'XTick',xt)
    set(gca,'XTickLabel',cell(1,length(xt)))
    legend(sp1,'show','Location','NorthWest')
    hold(sp1,'off')
    
sp2 = subplot(2,1,2,'YScale','log');
    box(sp2,'on');
    hold(sp2,'all');
    semilogy(t,phub,'Parent',sp2,'DisplayName','Data');
    semilogy([t(1) t(end)],[phubbar phubbar],'Parent',sp2, ... 
        'LineWidth',2,'Color',[1 0 0],'DisplayName','Mean');
    ylabel('KE Flux (W/m^2)')
    set(gca,'XTick',xt)
    set(gca,'XTickLabel',cell(1,length(xt)))
    datetick('x','mmm-yy','keeplimits','keepticks')
    legend(sp2,'show','Location','SouthWest')
    hold(sp2,'off')
linkaxes([sp1, sp2], 'x');

% [EOF]