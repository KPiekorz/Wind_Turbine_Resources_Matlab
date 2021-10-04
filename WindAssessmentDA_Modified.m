
%% Import Data

wind = readtable('winddata_new.txt');

% Initialize storage of hub velocities, instantaneous kinetic energy flux, and air density in the dataset
wind.vhub = zeros(size(wind,1),1);
wind.phub = zeros(size(wind,1),1);
wind.rho = zeros(size(wind,1),1);
                            
% Convert date to a Serial Date Number (1 = January 1, 0000 A.D.)
wind.t = datenum(char(wind.date),' dd-mmm-yyyy HH:MM:SS');

% Input Additional Analysis Information

% Additional information about proposed wind turbine
hhub = 80;                % hub height (m)

% Additional information about meteorological tower
hv = [49 49 38 38 60];    % vector of heights for velocity (m)
hvh = [hv hhub];          % vector of heights for velocity with hhub (m)
hd = [49 38 60];          % vector of heights for direction (m)
hT = 2;                   % vector of heights for temperature (m) 

nobs = size(wind,1);      % # of observations

% Air properties
patm = 101e3;             % atmospheric pressure (Pa)
Rair = 287;               % gas constant for air (J/kg K)

% air density (kg/m^3)
wind.rho = patm./(Rair*(wind.T3Avg+273.15));

clear patm Rair

% Sensor indices
ihub = 38;                % index of hub velocity estimates
iv = 2:4:18;              % indices of velocity measurements
ivh = [iv ihub];          % indices of velocity (with hub)
id = 22:4:30;             % indices of direction measurements
iT = 34;                  % indices of temperature measurements

% Expected ranges for sensors
vrange  = [0 100];        % min and max velocity (m/s)
drange = [0 360];         % expected range for direction (deg)
Trange = [-50 150];       % expected range for temperature measurements (C)

% Critical values for Icing Test
% vAvg > vice & dSD <= dstdice & TAvg < Tice
vice = 1;              % critical value for wind speed (m/s)
dstdice = 0.5;         % critical value for the std of wind direction (deg)
Tice = 2;              % critical value for temperature (C)

% indices of sensor sets for icing tests [vAvg dSD TAvg]
iice = [ 2 23 34;
         6 23 34;
        10 27 34;
        14 27 34;
        18 31 34];
    
% Critical values for stuck wind direction sensor
% dSD<dSDstuck & diff(d)<ddelta for at least ndt consecutive time samples
dSDstuck = 0.1;        % critical value for wind direction std (deg)
ddelta = 0.1;          % critical value for wind direction difference (deg)
ndt = 6;               % min # of time sampes for stuck conidition
istuck = [22 26 30];   % indices for wind direction sensor stuck test

% Create variable to store statistical analysis results
wresults = [];       % structure variable for results

% Visualize Data
clc
close all
figure
fcnvdttimeplot(wind)

%% Hub Height Wind Velocity Estimate

% Estimate the wind velocity at hub height using a power law model fitted
% to the measured wind velocities for each time sample.  

clc
npass = size(wind,1);
vhub = zeros(npass,1);

parfor ii = 1:100

    % Compute instantaneous power law shear models
    cfobj = fcnpowerlaw(hv, table2array(wind(ii, iv)))

    % Compute estimate of wind speed at the wind turbine hub height
    vhub(ii) = cfobj(hhub);
end

wind.vhub = vhub;

clear cfobj vhub

%% Compute Overall Averages

% Store overall averages (include hub height velocity with the velocity
% data)
wresults.overall.velocity = mean(table2array(wind(:,ivh)));
wresults.overall.direction = mean(table2array(wind(:,id)));
wresults.overall.temperature = mean(table2array(wind(:,iT))); 

% fprintf("Overall averages, velocity: %f, direction: %f, temperature: %f", ...
%         wresults.overall.velocity, wresults.overall.direction, wresults.overall.temperature);        

%% Wind Speed Distribution

close all

% Another view on the data is to compute and display the frequency the
% averaged wind speed was with in a certain range.  Let's create the wind
% speed distribution. 
vmax = max(max(table2array(wind(:,ivh))));

% Bin centers for histogram (m/s)
wresults.vdist.vbins = (0:1:ceil(vmax))';

% vnames = wind.Properties.VarNames(ivh);
vnames = ["v49Avg1" "v49Avg2" "v38Avg1" "v38Avg2" "v60Avg" "vhub"];

% Compute distributions for all averaged velocity data columns included the
% estimate at the hub height
for ii = 1:length(vnames)

    wresults.vdist.(vnames{ii}) = hist(wind.(vnames{ii}), wresults.vdist.vbins) / npass;

    figure(ii);
    fcnvdistplot(wresults, vnames{ii});
end

clear ii vmax vnames
    
%% Wind Rose

close all

% Create the wind rose plots where the direction represents the direction 
% the wind is blowing from. 

% Use wind_rose function from MATLAB Central with a small modification 
% regarding the meteorological angle conversion.
% (http://www.mathworks.com/matlabcentral/fileexchange/17748)
figure('color','white')
    fcnwindrose(wind.d49Avg, wind.v49Avg1,'dtype','meteo','n',16, ... 
         'labtitle','Height = 49 m, Sensor 1','lablegend','Velocity (m/s)')
figure('color','white')
    fcnwindrose(wind.d49Avg, wind.v49Avg2,'dtype','meteo','n',16, ... 
         'labtitle','Height = 49 m, Sensor 2','lablegend','Velocity (m/s)')
figure('color','white')
    fcnwindrose(wind.d38Avg, wind.v38Avg1,'dtype','meteo','n',16, ... 
         'labtitle','Height = 38 m, Sensor 1','lablegend','Velocity (m/s)')
figure('color','white')
    fcnwindrose(wind.d38Avg, wind.v38Avg2,'dtype','meteo','n',16, ... 
         'labtitle','Height = 38 m, Sensor 2','lablegend','Velocity (m/s)')
figure('color','white')
    fcnwindrose(wind.d60Avg, wind.v60Avg,'dtype','meteo','n',16, ... 
         'labtitle','Height = 60 m','lablegend','Velocity (m/s)')

%% Turbulence Intensity

% Compute the turbulence intensity for each observation and velocity 
% sensor and the distribution for each sensor.  The turbulence intensity is
% defined as the 10-minute standard deviation of the velocity divided by 
% the 10-minute average velocity. 

% Compute turbulence intensities
wresults.ti.data = table2array(wind(:,iv+1))./table2array(wind(:,iv));

% Display turbulence intensities versus wind speed for each sensor
timax = ceil(10*max(max(wresults.ti.data)))/10;
vmax = ceil(max(max(table2array(wind(:,iv)))));

% Visualize data
for ii = 1:length(iv)
    figure
        subplot(2,1,1);
            plot(table2array(wind(:,iv(ii))),wresults.ti.data(:,ii),'+')
            xlim([0 ceil(vmax)])
            ylim([0 ceil(10*timax)/10])
            box on
            %xlabel('Wind velocity (m/s)')
            ylabel('TI')
            title(['Turbulence Intensity for ' ...
                char(wind.Properties.VarNames(iv(ii)))])
        subplot(2,1,2);
            boxplot(wresults.ti.data(:,ii),round(table2array(wind(:,iv(ii)))))
            xlim([0 ceil(vmax)])
            xlabel('Wind velocity (m/s)')
            ylabel('TI')
end

clear ii timax vmax
        
%% Shear Profile

% Compute the shear exponent of the power law model for the atmospheric 
% boundary layer.  The power law model is of the form u = a*z^alpha.  The
% coefficient, a, and the exponent, alpha, are estimated using regression
% analysis using the data from all velocity sensors.  Remember to exclude
% the estimated velocity at hub height in the analysis of the shear or
% boundary layer profile.  

% Fit for overall data
[cfobj,cfgood] = fcnpowerlaw(hv,wresults.overall.velocity(1:length(hv)));
wresults.bl.overall.cfobj = cfobj;
wresults.bl.overall.cfgood = cfgood;
wresults.bl.overall.alpha = coeffvalues(cfobj);
wresults.bl.overall.alpha = wresults.bl.overall.alpha(2);

% Plot of overall data and fit (up to a height of 100 m)
x = logspace(-2,2,200);
y = cfobj(x);

figure
    plot(wresults.overall.velocity(1:length(hv)),hv,'o',y,x)
    xlabel('Wind velocity (m/s)')
    ylabel('Height (m)')
    legend('Data','Power Law','Location','Best')

clear x y cfobj cfgood

% % Fit for monthly average wind speed
% wresults.bl.monthly.date = wresults.monthavg.date;
% nmonths = length(wresults.bl.monthly.date);
% wresults.bl.monthly.cfobj = cell(length(nmonths),1);
% wresults.bl.monthly.cfgood = cell(length(nmonths),1);
% wresults.bl.monthly.alpha = zeros(length(nmonths),1);
% 
% for ii = 1:nmonths
%     [cfobj,cfgood] = fcnpowerlaw(hv, ... 
%                                   wresults.monthavg.data(ii,1:length(hv)));
%     alpha = coeffvalues(cfobj);
%     wresults.bl.monthly.cfobj{ii} = cfobj;
%     wresults.bl.monthly.cfgood{ii} = cfgood;
%     wresults.bl.monthly.alpha(ii) = alpha(2); 
% end
% 
% clear ii cfobj cfgood alpha
% 
% % Plot monthly and overal alpha values
% fcnalphaplot(1:nmonths,wresults.bl.monthly.alpha,wresults.bl.overall.alpha)
%     
% clear nmonths

%% Wind Power and Capacity Factor Estimate

% Compute the mean wind speed, kinetic energy flux, and the capacity factor
% using both the local, short term data and correlation to a data source
% with long term data accessible.  This report correlates the local site at
% Cohasset to the weather station at Boston Logan Internation Airport
% (KBOS).  

%% Short Term Kinetic Energy Flux

% Compute the power density in the wind as measured by the meteorological
% tower.  This will represent the local, short-term power that was 
% available to any wind turbine installed at this location during the
% measurement period.  

% Compute instantaneous kinetic energy flux
wind.phub = 0.5*wind.rho.*wind.vhub.^3;

% Compute and store mean wind kinetic energy flux at hub height
wresults.overall.phub = mean(wind.phub);
% display to published report
disp(['Mean KE flux (W/m^2): ' num2str(wresults.overall.phub,'%3.0f')])
disp(' ')

% Plot instantaneous hub wind speeds and KE flux
figure
fcnKEplot(wind,ivh,wresults)

%% Short Term Average Turbine Power and Capacity Factor

% Estimate the average turbine power and capacity factor for this site 
% using the short-term estimated hub height velocity distribution.  These
% calculations require knowledge of the proposed wind turbine model and its
% power curve.  For this demo, let's assume a 1.5 MW wind turbine with a
% power curve modelled in fcnpowercurve.

% Wind turbine rated power (W)
Prated = 1500e3;
wresults.short.Prated = Prated;

% Compute Pavgshort as the integral of the wind turbine power curve and the
% pdf of the wind speed at the hub height.  
dx = mean(diff(wresults.vdist.vbins));      % integral steps (m/s)
Pavgshort = sum(fcnpowercurve(wresults.vdist.vbins,Prated).* ... 
                wresults.vdist.vhub(:))*dx;
wresults.short.Pavgshort = Pavgshort;
            
% Compute short-term capacity factor
CFshort = Pavgshort/Prated;
wresults.short.CFshort = CFshort;

% Display results
disp(['Assumed wind turbine rated power (MW): ' num2str(Prated/1e6,'%3.1f')])
disp(['Short-term averaged power (kW): ' num2str(Pavgshort/1e3,'%3.0f')])
disp(['Short-term Capacity Factor (%): ' num2str(CFshort*100,'%2.0f')])
disp(' ')

clear Prated Pavgshort dx CFshort
