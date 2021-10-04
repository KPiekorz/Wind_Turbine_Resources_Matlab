function vbarlt = fcnGrabLongTermWind(WeatherLocation,beginyear,endyear)
%Calls function 'ScrapeDailyWeather' to scrape web site for averaged daily 
%   weather for all days in the range specified.  The years provided must
%   be full years, thus not the current year.  

% Copyright 2009 - 2011 MathWorks, Inc.

vm = zeros(12,endyear-beginyear+1);
years = beginyear:endyear;
% Loop through dates, scraping averaged monthly wind speeds
for iy = 1:length(years)
    for im = 1:12
        vm(im,iy)=ScrapeAvgMonthWeather(years(iy),im,WeatherLocation);
    end
end

% Compute long-term mean wind speed
vbarlt = nanmean(nanmean(vm));

end     % [EOF]

function vmonthly=ScrapeAvgMonthWeather(passyear,passmonth,WeatherLocation)
%Scrapes web site for daily weather and saves to dataset array
%   This function is used with the Weather demo.  It demonstrates data 
%   input by scrapping web sites with daily weather data.  
%   The url used is:
%       'http://www.wunderground.com/history/airport/KBOS/2009/1/1/
%                                             MonthlyHistory.html?format=1'
%   This function formats the url with the date and location passed as
%   arguments.

% Copyright The MathWorks, Inc. 2009

% Create and read the url
urlhead = 'http://www.wunderground.com/history/airport/';
urltail = '/MonthlyHistory.html?format=1';
url = [urlhead WeatherLocation '/' num2str(passyear) '/' ... 
       num2str(passmonth) '/1' urltail]; 

s = urlread(url);

% Remove newlines
expr = '\n';
s = regexprep(s, expr, '');

% Remove html tags
expr = '<(.*?)>';
s = regexprep(s, expr, ',');
 
% Remove last two commas
s = s(1:end-2);
 
% Convert to cell array
expr = ',';
tok = regexp(s, expr, 'split');
tok = reshape(tok,22,[])';

% Strip out wind speed
v = str2double(tok(2:end,18));
v = v*0.44704;      % 1 mph = 0.44704 meters/second

% Compute monthly mean wind speed (m/s)
vmonthly = nanmean(v);

end     % [EOF]