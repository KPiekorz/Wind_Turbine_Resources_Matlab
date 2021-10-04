function weatherds = fcnGrabWeatherData(WeatherLocation,begindate,enddate)
%Calls function 'ScrapeDailyWeather' to scrape web site for daily weather
%   for all days in the range specified.  The dates should be in one of the
%   standard formats described in the help for the datestr function.

% Copyright 2009 - 2011 MathWorks, Inc.

% Convert date strings to date numbers
begindatenum = floor(datenum(begindate));
enddatenum = floor(datenum(enddate));

weatherds = dataset;

% Loop through dates, scraping daily data and adding to the dataset array
for idate = begindatenum:enddatenum
    ttok = ScrapeDailyWeather (idate, WeatherLocation);
    weatherds = [weatherds; ttok];        %#ok
end

end     % [EOF]

function w = ScrapeDailyWeather(passdate,WeatherLocation)
%Scrapes web site for daily weather and saves to dataset array
%   This function is used with the Weather demo.  It demonstrates data 
%   input by scrapping web sites with daily weather data.  
%   The url used is:
%       'http://www.wunderground.com/history/airport/KBOS/2009/1/1/
%                                               DailyHistory.html?format=1'
%   This function formats the url with the date and location passed as
%   arguments.

% Copyright The MathWorks, Inc. 2009
 
% Convert date from string to yyyy/mm/dd format
passdatestr = datestr(passdate, 26);

% Create and read the url
urlhead = 'http://www.wunderground.com/history/airport/';
urltail = '/DailyHistory.html?format=1';
url = [urlhead WeatherLocation '/' passdatestr urltail]; 

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
tok = reshape(tok,12,[])';

% Prep for dataset array
fname = genvarname(tok(1,:));
fname{1} = 'Datenum'; 
fname{8} = 'WindSpeedMS';

time = tok(2:end,1);

daten = cellfun(@(x)[datestr(passdate) ' ' x],time,'UniformOutput', false);
daten = datenum(daten);

% Convert to dataset array

w = dataset;
for jj = 1:size(tok,2)
    switch jj
        case {2,3,4,5,6,8,9,10}
            w.(fname{jj})= str2double(tok(2:end,jj));
            
        case {7,11,12}
            w.(fname{jj})= nominal(tok(2:end,jj));
        otherwise
            w.(fname{jj}) = daten;
    end
end

% Convert wind speed from MPH to m/s
w.(fname{8}) = w.(fname{8})*0.44704;      % 1 mph = 0.44704 meters/second
% Replace sensor errors (-9999 in original text from webiste) with NaNs
temp = w.(fname{8});
temp(temp < 0) = NaN;
w.(fname{8}) = temp;

end     % [EOF]