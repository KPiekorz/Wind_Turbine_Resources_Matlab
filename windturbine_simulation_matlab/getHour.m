function hour   = getHour(dates)

% Copyright 2009 - 2011 MathWorks, Inc.

   parsedData   = datevec(dates);
   hour         = parsedData(:,4);
end
