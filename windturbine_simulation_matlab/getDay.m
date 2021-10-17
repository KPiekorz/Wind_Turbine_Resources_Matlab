function day    = getDay(dates)

% Copyright 2009 - 2011 MathWorks, Inc.

   parsedData   = datevec(dates);
   day          = parsedData(:,3);
end
