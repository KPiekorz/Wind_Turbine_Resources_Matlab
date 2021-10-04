function year   = getYear(dates)

% Copyright 2009 - 2011 MathWorks, Inc.

   parsedData   = datevec(dates);
   year         = parsedData(:,1);
end
