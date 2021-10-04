function month  = getMonth(dates)

% Copyright 2009 - 2011 MathWorks, Inc.

   parsedData   = datevec(dates);
   month        = parsedData(:,2);
end
