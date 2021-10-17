function wind = fcnImportDB

% Copyright 2009 - 2011 MathWorks, Inc.

% Set preferences with setdbprefs.
s.DataReturnFormat = 'structure';
s.ErrorHandling = 'store';
s.NullNumberRead = 'NaN';
s.NullNumberWrite = 'NaN';
s.NullStringRead = 'null';
s.NullStringWrite = 'null';
s.JDBCDataSourceFile = '';
s.UseRegistryForSources = 'yes';
s.TempDirForRegistryOutput = 'C:\Temp';
setdbprefs(s)

% Make connection to database.  Note that the password has been omitted.
% Using ODBC driver.
conn = database('winddata','','password');

% Read data from database.
e = exec(conn,'SELECT ALL date,v49Avg1,v49SD1,v49Max1,v49Min1,v49Avg2,v49SD2,v49Max2,v49Min2,v38Avg1,v38SD1,v38Max1,v38Min1,v38Avg2,v38SD2,v38Max2,v38Min2,v20Avg,v20SD,v20Max,v20Min,d49Avg,d49SD,d49Max,d49Min,d38Avg,d38SD,d38Max,d38Min,d20Avg,d20SD,d20Max,d20Min,T3Avg,T3SD,T3Max,T3Min FROM Winddata');
e = fetch(e);
close(e)

% Close database connection.
close(conn)

% Redefine output variable (added by TS)
wind = e.Data;