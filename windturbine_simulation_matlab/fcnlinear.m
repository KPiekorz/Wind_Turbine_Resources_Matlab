function [cfobj,cfgood] = fcnlinear(x,y)
% fcnpowerlaw.m 
% 
% This function fits a linear model, y = a*x+b, to the data provide in
% the variables x and y.  The outputs are the curve fit object and the
% goodness of fit structure.  Note that inputs x and y must be the same
% size vectors. 
% 
% Usage: [fobj,fgood] = fcnlinear(x,y)
% 
% Inputs:
% x = independent or predictor variable data as a vector
% y = dependent or response variable data as a vector
% 
% Outputs:
% cfobj = curve fit object with fit results
% cfgood = structure with goodness of fit information

% Copyright 2009 - 2011 MathWorks, Inc.
%   Author(s): T. Schultz, 7/2/2009

% Ensure column vectors
x = x(:);
y = y(:);

% Look for not finite numbers and ignore
Iok = isfinite(x) & isfinite(y);

% Look for outliers and ignore
thres = 30;             % critical value for outliers
Iout = x <= thres & y <= thres;

% Combine tests
I = Iok & Iout;

% Calculate fit
[cfobj,cfgood] = fit(x(I),y(I),'poly1');

% [EOF]