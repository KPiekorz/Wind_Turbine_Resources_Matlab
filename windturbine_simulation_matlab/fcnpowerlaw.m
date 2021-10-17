function [cfobj,cfgood] = fcnpowerlaw(x,y)
% fcnpowerlaw.m 
% 
% This function fits a power law model, y = a*x^b, to the data provide in
% the variables x and y.  The outputs are the curve fit object and the
% goodness of fit structure.  Note that inputs x and y must be the same
% size vectors. 
% 
% Usage: [fobj,fgood] = fcnpowerlaw(x,y)
% 
% Inputs:
% x = independent or predictor variable data as a vector
% y = dependent or response variable data as a vector
% 
% Outputs:
% cfobj = curve fit object with fit results
% cfgood = structure with goodness of fit information

% Copyright 2009 - 2011 MathWorks, Inc.
%   Author(s): T. Schultz, 6/10/2009

% Ensure column vectors
x = x(:);
y = y(:);

% Calculate fit
[cfobj,cfgood] = fit(x,y,'power1');

% [EOF]