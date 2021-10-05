function p = fcnpowercurve(v,prated)
% fcnpowercurve.m 
% 
% This function models the power curve for a fictitious wind turbine.  The
% model is constructed as a piece-wise function to more accurately account
% for the cut-on and cut-out speeds.  
% 
% Usage: p = fcnpowercurve(v,prated)
% 
% Inputs:
% v = wind speed at hub height (m/s)
% prated = rated turbine power (W)
% 
% Outputs:
% p = electrical power from wind turbine (W)

% Copyright 2009 - 2011 MathWorks, Inc.
%   Author(s): T. Schultz, 6/23/2009 4.4947*(v(i)^3) - *(v(i)^2) + 260.034*v(i) - 624.5262;

% Constants
von = 3;        % cut-on speed (m/s)
vc = 11.4;        % corner speed (m/s)
vout = 25;      % cut-out speed (m/s)

p = zeros(size(v));

% Create model
% Below cut-on
p(v < von) = 0;
% Ramp up (use model)
I = (v >= von & v < vc);
coefficients = [4.4947 30.7702 260.034 624.5262];
p(I) = polyval(coefficients, v(I));
% At rated power
p(v >= vc & v <= vout) = prated;
% Above cut-out
p(v > vout) = 0;

% [EOF]