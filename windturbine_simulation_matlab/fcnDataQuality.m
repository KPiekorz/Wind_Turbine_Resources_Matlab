classdef fcnDataQuality < handle
    %FCNDATAQUALITY Stores the data quality information
    %   There are 6 conditions for poor data that are caputured 
    %   independently.
    %     missingData         - data is not present
    %     clippedVelocity     - velocity  measurements saturated
    %     clippedDirection    - direction measurements saturated
    %     abnormalTemperature - temperature readings unreasonable
    %     icingConditions     - icing detected
    %     stuckWind           - sensor stuck
    %     
    %   In addition, the following are computed automatically.
    %     
    %     badData             - at least one condition above is true
    %     goodData            - everything else
    
    % Copyright 2009 - 2011 MathWorks, Inc.
    
    % Basic fields of the data structure
    properties
        xdata
        missingData
        clippedVelocity
        clippedDirection
        abnormalTemperature
        icingConditions
        stuckWind
    end
    
    % Automatically computed from the other properties
    properties (Dependent = true)
        goodData %1=data is good, 0=data is bad
        badData  %1=data is bad , 1=data is good
    end

    methods

        % Constructor, called when first created
        function obj = fcnDataQuality(xdata)
           obj.xdata               = xdata;%store the x-values
           obj.missingData         = zeros(size(xdata));%1=fault, 0=nofault
           obj.clippedVelocity     = zeros(size(xdata));%1=fault, 0=nofault
           obj.clippedDirection    = zeros(size(xdata));%1=fault, 0=nofault
           obj.abnormalTemperature = zeros(size(xdata));%1=fault, 0=nofault
           obj.icingConditions     = zeros(size(xdata));%1=fault, 0=nofault
           obj.stuckWind           = zeros(size(xdata));%1=fault, 0=nofault
        end

        % Recompute nofaults every time goodData is used
        function value = get.goodData(obj)
           value = ~obj.badData;
        end

        % Recompute faults every time badData is used
        function value = get.badData(obj)
           value =  obj.missingData         | obj.clippedVelocity     | ...
                    obj.clippedDirection    | obj.abnormalTemperature | ...
                    obj.icingConditions     | obj.stuckWind   ;%find all faults
        end
    end
end
