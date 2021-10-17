function fcnPlotDataQuality(data)

% Copyright 2009 - 2011 MathWorks, Inc.

   persistent figHandle

   if isempty(figHandle) || ~ishandle(figHandle)
      figHandle = figure;
   else
      figure(figHandle)    
   end

% Find begining of all months in date range
m = getMonth(data.xdata);
y = getYear(data.xdata);

dateGroups = unique([y m], 'rows');
dateGroups = sortrows(dateGroups,[1 2]);
dateGroups = dateGroups(1:2:length(dateGroups),:);
xt = datenum([dateGroups ones(size(dateGroups,1),1) ... 
    zeros(size(dateGroups,1),1) ones(size(dateGroups,1),1) ... 
    zeros(size(dateGroups,1),1)]);
   
   
   plot( data.xdata, toNaNandZero(data.goodData           )   , 'k.', ...
         data.xdata, toNaNandOne (data.missingData        )   , 'b.', ...
         data.xdata, toNaNandOne (data.clippedVelocity    )   , 'g.', ...
         data.xdata, toNaNandOne (data.clippedDirection   )   , 'y.', ...
         data.xdata, toNaNandOne (data.abnormalTemperature)   , 'r.', ...
         data.xdata, toNaNandOne (data.icingConditions    )   , 'c.', ...
         data.xdata, toNaNandOne (data.stuckWind          )   , 'm.'  );
   title('Data Quality Assurance Tests Results')
   %xlim([data.xdata(1) data.xdata(end)])
   ylim([-0.1 1.1])
   set(gca,'YTick',[0 1])
   set(gca,'YTickLabel',{'Passed','Failed'})
   set(gca,'XTick',xt)
   set(gca,'XTickLabel',cell(1,length(xt)))
   datetick('x','mmm-yy','keeplimits','keepticks')
   legend({'Good'     ,'Missing'     ,'Velocity', ...
           'Direction', 'Temperature', 'Icing','Stuck'}, ... 
           'Location', 'east');

end

%---------------------------------------------------------
% Find all the zeros and make them NaN so they won't plot
%---------------------------------------------------------
function outArray = toNaNandOne(inArray)
   outArray = double(inArray);
   outArray(outArray == 0) = NaN;
end

%---------------------------------------------------------
% Find all the zeros and make them NaN so they won't plot
% Find all the remaining ones and make them zeros
%---------------------------------------------------------
function outArray = toNaNandZero(inArray)
   outArray = double(inArray);
   outArray(outArray == 0) = NaN;
   outArray(outArray == 1) = 0  ;
end

