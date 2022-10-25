function [SlopeVector] = dABGdt(date,baselinemode, ABGTimes,pHvalues,pCO2values,pO2values,sO2values,time,PActData,DateTime)
%dABGdt generates the rate of change of different ABG values over time
%   date = the date of the experiment
%   baselinemode = 'avg' (reference the starting point as the avg baseline 
%       value for the animal) or 'prior' (reference the starting point
%       as the prior ABG to this set)
%   pHvalues = vector of pH values that correspond with ABGTimes
%   pCO2values = vector of pCO2 values that correspond with ABGTimes     
%   pO2values = vector of pO2 values that correspond with ABGTimes
%   sO2values = vector of sO2 values that correspond with ABGTimes
%   time = vector of timein s corresponding to powerlabdata
%   PActData = trimmed vector of data for the actuator waveform
%   DateTime = datetime of the start of the powerlab block


ABGfilename = [num2str(date) 'ABGData.mat'];
load(ABGfilename)
ABGslopefilename = [num2str(date) 'dABGdtData.mat'];
if isfile(ABGslopefilename)
   load(ABGslopefilename);
end

ABGandMinVentfilename = [num2str(date) 'ABGandMinVentData.mat'];

if isfile(ABGandMinVentfilename)
    load(ABGandMinVentfilename)
end

for i = 1:length(ABGTimes)
    abgIndex = find(modifiedABGDataCompiled{:, 1} == ABGTimes(i));     %%%%%%%%%%%%%%check what ABGtimes this is referencing, may be referencing wrong table
    
    if i==1% slope calculation needs to be determined by what is set as the 0 point

        [starttime] = findStartTime(time,PActData,DateTime,ABGTimes(1));
        dt = modifiedABGDataCompiled{abgIndex,1} - starttime;    
        switch baselinemode
            case 'avg'
                [baseline,SDbaseline] = calcBaseline(date);

            case 'prior'
                baseline = [ABGDataCompiled{abgIndex-1, 2} ABGDataCompiled{abgIndex-1, 3} ABGDataCompiled{abgIndex-1, 4} ABGDataCompiled{abgIndex-1, 5}]; 
        end
        dABG = [pHvalues(i) pCO2values(i) pO2values(i) sO2values(i)] - baseline;
    
    else
        dABG = [pHvalues(i) pCO2values(i) pO2values(i) sO2values(i)] - [pHvalues(i-1) pCO2values(i-1) pO2values(i-1) sO2values(i-1)];
        dt = modifiedABGDataCompiled{abgIndex,1}-modifiedABGDataCompiled{abgIndex-1,1};
 
    end
    

    slope = dABG./seconds(dt); 
        
    SlopeVector(i,:) = slope
    
    
    if exist('dABG_dtTable','var') == 0 %true if no table exists
        RateOfChangeavgbaseline = zeros(length(ABGDataCompiled{:, 1}),4);
        RateOfChangepriorbaseline = zeros(length(ABGDataCompiled{:, 1}),4);
        dABG_dtTable = addvars(modifiedABGDataCompiled,RateOfChangeavgbaseline,RateOfChangepriorbaseline);
    end

    switch baselinemode
        case 'avg'
            dABG_dtTable{abgIndex,7} = slope;

        case 'prior'
            dABG_dtTable{abgIndex,8} = slope;
    end
    

    save(ABGslopefilename,'dABG_dtTable')
    
end










end

