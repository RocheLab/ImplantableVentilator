function [ABGTimes,pHvalues,pCO2values,pO2values,sO2values,PriorPoint] = callABGData(date,DateTime, timeStart, timeDur)
%callABGData is modified fromalignAndPlotABGData, it is used to read
%already aligned data 
%   date =  date of the experiment
%   DateTime = is the datetime format of the absolute time that the block
%   of PowerLab data began recording
%   timeStart = number of seconds after DateTime that marks the beginning
%   of the selected block of data to analyze
%   timeDur = nlength of this selection of data in umber of seconds 
%   time = a vector of time series data from the arterial pressure
%   PArtData = data from the arterial pressure from the specific block we
%   are looking at, this will be an aligning source (if there is no
%   corresponding alignment piece of data, we will use the excel sheet
%   printing time as the reference unless otherwise specified. 

%load ABG data that was imported from excel file
ABGfilename = [num2str(date) 'ABGData.mat'];
load(ABGfilename)

%find span of ABG data points that fall within the span of this selection
selectionStart = DateTime + seconds(timeStart);
selectionEnd = selectionStart+seconds(timeDur);
%find the first ABG in range
logicalStart = find(modifiedABGDataCompiled{:, 1} >= (selectionStart)); 
%find the last ABG in range
logicalEnd = find(modifiedABGDataCompiled{:, 1} <= (selectionEnd));
if isempty(logicalStart) || isempty(logicalEnd) %no paper ABGs in this time range
    disp('No ABGs taken for this block.')
    ABGTimes = [];
    pHvalues = [];
    pCO2values = [];
    pO2values = [];
    sO2values = [];
    return
end
startInd = logicalStart(1);
endInd = logicalEnd(end); 
if isempty(startInd) && isempty(endInd) %this means no ABGs are found for this selection
    disp('No ABGs found for this selection')
    return
elseif isempty(startInd) && ~isempty(endInd)
    endInd = length(ABGDataCompiled{:,1});
elseif ~isempty(startInd) && isempty(endInd)
    startInd = 1;
end

%ABG data in range
ABGTimes = modifiedABGDataCompiled{startInd:endInd, 1};%column 1 has slightly better alignment with powerlab times
pHvalues = ABGDataCompiled{startInd:endInd, 3};
pCO2values = ABGDataCompiled{startInd:endInd, 4};
pO2values = ABGDataCompiled{startInd:endInd, 5};
sO2values = ABGDataCompiled{startInd:endInd, 6};

PriorPoint = [ABGDataCompiled{startInd-1,3} ABGDataCompiled{startInd-1,4} ABGDataCompiled{startInd-1,5} ABGDataCompiled{startInd-1,6}];
end

