function analyzeMinVentAndABGMODIFIED(date,TminVent,PminVent,NminVent, T,ABGTimes,pHvalues,pCO2values,pO2values,sO2values, baseline, SDbaseline,time, PActData,DateTime,plotsON)
%analyzeMinVentAndABG This function pulls data for a period of T prior to
%each ABG taken and fine the average minute ventilation. It will then
%process the ABG data as related to this average minute ventilation and
%generate a correlation graph. 
%   date = the date of the experiment
%   TminVent = a time vector for the minute vent calculation with each T
%   at end of each dt bucket (calculated from minVentOverTime)
%   PminVent = a vector of the minute ventilation calculated off of
%   inspiration over timesteps of dt (calculated from minVentOverTime)
%   NminVent = a vector of the minute ventilation calculated off of
%   exspiration over timesteps of dt (calculated from minVentOverTime)
%   T = the period of time (in seconds) over which to calculate an average minute
%   ventilation
%   ABGTimes = a vector of times that the selected ABGs were take 
%   pHvalues = vector of pH values that correspond with ABGTimes
%   pCO2values = vector of pCO2 values that correspond with ABGTimes     
%   pO2values = vector of pO2 values that correspond with ABGTimes
%   sO2values = vector of sO2 values that correspond with ABGTimes
%   baseline = vector of average baseline values [pH pCO2 pO2 sO2] for
%       given date
%   SDbaseline = vector of stdev of baseline values [pH pCO2 pO2 sO2]
%   time = vector of time in s corresponding to the powerlab data
%   PActData = trimmed actuator data from powerlab
%   DateTime = datetime of the start of the data from the powerlab
%   plotON = whether or not to plot the min vent x ABG values


%load ABG data that was imported from excel file
ABGfilename = [num2str(date) 'ABGData.mat'];
load(ABGfilename)




%% Average minute ventilation over period T

ABGandMinVentfilename = [num2str(date) 'ABGandMinVentData.mat'];

if isfile(ABGandMinVentfilename)
    load(ABGandMinVentfilename)
end

minventTname = sprintf('Avg Min Vent with Period T = %d[s]',T);
MinVentVector = [];
for i = 1:length(ABGTimes)
    abgIndex = find(modifiedABGDataCompiled{:, 1} == ABGTimes(i)); 
    timeindlogical = find(TminVent >= ABGTimes(i)); 
    closestInd = timeindlogical(1);
    PeriodStartIndLogical = find(TminVent >= ABGTimes(i)-seconds(T)); 
    PeriodStartInd = PeriodStartIndLogical(1);
    indices = (PeriodStartInd:closestInd);
    AvgMinVent = mean(PminVent(indices));
    
    

    if exist('ABGMinVentTable','var') == 0 %true if no table exists
        MinVent = zeros(length(ABGDataCompiled{:, 1}),1);
        ABGMinVentTable = addvars(modifiedABGDataCompiled,MinVent);
        ABGMinVentTable.Properties.VariableNames{end} = minventTname;
       %the last column of ABGMinVentTable will be the appropriate location to
       %place min vent data for this T value

        col = width(ABGMinVentTable);
    elseif ~any(strcmp(ABGMinVentTable.Properties.VariableNames,minventTname)) %true if no columns are titled with correct period
        MinVent = zeros(length(ABGDataCompiled{:, 1}),1);
        ABGMinVentTable = addvars(ABGMinVentTable,MinVent);
        ABGMinVentTable.Properties.VariableNames{end} = minventTname;
       %the last column of ABGMinVentTable will be the appropriate location to
       %place min vent data for this T value
        col = width(ABGMinVentTable);

    elseif any(strcmp(ABGMinVentTable.Properties.VariableNames,minventTname)) %the header already exists
        col = find(strcmp(ABGMinVentTable.Properties.VariableNames,minventTname));

    end
    ABGMinVentTable{abgIndex,col} = AvgMinVent;
    MinVentVector(i) = AvgMinVent;


%% Normalized to avg baseline

if exist('ABGAdjBaseline','var') == 0 %true if this variable does not exist
    rel_pH = zeros(length(ABGDataCompiled{:, 1}),1);
    rel_pCO2 = zeros(length(ABGDataCompiled{:, 1}),1);
    rel_pO2 = zeros(length(ABGDataCompiled{:, 1}),1);
    rel_sO2 = zeros(length(ABGDataCompiled{:, 1}),1);
    ABGAdjBaseline = addvars(modifiedABGDataCompiled,rel_pH,rel_pCO2,rel_pO2,rel_sO2);
end
ABGAdjBaseline{abgIndex,7} = pHvalues(i)-baseline(1);
ABGAdjBaseline{abgIndex,8} = pCO2values(i)-baseline(2);
ABGAdjBaseline{abgIndex,9} = pO2values(i)-baseline(3);
ABGAdjBaseline{abgIndex,10} = sO2values(i)-baseline(4);

%% Slope Analysis

end

%% 





save(ABGandMinVentfilename, 'ABGMinVentTable','ABGAdjBaseline','')


%% Process ABG Data
%this is redundant but useful for plotting
relpH = pHvalues-baseline(1);
relpCO2 = pCO2values-baseline(2);
relpO2 = pO2values-baseline(3);
relsO2 = sO2values-baseline(4);




%% Slope Analysis
% baselinemode = 'avg';
% [SlopeVector] = dABGdt(date,baselinemode, ABGTimes,pHvalues,pCO2values,pO2values,sO2values,time,PActData,DateTime)


%% Plot ABG Data

if plotsON == 1
    figure
    title(sprintf('Minute ventilation over a period T = %d[s] vs absolute ABG values',T))
    subplot (2,2,1)
    plot(MinVentVector,pHvalues,'x')
    xlabel('Average Minute Ventilation [L]')
    ylabel('pH')
    
    subplot (2,2,2)
    plot(MinVentVector,pCO2values,'x')
    xlabel('Average Minute Ventilation [L]')
    ylabel('pCO2 [mmHg]')
    
    subplot (2,2,3)
    plot(MinVentVector,pO2values,'x')
    xlabel('Average Minute Ventilation [L]')
    ylabel('pO2 [mmHg]')
    
    subplot (2,2,4)
    plot(MinVentVector,sO2values,'x')
    xlabel('Average Minute Ventilation [L]')
    ylabel('sO2 [%]')
    
    figure
    title(sprintf('Minute ventilation over a period T = %d[s] vs ABG values relative to baseline',T))
    subplot (2,2,1)
    plot(MinVentVector,relpH,'x')
    xlabel('Average Minute Ventilation [L]')
    ylabel('pH Rate of Change [/s]')
    
    subplot (2,2,2)
    plot(MinVentVector,relpCO2,'x')
    xlabel('Average Minute Ventilation [L]')
    ylabel('pCO2 Rate of Change [mmHg/s]')
    
    subplot (2,2,3)
    plot(MinVentVector,relpO2,'x')
    xlabel('Average Minute Ventilation [L]')
    ylabel('pO2 Rate of Change [mmHg/s]')
    
    subplot (2,2,4)
    plot(MinVentVector,relsO2,'x')
    xlabel('Average Minute Ventilation [L]')
    ylabel('sO2 Rate of Change [%/s]')
    
%     figure
%     subplot (2,2,1)
%     plot(MinVentVector,SlopeVector(:,1),'x')
%     xlabel('Average Minute Ventilation [L]')
%     ylabel('\Delta pH')
%     
%     subplot (2,2,2)
%     plot(MinVentVector,SlopeVector(:,2),'x')
%     xlabel('Average Minute Ventilation [L]')
%     ylabel('\Delta pCO2 [mmHg]')
%     
%     subplot (2,2,3)
%     plot(MinVentVector,SlopeVector(:,2),'x')
%     xlabel('Average Minute Ventilation [L]')
%     ylabel('\Delta pO2 [mmHg]')
%     
%     subplot (2,2,4)
%     plot(MinVentVector,SlopeVector(:,4),'x')
%     xlabel('Average Minute Ventilation [L]')
%     ylabel('\Delta sO2 [%]')
    
end







end

