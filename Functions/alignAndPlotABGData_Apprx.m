function [ABGTimes,pHvalues,pCO2values,pO2values,sO2values] = alignAndPlotABGData_Apprx(date,DateTime, timeStart, timeDur, time,PArtData)
%alignAndPlotABGData takes information from one section of data, and aligns
%the ABGs with the appropriate timestamp of the data, and then blots the
%ABG data as a scatter plot over time relative to that section of data. 
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

load(ABGfilename)

%find span of ABG data points that fall within the span of this selection
selectionStart = DateTime + seconds(timeStart);
selectionEnd = selectionStart+seconds(timeDur);
%find the first ABG in range
logicalStart = find(ABGDataCompiled{:, 8} >= (selectionStart-minutes(1))); %gives a 1 minute buffer just in case
%find the last ABG in range
logicalEnd = find(ABGDataCompiled{:, 8} <= (selectionEnd+minutes(1))); %gives a 1 minute buffer just in case
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
indices = [startInd:endInd];
if isempty(startInd) && isempty(endInd) %this means no ABGs are found for this selection
    disp('No ABGs found for this selection')
    return
elseif isempty(startInd) && ~isempty(endInd)
    endInd = length(ABGDataCompiled{:,8});
elseif ~isempty(startInd) && isempty(endInd)
    startInd = 1;
end

%ABG data in range
ABGTimes = ABGDataCompiled{startInd:endInd, 8}%column 8 has approx times that were manually derived
pHvalues = ABGDataCompiled{startInd:endInd, 3};
pCO2values = ABGDataCompiled{startInd:endInd, 4};
pO2values = ABGDataCompiled{startInd:endInd, 5};
sO2values = ABGDataCompiled{startInd:endInd, 6};

AbsTime = DateTime + seconds(time);



%% Adjust ABG values to match with the PArt Data
%Find where the ABG is shown in the PArt Data
SatInd = find(PArtData>1.92); %threshold indices to 1.92
GapInd = find((SatInd(2:end)-SatInd(1:end-1))>60000);
for i=1:(length(GapInd)+1)%generates avgABGtime
    cluster = [];
    if (i == 1) && (isempty(GapInd)==0) %first cluster
        cluster = SatInd(1:GapInd(1));
    elseif (i == 1) && (isempty(GapInd) == 1) %only 1 ABG
        cluster = SatInd;
    elseif i == (length(GapInd)+1) %last cluster
        cluster = SatInd(GapInd(end)+1:end);
    else
        cluster = SatInd(GapInd(i-1)+1:GapInd(i));
    end
    
    avgABGtime(i)= mean(cluster) %average time of clusters in seconds
end 

%avgABGtime is a vector of the time (in ms) of the average location of each
%cluster of PArt data indicating an ABG was taken, length(avgABGtime) is
%equal to the number of ABGs taken in this data set
PArtABGTimes = selectionStart + seconds(0.001.*avgABGtime)
    
    NPArt = length(avgABGtime) ;
    NPaper = length(ABGTimes);
%Align new ABG times
%Adjust number of NPArt and NPaper if there is imsmatch, if there is no
%mismatch proceed forward without this "if" block


if NPArt > NPaper
    disp('A mismatch is found within the number of ABGs in this section.')
    fprintf('Within this time range +/-3 minutes, there are %d values found from the paper times and %d values found from the PArt waveform.',NPaper,NPArt)

    figure
    plot(AbsTime, PArtData,'-',PArtABGTimes,ones(1,NPArt),'xr') 

    prompt = 'Would you like to ignore 1 (or more) of the PArt ABG time indices? [Y/N] (Y = eg. if there was an injection of medication, it can look like an ABG waveform, N = eg. if the timespan is not wide enough to catch all the paper reads)';
    response = input(prompt,'s');
        if response == 'Y'
            prompt2 = 'Which ABG time index would you like to ignore? (input a number)';
            response2 = input(prompt2);
            PArtABGTimes(response2) = [];
        elseif response == 'N'
            prompt2 = 'Do you want to add in the prior paper read?(paper times are often ahead of PArt times) [Y/N]';
            response2 = input(prompt2,'s');
            if response2 == 'Y'
                indices = [startInd-1:endInd];
                ABGTimes = ABGDataCompiled{indices, 8};%column 1 has slightly better alignment with powerlab times
                pHvalues = ABGDataCompiled{indices, 3};
                pCO2values = ABGDataCompiled{indices, 4};
                pO2values = ABGDataCompiled{indices, 5};
                sO2values = ABGDataCompiled{indices, 6};

            end
        end

elseif NPArt < NPaper
    disp('A mismatch is found within the number of ABGs in this section.')
    fprintf('Within this time range +/-1 minutes, there are %d values found from the paper times and %d values found from the PArt waveform.\n',NPaper,NPArt)
    
    figure
    plot(AbsTime, PArtData,'-',ABGTimes,ones*(NPaper),'xr') 

    prompt = 'Would you like to ignore 1 (or more) of the Paper ABG time indices? [Y/N] (Y=eg. if the time buffer of +/- three minutes was too large and caught an extra paper read, N=eg. powerlab did not start recording and did not capture the PArt data of the ABG)';
    response = input(prompt,'s');
        if response == 'Y'
            prompt2 = 'Which paper time would you like to ignore? (input a number)';
            response2 = input(prompt2);
            ABGTimes(response2) = [];
            pHvalues(response2) = []; 
            pCO2values(response2) = []; 
            pO2values(response2) = []; 
            sO2values(response2) = [];
            indices(response2) = [];
        elseif response == 'N'
            prompt2 = 'Would you like to artificially add a "PArt" time? [Y/N]';
            response2 = input(prompt2);
            if response2 == 'Y'
                subprompt1 = 'What time should be added? (input an array of [HH mm ss])';
                subresponse1 = input(subprompt1);
                newtime = PArtABGTimes(1);
                newtime.Hour = subresponse1(1);
                newtime.Minute = subresponse1(2);
                newtime.Second = subresponse1(3);
                PArtABGTimes = sort([PArtABGTimes newtime]);
            elseif response2 == 'N'
            end
        end

end

figure('Position', [100 100 500 500])
title('Using ABG Paper Print Times')
yyaxis left
plot(AbsTime, PArtData,'-') 
ylim([0 2.5])
yyaxis right
plot(ABGTimes,pCO2values,'o')
% ylim([30 100])

figure('Position', [700 100 500 500])
title('Using Average Arterial Pressure ABG Times')
yyaxis left
plot(AbsTime, PArtData,'-') 
ylim([0 2.5])
yyaxis right
plot(PArtABGTimes,pCO2values,'o')
% ylim([30 100])

prompt3 = 'Is it appropriate to use these new ABG times derived from the arterial pressure data? [Y/N]';
response3 = input(prompt3,'s');
    if response3 == 'Y'
        if exist('modifiedABGDataCompiled','var') == 0
            modifiedABGDataCompiled = removevars(ABGDataCompiled,{'TestTime'});
            modifiedABGDataCompiled.Properties.VariableNames(1)={'Time'};
        end
        modifiedABGDataCompiled{indices,1} = PArtABGTimes';
        save(ABGfilename,'ABGDataCompiled', 'modifiedABGDataCompiled')
        ABGTimes = PArtABGTimes;
    else
        disp('Keep paper ABG times.')

    end

%%Plot the ABG Data for this section
figure 
title(strcat('ABG Data for',datestr(timeStart)))

subplot(4,1,1)
plot(ABGTimes,pHvalues,'o')
ylabel('pH')

subplot(4,1,2)
plot(ABGTimes,pCO2values,'o')
ylabel('pCO2')

subplot(4,1,3)
plot(ABGTimes,pO2values,'o')
ylabel('pO2')

subplot(4,1,4)
plot(ABGTimes,sO2values,'o')
ylabel('SO2')


end

