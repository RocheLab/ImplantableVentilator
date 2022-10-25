function apneaVSmaintenence(date,time,PActData,timeStart,ABGTimes)
%UNTITLED9 Summary of this function goes here
%   time
%   PActData = trimmed actuator pressure data
%   timeStart
%   ABGTimes = from modifiedABGDataCompiled

ABGfilename = [num2str(date) 'ABGData.mat'];
load(ABGfilename)


ind = find(PActData>5); %may be off by 1-2 s

timefromsegstart = time(ind(1))-timeStart

for i = 1:length(ABGTimes)
    
    abgIndex = find(modifiedABGDataCompiled{:, 1} == ABGTimes(i));     
    if ~any(strcmp(modifiedABGDataCompiled.Properties.VariableNames,'ExperimentType')) 
        ExperimentType = string(zeros(length(ABGDataCompiled{:, 1}),1));
        modifiedABGDataCompiled = addvars(modifiedABGDataCompiled,ExperimentType);
    end
    
    if timefromsegstart > 60 %if the delay is greater than one minute, we will consider that an apnea period, and will use the start of the trimmed block as the starttime
        modifiedABGDataCompiled{abgIndex,7}="Apnea";
    else %if the delay is less than 60 seconds, assume no apnea and assert the start time as the start of the Actuator start block    
        modifiedABGDataCompiled{abgIndex,7}="Maintenence";
    end
     
end

        save(ABGfilename,'ABGDataCompiled', 'modifiedABGDataCompiled')


end

