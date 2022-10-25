function [starttime] = findStartTime(time,PActData,DateTime,ABGRefTime)
%findStartTime will output the "starttime" of a block of data by looking
%for the start of the actuation data
%   time = trimmed vector of time data
%   PActData = trimmed vector of actuator data
%   DateTime = datetime of the start of the powerlab block of data that
%       time and PActData comefrom 
%   ABGRefTime = 1 ABG time that this function will reference


ind = find(PActData>3); %may be off by 1-2 s
ActStart = DateTime + seconds(time(ind(1)));

timefromsegstart = ABGRefTime-ActStart;

if timefromsegstart > 60 %if the delayis greater than one minute, we will consider that an apnea period, and will use the start of the trimmed block as the starttime
    starttime = DateTime + seconds(time(1));
else %if the delay is less than 60 seconds, assume no apnea and assert the start time as the start of the Actuator start block    
    starttime = ActStart;
end



end

