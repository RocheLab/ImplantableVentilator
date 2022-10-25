function [newVolAutoData] = spirometryNormalization(time,newVolAbsData)
%spirometryNormalization This will break a segment down breath by breath
%and normalize the breath to reset at zero 

Correction = [];
%find minima

Vminind = islocalmin(newVolAbsData,'FlatSelection','last','MinSeparation',500,'MinProminence',0.01);
Vminbounds = find(Vminind);

for i = 2:(length(Vminbounds))
    if i==2
       y = [newVolAbsData(Vminbounds(1)) newVolAbsData(Vminbounds(2))]';
       x = [time(Vminbounds(1)) time(Vminbounds(2))]';
       X = [ones(length(x),1) x]; %X for the linear regression, ones column is for y intercept
       B = X\y;%B is the regression coefficient (slope)inter = X*B;
       T = [time(1:Vminbounds(2))]';
       timeint = [ones(length(T),1) T];
       interp = timeint*B;
       
       Correction = [Correction interp'];
    elseif i>2 && i<(length(Vminbounds))
       y = [newVolAbsData(Vminbounds(i-1)) newVolAbsData(Vminbounds(i))]';
       x = [time(Vminbounds(i-1)) time(Vminbounds(i))]';
       X = [ones(length(x),1) x]; %X for the linear regression, ones column is for y intercept
       B = X\y;%B is the regression coefficient (slope)inter = X*B;
       T = [time((Vminbounds(i-1)+1):Vminbounds(i))]';
       timeint = [ones(length(T),1) T];
       interp = timeint*B;
       
       Correction = [Correction interp'];
    elseif i==(length(Vminbounds))
       y = [newVolAbsData(Vminbounds(i-1)) newVolAbsData(Vminbounds(i))]';
       x = [time(Vminbounds(i-1)) time(Vminbounds(i))]';
       X = [ones(length(x),1) x]; %X for the linear regression, ones column is for y intercept
       B = X\y;%B is the regression coefficient (slope)inter = X*B;
       T = [time((Vminbounds(i-1)+1):end)]';      
       timeint = [ones(length(T),1) T];
       interp = timeint*B;
       Correction = [Correction interp'];
    end    
end


if length(Correction) == length(newVolAbsData)
    newVolAutoData = newVolAbsData-Correction;
end

% plot(time,newVolAutoData,time,newVolAbsData)

end

