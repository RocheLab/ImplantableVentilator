function [PVT,NVT] = calcMinuteVentilation(FlowData,timeStart,timeDur,plotON,zeroshift,scale)
%calcMinuteVentilation calculates the minute ventilation for a given
%segment of flow data
%   FlowData = full section of data from flow channel
%   timeStart = the starting point of the selection in the data (in seconds)
%   timeDur = the slength of data selected (in seconds)
%   plotON = binary, specifies whether or not you want the function to
%   output a plot



timerange = [(timeStart*1000):((timeStart+timeDur)*1000)-1];

flow = FlowData(timerange);
smoothed = smooth(flow);
shifted=flow-zeroshift;
pos = [];
neg = [];
for i = [1:length(flow)]
    if shifted(i)>=0
        pos(i) = shifted(i);
        neg(i) = 0;
    else
        pos(i) = 0;
        neg(i) = shifted(i);
    end
end
Pint = cumtrapz(pos);
PVT = Pint(end)*60*scale/(timeDur);
Nint = cumtrapz(neg);
NVT = Nint(end)*60*scale/(timeDur);

t=[timeStart:0.001:timeStart+timeDur-0.001];



if plotON == 1
    figure
subplot(5,1,1)
plot(t,flow)
xlabel('Time [s]')
ylabel('Flow [L/s]')

subplot(5,1,2)
plot(t,smoothed)
xlabel('Time [s]')
ylabel('Smoothed Flow [L/s]')

subplot(5,1,3)
plot(t,shifted)
xlabel('Time [s]')
ylabel('Shifted Flow [L/s]')

subplot(5,1,4)
plot(t,pos)
xlabel('Time [s]')
ylabel('Inspiratory Flow [L/s]')

subplot(5,1,5)
plot(t,neg)
xlabel('Time [s]')
ylabel('Expiratory Flow [L/s]')
    
end

end

