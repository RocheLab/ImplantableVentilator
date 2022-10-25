function [TminVent,PminVent,NminVent] = minVentOverTimeTRIMMED(FlowData,timeStart,timeDur,dt,plotON,smoothON, shiftON,calPlotON, zeroshift,scale)
%minVentOverTime calculates the minute ventilation for a given
%segment of flow data for specific divisions of time (dt)
%INPUTS
%   FlowData = full section of data from flow channel
%   timeStart = the starting point of the selection in the data (in seconds)
%   timeDur = the length of data selected (in seconds)
%   dt = timestep unit the dynamic minute ventilation is over (s)
%   plotON = binary, specifies whether function will plot flow and minute
%   ventilation plot
%   smoothON = binary, default is 0, will determine whether to apply a
%   smoothing function to data
%   shiftON = binary, default is 0, will determine whether to apply a
%   zeroshift to data
%   calPlotON = binary, default is 0, useful plot for 
%   zeroshift = default is 0, if shiftON is 1, then we will adjust the zero
%   point based off of this value
%   scale = default is 0.001 to comvert from mL to L, useful for calibration
%OUTPUTS
%   TminVent = a time vector for the minute vent calculation with each T
%   at end of each dt bucket
%   PminVent = a vector of the minute ventilation calculated off of
%   inspiration over timesteps of dt
%   NminVent = a vector of the minute ventilation calculated off of
%   exspiration over timesteps of dt

Ndt = ceil(timeDur/dt);
overTime = dt*Ndt;
limit = roundn(length(FlowData),4);
data = FlowData;

%% data modification
Nsubplots = 3;
if smoothON == 1%smoothing
    data = smooth(data);
    Nsubplots = Nsubplots + 1;
end

if shiftON == 1%apply a zeroshift
    data = data - zeroshift;
    Nsubplots = Nsubplots + 1;
end

%% Cumulative Trapezoidal integration
%Because in terms of vntilation, the delivered air is what is of concern,
%we focus on the positive minute ventilation, neg ventilation is useful for
%troubleshooting and calibration purposes
pos = [];
neg = [];
for i = [1:length(data)]
    if data(i)>=0
        pos(i) = data(i);
        neg(i) = 0;
    else
        pos(i) = 0;
        neg(i) = data(i);
    end
end

trimfactor = floor(length(pos)/(dt*1000));

dataIntP = reshape(pos(1:trimfactor*dt*1000),dt*1000,[]);%data reshaped into an array to perform a cumulative trapezoidal integration over
% dataIntP = reshape(pos,dt*1000,[]);%data reshaped into an array to perform a cumulative trapezoidal integration over
PInt = cumtrapz(dataIntP);
PminVent = PInt(end,:).*(60*scale/(dt));

% dataIntN = reshape(neg(1:trimfactor*dt*1000),[dt*1000,Ndt]);
dataIntN = reshape(neg(1:trimfactor*dt*1000),dt*1000,[]);%data reshaped into an array to perform a cumulative trapezoidal integration over
NInt = cumtrapz(dataIntN);
NminVent = NInt(end,:).*(60*scale/(dt));

TminVent = [timeStart+(dt):dt:timeStart+overTime]; %time vector for plotting purposes

t=[timeStart:0.001:timeStart+overTime-0.001]; %for plotting against data
xLim = [timeStart,timeStart+timeDur];
%% Plotting


if plotON == 1
figure
%plot to see Min vent per dt plotted with flow
subplot(2,1,1)
plot(t,data)
xlabel('Time [s]')
ylabel('Flow [L/s]')
xlim(xLim)

subplot(2,1,2)
plot(TminVent,PminVent)
xlabel('Time [s]')
ylabel('Minute Ventilation [L/min]')
xlim(xLim)
end

if calPlotON ==1
% Plot for troubleshooting and calibrating    
figure
subplot(Nsubplots,1,1)
plot(t,data)
xlabel('Time [s]')
ylabel('Flow [L/s]')
xlim(xLim)

subplot(Nsubplots,1,2)
plot(t,pos)
xlabel('Time [s]')
ylabel('Inspiratory Flow [L/s]')
xlim(xLim)

subplot(Nsubplots,1,3)
plot(t,neg)
xlabel('Time [s]')
ylabel('Expiratory Flow [L/s]')
xlim(xLim)

if smoothON == 1
subplot(Nsubplots,1,4)
plot(t,smoothed)
xlabel('Time [s]')
ylabel('Smoothed Flow [L/s]')
xlim(xLim)
end

if shiftON == 1
subplot(Nsubplots,1,5)
plot(t,shifted)
xlabel('Time [s]')
ylabel('Shifted Flow [L/s]')
xlim(xLim)
end

    
end



end

