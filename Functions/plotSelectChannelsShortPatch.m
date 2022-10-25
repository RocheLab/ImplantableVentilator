function plotSelectChannelsShortPatch(EKGON,SpO2ON,PArtON,CapnoON, FlowON,PActON,PPlON, PAbON,VolAutoON,VolAbsON,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,VolAbsData, time, Title)
%plotSelectChannels Will plot the selected channels for a given selected
%timerange
%   EKGON...VolAbsON are the different channels of data. Their value is 1
%   if the channel is desired to be plotted and 0 if the channel is not
%   to be output in the plot.
%   
%   EKGData...VolAbsData is the full data set from a given block of data
%   timeStart = the starting point of the selection in the data (in seconds)
%   timeDur = the slength of data selected (in seconds)
%   Title = string of the title for the figure


index = [PActON EKGON SpO2ON PArtON CapnoON FlowON  PPlON PAbON VolAutoON VolAbsON]; 
refData = [PActData;EKGData;SpO2Data;PArtData;CapnoData;FlowData;PPlData;PAbData;VolAutoData.*1000;VolAbsData.*1000];
Labels = {'Actuator Pressure [psi]' 'EKG' 'SpO2' 'Arterial Pressure' 'Capnography' 'Flow [L/s]' 'Pleural Pressure' 'Abdominal Pressure' 'Volume (AutoCorrection) [L]' 'Volume (Absolute) [L]'};
S = sum(index); %total number of channels turned on
ON = find(index); %indices of which channels are turned on
timerange = [(time(1)):(time(end))]; %time index range specified by timeStart and timeDur
xdata=[time]; %time x data for plotting purposes

timepatch = [94.07 714.8];%in s
% v = [0 ymin; timepatch(1) ymin; timepatch(1) ypk; 0 ypk; timepatch(2) ymin; time(end) ymin; time(end) ypk; timepatch(2) ypk];
% f = [1 2 3 4; 5 6 7 8];
% patch('Faces',f,'Vertices',v,'FaceColor',[0.9 0.9 0.9])


% figure
hold on
title(Title)
    for n = 1:S 
        ydata = refData(ON(n),:);
        
        subplot(S,1,n)
%         plot(xdata,ydata);
        
        ylabel(Labels{ON(n)});
        xlim([time(1) time(end)])
        switch Labels{ON(n)} 
            case 'Actuator Pressure [psi]'
                ymin= min(ydata);
                ypk= max(ydata);
                v = [0 ymin; timepatch(1) ymin; timepatch(1) ypk; 0 ypk; timepatch(2) ymin; time(end) ymin; time(end) ypk; timepatch(2) ypk];
                f = [1 2 3 4; 5 6 7 8];
                patch('Faces',f,'Vertices',v,'FaceColor',[0.9 0.9 0.9])
                plot(xdata,ydata);
            case 'Flow [L/s]'
                ymin= min(ydata);
                ypk= max(ydata);
%                 ylim([-0.5 0.4])
                v = [0 ymin; timepatch(1) ymin; timepatch(1) ypk; 0 ypk; timepatch(2) ymin; time(end) ymin; time(end) ypk; timepatch(2) ypk];
                f = [1 2 3 4; 5 6 7 8];
                patch('Faces',f,'Vertices',v,'FaceColor',[0.9 0.9 0.9])
                plot(xdata,ydata);
            case 'Volume (AutoCorrection) [mL]'
                ymin= min(ydata);
                ypk= max(ydata);
                v = [0 ymin; timepatch(1) ymin; timepatch(1) ypk; 0 ypk; timepatch(2) ymin; time(end) ymin; time(end) ypk; timepatch(2) ypk];
                f = [1 2 3 4; 5 6 7 8];
                patch('Faces',f,'Vertices',v,'FaceColor',[0.9 0.9 0.9])
                plot(xdata,ydata);
        end

        
        if n == S
            xlabel('Time [s]')
        end
    end
hold off

end