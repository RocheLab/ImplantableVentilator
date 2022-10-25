function [dataSAVE] = plotSelectChannelsMS(EKGON,SpO2ON,PArtON,CapnoON, FlowON,PActON,PPlON, PAbON,VolAutoON,VolAbsON,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,VolAbsData, time, Title)
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
Labels = {'Actuator Pressure [psi]' 'EKG' 'SpO2' 'Arterial Pressure' 'Capnography' 'Flow [L/s]' 'Pleural Pressure' 'Abdominal Pressure' 'Volume (AutoCorrection) [mL]' 'Volume (Absolute) [L]'};
S = sum(index); %total number of channels turned on
ON = find(index); %indices of which channels are turned on
timerange = [(time(1)):(time(end))]; %time index range specified by timeStart and timeDur
xdata=[time].*1000; %time x data for plotting purposes

% figure
hold on
title(Title)
    for n = 1:S 
        ydata = refData(ON(n),:);
        subplot(S,1,n)
%         plot(xdata,ydata);

        switch Labels{ON(n)} 
            case 'Actuator Pressure [psi]'
                %Pstart
                Pmin_ind = islocalmin(PActData,'MaxNumExtrema',2,'MinSeparation',1500,'FlatSelection','last');
%                 subplot(S,1,n)
                plot(xdata,ydata,'-',time(Pmin_ind).*1000,PActData(Pmin_ind),'o')
%                 ylabel(Labels{ON(n)});
                xlim([time(1)*1000 time(end)*1000])
                dataSAVE{1,1} = xdata;
                dataSAVE{1,2} = ydata;
                dataSAVE{1,3} = time(Pmin_ind).*1000;
                dataSAVE{1,4} = PActData(Pmin_ind);
            case 'Flow [L/s]'
%                 ylim([-0.5 0.4])
                %fpk
                [F_M,F_I] = max(FlowData);
%                 subplot(S,1,n)
                plot(xdata,ydata,'-',time(F_I)*1000,F_M,'o')
%                 ylabel(Labels{ON(n)});
                xlim([time(1)*1000 time(end)*1000])
                dataSAVE{2,1} = xdata;
                dataSAVE{2,2} = ydata;
                dataSAVE{2,3} = time(F_I)*1000;
                dataSAVE{2,4} = F_M;
            case 'Volume (AutoCorrection) [mL]'
%                                 %Vmins
                Vmin_ind = islocalmin(VolAutoData,'MaxNumExtrema',2,'MinSeparation',1500);
%                 Vpk
                [V_M,V_I] = max(VolAutoData);
% %                 subplot(S,1,n)
                plot(xdata,ydata,'-',time(V_I)*1000,V_M*1000,'o',time(Vmin_ind).*1000,VolAutoData(Vmin_ind).*1000,'o')
%                 ylabel(Labels{ON(n)});
                xlim([time(1)*1000 time(end)*1000])
                dataSAVE{3,1} = xdata;
                dataSAVE{3,2} = ydata;
                dataSAVE{3,3} = time(V_I)*1000;
                dataSAVE{3,4} = V_M*1000;
                dataSAVE{3,5} = time(Vmin_ind).*1000;
                dataSAVE{3,6} = VolAutoData(Vmin_ind).*1000;
        end
        if n == S
            xlabel('Time [ms]')
        end
    end
hold off

end