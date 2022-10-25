function plotSelectChannels(Title,EKGON,SpO2ON,PArtON,CapnoON, FlowON,PActON,PPlON, PAbON,VolAutoON,VolAbsON,minVentON,fftPActON,fftFlowON, pHON, pCO2ON,pO2ON,sO2ON, EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,VolAbsData,TminVent,PminVent,xFreqFlow,yMagFlow,xFreqPAct,yMagPAct, ABGTimes,pHvalues,pCO2values,pO2values,sO2values,timeStart,timeDur)
%plotSelectChannels Will plot the selected channels for a given selected
%timerange
%   Title = title for that figure
%   EKGON...VolAbsON are the different channels of data. Their value is 1
%   if the channel is desired to be plotted and 0 if the channel is not
%   to be output in the plot.
%   
%   EKGData...VolAbsData is the full data set from a given block of data
%   TminVent = time indices for the minute ventilation calculated via
%   minVentOverTime
%   PminVent = inspiratory minute ventilation values calculated via
%   minVentOverTime
%   timeStart = the starting point of the selection in the data (in seconds)
%   timeDur = the slength of data selected (in seconds)
%   Title = string of the title for the figure


index = [PActON EKGON SpO2ON PArtON CapnoON FlowON  PPlON PAbON VolAutoON VolAbsON minVentON  pHON pCO2ON pO2ON sO2ON fftPActON fftFlowON]; 
refData = [PActData;EKGData;SpO2Data;PArtData;CapnoData;FlowData;PPlData;PAbData;VolAutoData;VolAbsData];
Labels = {'Actuator Pressure [psi]' 'EKG' 'SpO2' 'Arterial Pressure' 'Capnography' 'Flow [L/s]' 'Pleural Pressure' 'Abdominal Pressure' 'Volume (AutoCorrection) [L]' 'Volume (Absolute) [L]' 'Minute Ventilation [L/min]' 'pH' 'pCO2 [mmHg]' 'pO2 [mmHg]' 'sO2 [%]' 'Actuator Pressure FFT Magnitude' 'Flow FFT Magnitude' };
S = sum(index); %total number of channels turned on
ON = find(index); %indices of which channels are turned on
timerange = [(timeStart*1000):((timeStart+timeDur)*1000)]; %time index range specified by timeStart and timeDur
tdata=[timeStart:0.001:timeStart+timeDur]; %time x data for plotting purposes

figure
hold on
title(Title)
    for n = 1:S 
        subplot(S,1,n)
        switch Labels{ON(n)} 
            case'Minute Ventilation [L/min]'
                xdata = TminVent;
                ydata = PminVent;
            case 'Actuator Pressure FFT Magnitude'
                xdata = xFreqPAct;
                ydata = abs(yMagPAct);
            case 'Flow FFT Magnitude'
                xdata = xFreqFlow;
                ydata = abs(yMagFlow);
            case 'pH'
                xdata = ABGTimes;
                ydata = pHvalues;
            case 'pCO2 [mmHg]'
                xdata = ABGTimes;
                ydata = pCO2values;
            case 'pO2 [mmHg]' 
                xdata = ABGTimes;
                ydata = pO2values;
            case 'sO2 [%]'
                xdata = ABGTimes;
                ydata = sO2values;
            otherwise
                xdata = tdata;
                ydata = refData(ON(n),timerange);
        end
        plot(xdata,ydata);
        ylabel(Labels{ON(n)});
        switch Labels{ON(n)} 
            case 'Actuator Pressure [psi]'
                ylim([-0.5 22])                
        end
        

          
        if fftPActON+fftFlowON>0
            if n < S-(fftPActON+fftFlowON)
                xlim([timeStart timeStart+timeDur])
            elseif n == S-(fftPActON+fftFlowON)
                xlabel('Time [s]')
                xlim([timeStart timeStart+timeDur])
            elseif n > S-(fftPActON+fftFlowON) && n<S
                xlim([0 5])
            elseif n == S
                xlabel('Frequency [Hz]')
                xlim([0 5])
            end
        else%none of the ffts are on
            xlim([timeStart timeStart+timeDur])
            if n==S 
                xlabel('Time [s]')
            end
        end
    hold off

    end

