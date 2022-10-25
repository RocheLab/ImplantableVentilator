function [fig, dataSAVE]=plotSelectChannelsABG(Title,EKGON,SpO2ON,PArtON,CapnoON, FlowON,PActON,PPlON, PAbON,VolAutoON,VolAbsON,minVentON,fftPActON,fftFlowON, pHON, pCO2ON,pO2ON,sO2ON, EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,VolAbsData,TminVent,PminVent,xFreqFlow,yMagFlow,xFreqPAct,yMagPAct, ABGTimes,pHvalues,pCO2values,pO2values,sO2values,priorPoint,timeStart,timeDur,blockStart)
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
Labels = {'Actuator Pressure [psi]' 'EKG' 'SpO2' 'Arterial Pressure' 'Capnography' 'Flow [L/s]' 'Pleural Pressure' 'Abdominal Pressure' 'Volume [mL]' 'Volume (Absolute) [L]' 'Minute Ventilation [L/min]' 'pH' 'pCO2 [mmHg]' 'pO2 [mmHg]' 'sO2 [%]' 'Actuator Pressure FFT Magnitude' 'Flow FFT Magnitude' };
S = sum(index); %total number of channels turned on
ON = find(index); %indices of which channels are turned on

if strcmp(datestr(blockStart), '06-Oct-2020 11:18:24') %trims the manual high breaths included at the end of this trial
    timeDur = timeDur-40;
end

timerange = [(timeStart*1000):((timeStart+timeDur)*1000)]; %time index range specified by timeStart and timeDur
tdata=[timeStart:0.001:timeStart+timeDur]; %time x data for plotting purposes
tdata = tdata - timeStart; %normalize time to zero
% Tdata = blockStart + seconds(tdata);%convert to datetimes
ABGTimes = seconds(ABGTimes-blockStart)-timeStart;

fig=figure('Position',[00 50 280 600]);
hold on
title(Title)
    for n = 1:S 
        subplot(S,1,n)
        switch Labels{ON(n)} 
            case'Minute Ventilation [L/min]'
                xdata = TminVent;
                ydata = PminVent;
                dataSAVE{n,1} = xdata;
                dataSAVE{n,2} = ydata;

                plot(xdata,ydata,'o');
            case 'Actuator Pressure FFT Magnitude'
                xdata = xFreqPAct;
                ydata = abs(yMagPAct);
                dataSAVE{n,1} = xdata;
                dataSAVE{n,2} = ydata;
                plot(xdata,ydata,'o');
            case 'Flow FFT Magnitude'
                xdata = xFreqFlow;
                ydata = abs(yMagFlow);
                dataSAVE{n,1} = xdata;
                dataSAVE{n,2} = ydata;
                plot(xdata,ydata,'o');
            case 'pH'
%                 xdata = [Tdata(1); ABGTimes]; %in datetimes
                xdata = [tdata(1); ABGTimes]; %in seconds
                ydata = [priorPoint(1); pHvalues];
                dataSAVE{n,1} = xdata;
                dataSAVE{n,2} = ydata;
                plot(xdata,ydata,'o');
                v = [0 7.35; xdata(end)+100 7.35; xdata(end)+100 7.45 ; 0 7.45];%normal range
                f = [1 2 3 4];
                patch('Faces',f,'Vertices',v,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
%                 text(xdata+0.0001,ydata,cellstr(num2str(ydata)))
                text(xdata+15,ydata,cellstr(num2str(ydata)))
                ylim([7.2 7.5])
       
            case 'pCO2 [mmHg]'
%                 xdata = [Tdata(1); ABGTimes]; %in datetimes
                xdata = [tdata(1); ABGTimes]; %in seconds
                ydata = [priorPoint(2); pCO2values];
                dataSAVE{n,1} = xdata;
                dataSAVE{n,2} = ydata;
                plot(xdata,ydata,'o');
                v = [0 35; xdata(end)+100 35; xdata(end)+100 45 ; 0 45];%normal range
                f = [1 2 3 4];
                patch('Faces',f,'Vertices',v,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
%                 text(xdata+0.0001,ydata,cellstr(num2str(ydata)))
                text(xdata+15,ydata,cellstr(num2str(ydata)))
                ylim([30 90])
            case 'pO2 [mmHg]' 
%                 xdata = [Tdata(1); ABGTimes]; %in datetimes
                xdata = [tdata(1); ABGTimes]; %in seconds
                ydata = [priorPoint(3); pO2values];
                dataSAVE{n,1} = xdata;
                dataSAVE{n,2} = ydata;
                plot(xdata,ydata,'o');

            case 'sO2 [%]'
%                 xdata = [Tdata(1); ABGTimes]; %in datetimes
                xdata = [tdata(1); ABGTimes]; %in seconds
                ydata = [priorPoint(4); sO2values];
                dataSAVE{n,1} = xdata;
                dataSAVE{n,2} = ydata;
                plot(xdata,ydata,'o');
                v = [0 98; xdata(end)+100 98; xdata(end)+100 100 ; 0 100];%normal range
                f = [1 2 3 4];
                patch('Faces',f,'Vertices',v,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
%                 text(xdata+0.0001,ydata,cellstr(num2str(ydata)))
                text(xdata+15,ydata,cellstr(num2str(ydata)))

            otherwise 
                Xdata = tdata;
                Ydata = refData(ON(n),timerange);

%                 [pks,locs]=findpeaks(Ydata,Xdata,'MinPeakDistance',0.000019);%in datetimes
                [pks,locs]=findpeaks(Ydata,Xdata,'MinPeakDistance',1.5);%in seconds
                xdata = locs;
                ydata = pks;
                dataSAVE{n,1} = xdata;
                dataSAVE{n,2} = ydata;
                p=plot(xdata,ydata,'.');
        end

        ylabel(Labels{ON(n)});
        switch Labels{ON(n)} 
            case 'Actuator Pressure [psi]'
                ylim([-0.5 22])                
        end
        switch Labels{ON(n)} 
            case 'Actuator Pressure [psi]'
                ylabel({'Actuator';'Pressure [psi]'})
                p.Color = hex2rgb('60B893');
            case 'Flow [L/s]'
                 ylim([-0.5 0.5])
                p.Color = hex2rgb('4B7CB7');
            case 'Volume [mL]'
                p.Color = hex2rgb('F6C433');
                 ylim([0 200])
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
            xlim([tdata(1) tdata(end)])
            if n==S 
                xlabel('Time [s]')
            end
        end
    hold off

    
    end

