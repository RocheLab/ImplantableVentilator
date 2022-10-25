function plotSelectChannelsTRIMMED(Title,EKGON,SpO2ON,PArtON,CapnoON, FlowON,PActON,PPlON, PAbON,VolAutoON,VolAbsON,minVentON,fftPActMagON,fftPActPhaseON,fftFlowMagON,fftFlowPhaseON, coherenceON, pHON, pCO2ON,pO2ON,sO2ON, EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,VolAbsData,TminVent,PminVent,xFreqFlow,yDFTFlow,xFreqPAct,yDFTPAct, ABGTimes,pHvalues,pCO2values,pO2values,sO2values,DateTime,time)
%plotSelectChannelsTRIMMED Will plot the selected channels for a given selected
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
%   DateTime = datetime of the start of thepowerlab block this data was
%   derived from
%   time = a trimmed time vector of what portion of the powerlab data to use

%   Title = string of the title for the figure


index = [PActON EKGON SpO2ON PArtON CapnoON FlowON  PPlON PAbON VolAutoON VolAbsON minVentON  pHON pCO2ON pO2ON sO2ON fftPActMagON fftPActPhaseON fftFlowMagON fftFlowPhaseON coherenceON]; 
refData = [PActData;EKGData;SpO2Data;PArtData;CapnoData;FlowData;PPlData;PAbData;VolAutoData;VolAbsData];
Labels = {'Actuator Pressure [psi]' 'EKG' 'SpO2' 'Arterial Pressure' 'Capnography' 'Flow [L/s]' 'Pleural Pressure' 'Abdominal Pressure' 'Volume (AutoCorrection) [L]' 'Volume (Absolute) [L]' 'Minute Ventilation [L/min]' 'pH' 'pCO2 [mmHg]' 'pO2 [mmHg]' 'sO2 [%]' 'Actuator Pressure FFT Power [dB]' 'Actuator Pressure FFT Phase' 'Flow FFT Power [dB]' 'Flow FFT Phase' 'Actuation & Flow Coherence'};
S = sum(index); %total number of channels turned on
ON = find(index); %indices of which channels are turned on
% timerange = [(timeStart*1000):((timeStart+timeDur)*1000)]; %time index range specified by timeStart and timeDur
tdata=DateTime + seconds(time); %time x data for plotting purposes
T = length(time)/1000; %period of sampled data
dt = 0.001;

figure('Position',[100 5 500 980])
hold on

    for n = 1:S 
        subplot(S,1,n)
        switch Labels{ON(n)} 
            case'Minute Ventilation [L/min]'
                xdata = TminVent;
                ydata = PminVent;
            case 'Actuator Pressure FFT Power [dB]'
                Sxx = 2*dt^2/T * yDFTPAct.*conj(yDFTPAct);			    % Compute the power spectrum. 
                Sxx=10*log10(Sxx);                                  % ... convert to a decibel scale.
                Sxx = Sxx(1:length(time)/2+1);                         % ... ignore negative frequencies.	
                ydata = 2*0.001^2/T*yDFTPAct.*conj(yDFTPAct);  
%                 df = 1/T;                                           % Determine the frequency resolution. 
%                 fNQ=1/dt/2;                                         % Determine the Nyquist frequency. 
%                 xdata = (0:df:fNQ);  
                xdata = xFreqPAct;
            case 'Flow FFT Power [dB]'
                Sxx = 2*dt^2/T * yDFTPAct.*conj(yDFTPAct);			    % Compute the power spectrum. 
                Sxx=10*log10(Sxx);                                  % ... convert to a decibel scale.
                Sxx = Sxx(1:length(time)/2+1);                         % ... ignore negative frequencies.	
                ydata = 2*0.001^2/T*yDFTFlow.*conj(yDFTFlow);  
                xdata = xFreqFlow;
            case 'Actuator Pressure FFT Phase'
                xdata = xFreqPAct;
                ydata = atan2(imag(yDFTPAct),real(yDFTPAct)).*180/pi;
            case 'Flow FFT Phase'
                xdata = xFreqFlow;
                ydata = atan2(imag(yDFTFlow),real(yDFTFlow)).*180/pi;
            case 'Actuation & Flow Coherence'

                Sxx = 2*dt^2/T * yDFTFlow .* conj(yDFTFlow);          % Power spectrum of x.
                Syy = 2*dt^2/T * yDFTPAct .* conj(yDFTPAct);          % Power spectrum of y.
                Sxy = 2*dt^2/T * yDFTFlow .* conj(yDFTPAct);          % Cross spectrum between x and y.

                Sxx = Sxx(1:length(time)/2+1);                                 % Ignore the negative frequencies.
                Syy = Syy(1:length(time)/2+1);
                Sxy = Sxy(1:length(time)/2+1);

                ydata = abs(Sxy) ./ (sqrt(Sxx) .* sqrt(Syy));        % Compute the coherence.

                df = 1/T;                                           % Determine the frequency resolution.
                fNQ = 1/dt/2;                                       % Determine the Nyquist frequency.
                xdata = (0:df:fNQ);    
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
                ydata = refData(ON(n),:);
        end
        tlim = [tdata(1) tdata(end)];
        switch Labels{ON(n)} 
        case 'pH'
            plot(tlim,(7.35)*ones(1,2),'c-',tlim,(7.45)*ones(1,2),'c-',xdata,ydata,'ok');
        case 'pCO2 [mmHg]'
           plot(tlim,(35)*ones(1,2),'c-',tlim,(45)*ones(1,2),'c-',xdata,ydata,'ok');
        case 'pO2 [mmHg]'
           plot(xdata,ydata,'ok')
        case 'sO2 [%]'
            plot(tlim,(95)*ones(1,2),'c-',tlim,(100)*ones(1,2),'c-',xdata,ydata,'ok');
        otherwise
            plot(xdata,ydata);
        end
        
        sgtitle(Title)
        
        ylabel(Labels{ON(n)});
        
        switch Labels{ON(n)} 
            case 'Actuator Pressure [psi]'
                ylim([-0.5 22])                
        end
        

          
        if fftPActMagON+fftFlowMagON+fftPActPhaseON+fftFlowPhaseON+coherenceON>0
            if n < S-(fftPActMagON+fftFlowMagON+fftPActPhaseON+fftFlowPhaseON+coherenceON)
                xlim([tdata(1) tdata(end)])
            elseif n == S-(fftPActMagON+fftFlowMagON+fftPActPhaseON+fftFlowPhaseON+coherenceON)
                xlabel('Time [s]')
                xlim([tdata(1) tdata(end)])
            elseif n > S-(fftPActMagON+fftFlowMagON+fftPActPhaseON+fftFlowPhaseON+coherenceON) && n<S
                xlim([0 4])                 
            elseif n == S
                xlabel('Frequency [Hz]')
                xlim([0 4])
            end
            
        else%none of the ffts are on
            xlim([tdata(1) tdata(end)])
            if n==S 
                xlabel('Time and Date')
            end
        end
        
       
        
        hold off

    end
    
  saveas(gcf,strcat('FFTTest',Title,'.png'))   
    
end

