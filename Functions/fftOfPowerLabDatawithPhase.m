function [xFreqFlow,yDFTFlow,xFreqPAct,yDFTPAct] = fftOfPowerLabDataTRIMMED(FlowData,PActData)
%fftOfPowerLabData Conduct a fft of a snippet of data bound by timeStart
%and timeStart+timeDur. 
%   FlowData = flow data from powerlab
%   PActData = actuator pressurization data on powerlab


yDFTFlow = fft(FlowData);
fs = 1/(0.001);
xFreqFlow = (0:length(yDFTFlow)-1)*fs/length(yDFTFlow);

yDFTPAct = fft(PActData);
xFreqPAct = (0:length(yDFTPAct)-1)*fs/length(yDFTPAct);





end

