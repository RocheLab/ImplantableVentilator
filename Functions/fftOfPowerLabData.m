function [xFreqFlow,yDFTFlow,xFreqPAct,yDFTPAct] = fftOfPowerLabData(FlowData,PActData,timeStart,timeDur)
%fftOfPowerLabData Conduct a fft of a snippet of data bound by timeStart
%and timeStart+timeDur. 
%   FlowData = flow data from powerlab
%   PActData = actuator pressurization data on powerlab

timerange = [(timeStart*1000)+1:((timeStart+timeDur)*1000)];

yDFTFlow = fft(FlowData(timerange));
fs = 1/(0.001);
xFreqFlow = (0:length(yDFTFlow)-1)*fs/length(yDFTFlow);

yDFTPAct = fft(PActData(timerange));
xFreqPAct = (0:length(yDFTPAct)-1)*fs/length(yDFTPAct);





end

