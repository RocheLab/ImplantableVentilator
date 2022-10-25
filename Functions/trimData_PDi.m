function [time,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,PDiData,VolAutoData,VolAbsData,newVolAbsData,newVolAutoData] = trimData_PDi(timeStart,timeDur,time,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,PDiData,VolAutoData,VolAbsData,newVolAbsData,newVolAutoData)
%trimData Will trim data channels to specified timeStart and timeDur
%parameters
timerange = [(round(timeStart*1000))+1:floor((timeStart+timeDur)*1000)];
time = time(timerange);
EKGData = EKGData(timerange);
SpO2Data = SpO2Data(timerange);
PArtData = PArtData(timerange);
CapnoData = CapnoData(timerange);
FlowData = FlowData(timerange);
PActData = PActData(timerange);
PPlData = PPlData(timerange);
PAbData = PAbData(timerange);
PDiData = PDiData(timerange);
VolAutoData = VolAutoData(timerange);
VolAbsData = VolAbsData(timerange);
newVolAbsData = newVolAbsData(timerange);
newVolAutoData = newVolAutoData(timerange);


end

