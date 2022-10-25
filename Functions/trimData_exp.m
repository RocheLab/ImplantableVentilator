function [time,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,VolAbsData,newVolAbsData,newVolAutoData,DateTime] = trimData_DT(timeStart,timeDur,time,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,VolAbsData,newVolAbsData,newVolAutoData,DateTime)
%trimData Will trim data channels to specified timeStart and timeDur
%parameters
timerange = [(round(timeStart*1000))+1:floor((timeStart+timeDur)*1000)];
time = time(timerange);
EKGData = EKGData(timerange);
if isempty(SpO2Data)
    SpO2Data = [];
else
SpO2Data = SpO2Data(timerange);
end
PArtData = PArtData(timerange);
CapnoData = CapnoData(timerange);
FlowData = FlowData(timerange);
PActData = PActData(timerange);
PPlData = PPlData(timerange);
PAbData = PAbData(timerange);
VolAutoData = VolAutoData(timerange);
VolAbsData = VolAbsData(timerange);
newVolAbsData = newVolAbsData(timerange);
newVolAutoData = newVolAutoData(timerange);


end

