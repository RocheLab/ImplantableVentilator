function [time,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,VolAbsData,PDiData,DateTime] = blockAnalysis(data,num,datastart,dataend,blocktimes,titles)
%blockAnalysis This function parses data for a specified block
%   data = data from exported powerlab file
%   num = which block to select from this exported powerlab file
%   datastart = data start indices from exported powerlab file
%   dataend = data end indices from exported powerlab file
%   blocktimes = time indices for blocks from exported powerlab file
%   titles = channel names from exported powerlab file

    Channel(1,:) = data(datastart(1,num):dataend(1,num));
    Channel(2,:) = data(datastart(2,num):dataend(2,num));
    Channel(3,:) = data(datastart(3,num):dataend(3,num));
    Channel(4,:) = data(datastart(4,num):dataend(4,num));
    Channel(5,:) = data(datastart(5,num):dataend(5,num));
    Channel(6,:) = data(datastart(6,num):dataend(6,num));
    Channel(7,:) = data(datastart(7,num):dataend(7,num));
    Channel(8,:) = data(datastart(8,num):dataend(8,num));
    Channel(9,:) = data(datastart(9,num):dataend(9,num));
    Channel(10,:) = data(datastart(10,num):dataend(10,num));
    if size(datastart,1)>10
        if datastart(11,1) > 0
    Channel(11,:) = data(datastart(11,num):dataend(11,num));
        end
        if size(datastart,1)>11
            if datastart(12,1)>0
                Channel(12,:) = data(datastart(12,num):dataend(12,num));
            end
        end
    end
    time = 0:(length(Channel(1,:))-1);
    time = time.*0.001; %time in seconds
% for loop to match all the channels to their appropriate data names
    for i=1:size(titles,1)
        name = strtrim(titles(i,:));
        switch name
            case {'EKG','ECG','Channel5'}
                EKGData = Channel(i,:);
            case {'SpO2','SPO2','spO2','spo2','Spo2'}
                SpO2Data = Channel(i,:);
            case {'Arterial Pressure','ArterialPressure','P_art','P_Art','P_Arterial','ArtPressure','Arterial Line'}
                PArtData = Channel(i,:);
            case {'Capnography','capnography','etco2','EtCO2','etCO2','ETCO2','Capno','CO2Capno'}
                CapnoData = Channel(i,:);
            case {'Flow','flow'}
                FlowData = Channel(i,:);
            case {'Actuator Pressure','P_Actuator','P_actuator','ActuatorPressure','P_act'}
                PActData = Channel(i,:);
            case {'P_pl', 'P_Pl','P_pleural'}
                PPlData = Channel(i,:);
            case {'P_ab', 'P_Ab','P_abdominal'}
                PAbData = Channel(i,:);
            case {'Volume (Corrected)'}
                VolAutoData = Channel(i,:);
            case {'Volume with drift'}
                VolAbsData = Channel(i,:);
            case {'P_di', 'Channel 11'}
                if size(Channel,1) >10
                PDiData = Channel(i,:);
                else
                    PDiData = [];
                end
            case {'Trigger', 'Channel 12'} %this is saved but not output in the function, you can add that if you want
                if size(Channel,1) >11
                TriggerData = Channel(i,:);
                else
                    TriggerData = [];
                end
            otherwise
                warning('Label does not match any expected category')
                disp(name)
        end
    end 
    
    if ~exist('SpO2Data','var')
        SpO2Data = [];
    end
    
    DateTime = datetime(blocktimes(num),'ConvertFrom','datenum');
    %% If you use more than 12 variables, you may need to edit this code
    if size(datastart,1) == 11
        save('blockdata.mat','time','EKGData','SpO2Data','PArtData','CapnoData','FlowData','PActData','PPlData','PAbData','VolAutoData','VolAbsData','PDiData','DateTime')
    elseif size(datastart,1) == 10
        PDiData = [];
        save('blockdata.mat','time','EKGData','SpO2Data','PArtData','CapnoData','FlowData','PActData','PPlData','PAbData','VolAutoData','VolAbsData','DateTime')
    elseif size(datastart,1) == 12
        save('blockdata.mat','time','EKGData','SpO2Data','PArtData','CapnoData','FlowData','PActData','PPlData','PAbData','VolAutoData','VolAbsData','PDiData','TriggerData','DateTime')

    end
    
end

