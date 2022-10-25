%generateFig_Cambpbell_PeakPressures_Revision

clear all
close all

% load('Block_16_USL_Shapes_Pressure_Trial.mat')
% load('Block_19_USR_Pressure_Trial.mat')
% 
load('listRevisionCharacterization.mat')

f1=figure('Position',[100 100 250 250]);

l=0;

% experiment = 12;
% load('Block_19_USR_Pressure_Trial.mat')
% Block = 19;
% for FileNum = [35 37 39 41 33]; %10 s segments
% for FileNum = [34 36 38 40 32]; %20 s segments
%%%%%%%%%%%%OR%%%%%%%%%%%%%
experiment = 12;
load('Block_16_USL_Shapes_Pressure_Trial.mat')
Block = 16;
for FileNum = [31 39 37 35 33]; %10 s segments %for Block 16 % spont 5 10 15 20
% for FileNum = [21 23 25 27 19 ]; %10 s segments % for Block 16

%% Load data
filepath = strcat('Revision Characterization Files/',list.name{FileNum});
load(filepath);


[time,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,VolAbsData,DateTime] = blockAnalysis(data,1,datastart,dataend,blocktimes,titles);
[time,VolAutoData,VolAbsData,newVolAbsData,Correction] = calibrateVolumeRemoveIntegralNoise(time,VolAutoData,VolAbsData);
[newVolAutoData] = spirometryNormalization(time,newVolAbsData);
%Trim data into the selected portion for the file we are loading
[time,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,VolAbsData,newVolAbsData,newVolAutoData] = trimData_exp(timeStart,timeDur,time,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,VolAbsData,newVolAbsData,newVolAutoData);

%%%%%%%%%%%%%%%%%%Filter

%% find HR
%use spo2 OR art pressure to calc heart rate
if isempty(find(SpO2Data>1.75, 1)) %if there are no points where spo2 saturates out
    %use spo2 to determine heartrate
    [pks,locs]=findpeaks(SpO2Data,time,'MinPeakDistance', 0.4);
    HR = length(locs)/(locs(end)-locs(1)); %beats per second
    bpm = HR*60;
elseif isempty(find(PArtData>1.75, 1)) %if there are no points where art pressure saturates out
    [pks,locs]=findpeaks(PArtData,time,'MinPeakDistance', 0.4);
    HR = length(locs)/(locs(end)-locs(1)); %beats per second
    bpm = HR*60;
else
    subrange = 10*1000:20*1000;
    if isempty(find(SpO2Data(subrange)>1.75, 1)) %if there are no points where spo2 saturates out
    %use spo2 to determine heartrate
    [pks,locs]=findpeaks(SpO2Data(subrange),time(subrange),'MinPeakDistance', 0.4);
    HR = length(locs)/(locs(end)-locs(1)); %beats per second
    bpm = HR*60;
    elseif isempty(find(PArtData(subrange)>1.75, 1)) %if there are no points where art pressure saturates out
    [pks,locs]=findpeaks(PArtData(subrange),time(subrange),'MinPeakDistance', 0.4);
    HR = length(locs)/(locs(end)-locs(1)); %beats per second
    bpm = HR*60;
    end
end



range = 1:length(time);

% Filter Signals
% range = 100*1000:800*1000;

fpass = 0.01; %Hz



% filtPPl = lowpass(PPlData,fpass,1000,'Steepness',1); 
fc = HR-0.1;
fs = 1000;

[b,a] = butter(6,fc/(fs/2));
filtPPl = filtfilt(b,a,PPlData);
filtPAb = filtfilt(b,a,PAbData);
altPDi=filtPAb-filtPPl;

% Find breath bounds

[TF_Vol] = islocalmin(newVolAutoData(range),'MinSeparation',1400);
Vbnds = time(range(TF_Vol));


%% Normalize pressures

norm_filtPPl = filtPPl-mean(filtPPl(range(TF_Vol)));
norm_filtPAb = filtPAb-mean(filtPAb(range(TF_Vol)));
norm_filtPDi = norm_filtPAb-norm_filtPPl;










time=(time-time(1));%in s

% plotSelectChannelsShort(EKGON,SpO2ON,PArtON,CapnoON, FlowON,PActON,PPlON, PAbON,VolAutoON,VolAbsON,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,newVolAutoData,newVolAbsData, time, Title)
% plotSelectChannelsShort(EKGON,SpO2ON,PArtON,CapnoON, FlowON,PActON,PPlON, PAbON,VolAutoON,VolAbsON,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,newVolAbsData, timeStart,timeDur, Title)

% plotSelectChannels(Title,EKGON,SpO2ON,PArtON,CapnoON, FlowON,PActON,PPlON, PAbON,VolAutoON,VolAbsON,minVentON,fftPActON,fftFlowON, pHON, pCO2ON,pO2ON,sO2ON, EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,newVolAbsData,TminVent,PminVent,xFreqFlow,yMagFlow,xFreqPAct,yMagPAct, ABGTimes,pHvalues,pCO2values,pO2values,sO2values,timeStart,timeDur)
% end

%dev vs spont FN2

% cutoff = 20*1000;

range=[3*1000:9*1000];

% ColorDev = hex2rgb('#7262C1');
% ColorMech = ColorDev.*1.3;
% ColorSpont = ColorDev.*0.6;
l = l+1;

switch FileNum %[31 39 37 35 33]; %10 s segments %for Block 16 % spont 5 10 15 20
    case {31}%spont
        ColorCode =  [0.3 0.3 0.3];
    case {39}%5psi
        ColorCode = hex2rgb('#C9DEDB');
    case {37}%10
        ColorCode = hex2rgb('#83ABA4');
    case {35}%15
        ColorCode = hex2rgb('#428C80');
    case {33}%20
        ColorCode = hex2rgb('#146759');

end


hold on
plot(norm_filtPPl(range),newVolAutoData(range).*1000,'Color',ColorCode,'LineWidth',1)
xlabel('\Delta Pleural Pressure [cmH_2O]')
ylabel('Tidal Volume [mL]')
P{l} = norm_filtPPl(range);
V{l} = newVolAutoData(range).*1000;


end
hold off

% 
% %dev vs mech FN5
% 
% cutoff = 15.5*1000;
% 
% range1=[1:cutoff];
% range2=[cutoff:length(time)];
% 
% ColorDev = hex2rgb('#7262C1');
% ColorMech = ColorDev.*1.3;
% ColorSpont = ColorDev.*0.6;
% 
% f1=figure('Position',[100 100 300 300]);
% hold on
% plot(norm_filtPPl(range1),newVolAutoData(range1).*1000,'Color',ColorDev)
% plot(norm_filtPPl(range2),newVolAutoData(range2).*1000,'Color',ColorMech)
% xlabel('\Delta Pleural Pressure [cmH_2O]')
% ylabel('Tidal Volume [mL]')
% hold off


%%
%% Save screencapture


pathWithFolderName =  strcat(pwd,'\Figures For Paper\');
Condition = 'Campbell_Pressures_Revision';
figCondition = strcat('',Condition);
figPrefix = strcat('Block',num2str(Block));

figName = strcat(figPrefix,figCondition);
figFileName = strcat(pathWithFolderName,figName);

exportgraphics(f1,strcat(figFileName,'.png')) %will save figure f as a .png
exportgraphics(f1,strcat(figFileName,'.eps'),'ContentType','vector') %will save figure f as a .png

