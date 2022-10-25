%generateFig_SampleTransDiaphragmaticPressure

%test_transdiaphragm


close all 
clear all

for FileNum = [31 32 35];
% FileNum = 31;
%31: Fig7e
%32: Fig7f
%35: Fig7d

load('listPressureSelection.mat') %read this file to look at names of files and the associated notes
load('fileBlockInfo_C1.mat')



%% Load data
filepath = strcat('Pressure Selection Files/',list.name{FileNum});

load(filepath);
experiment = ExpN;
save('targetBlock.mat','experiment','Block')
[date,x,y,num] = selectBlock(Block,experiment,fileBlockInfo);
[blocktimes,com,comtext,data,dataend,datastart,firstsampleoffset,rangemax,rangemin,samplerate,tickrate,titles,unittext,unittextmap,filename] = loadPowerLabFile(date,x,y);
[time,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,VolAbsData,PDiData,DateTime] = blockAnalysis(data,num,datastart,dataend,blocktimes,titles);
[time,VolAutoData,VolAbsData,newVolAbsData,Correction] = calibrateVolumeRemoveIntegralNoise(time,VolAutoData,VolAbsData);
[newVolAutoData] = spirometryNormalization(time,newVolAbsData);
%Trim data into the selected portion for the file we are loading
[time,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,PDiData,VolAutoData,VolAbsData,newVolAbsData,newVolAutoData] = trimData_PDi(timeStart,timeDur,time,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,PDiData,VolAutoData,VolAbsData,newVolAbsData,newVolAutoData);


%% find HR
%use spo2 OR art pressure to calc heart rate
if isempty(find(SpO2Data>1.75)) %if there are no points where spo2 saturates out
    %use spo2 to determine heartrate
    [pks,locs]=findpeaks(SpO2Data,time,'MinPeakDistance', 0.4);
    HR = length(locs)/(locs(end)-locs(1)); %beats per second
    bpm = HR*60;
elseif isempty(find(PArtData>1.75)) %if there are no points where art pressure saturates out
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





%% Filter Signals
% range = 100*1000:800*1000;
% range = 600*1000:610*1000;
range = 11*1000:23*1000;
fpass = 0.01; %Hz



% filtPPl = lowpass(PPlData,fpass,1000,'Steepness',1); 
fc = HR-0.1;
fs = 1000;

[b,a] = butter(6,fc/(fs/2));
filtPPl = filtfilt(b,a,PPlData);
filtPAb = filtfilt(b,a,PAbData);
altPDi=filtPAb-filtPPl;

% 
% figure('Position',[100 100 800 800]);
% subplot(5,1,1)
% plot(time(range),PPlData(range),time(range),filtPPl(range));
% 
% subplot(5,1,2)
% plot(time(range),PAbData(range),time(range),filtPAb(range));
% 
% subplot(5,1,3)
% plot(time(range),PDiData(range),time(range),filtPDi(range));
% 
% subplot(5,1,4)
% plot(time(range),filtPDi(range),time(range),filtPAb(range)-filtPPl(range))
% 
% subplot(5,1,5)
% plot(time(range),FlowData(range));



%% Normalize pressures

[pks]=findpeaks(filtPPl(range),time(range),'MinPeakDistance', 1);
norm_filtPPl = filtPPl-mean(pks);

[TF_Ab] = islocalmin(filtPAb(range),'MinSeparation',1400);
norm_filtPAb = filtPAb-mean(filtPAb(range(TF_Ab)));
figure

plot(time(range),filtPAb(range), time(range(TF_Ab)),filtPAb(range(TF_Ab)),'*')



%% Find breath bounds for patches

[TF_Vol] = islocalmin(newVolAutoData(range),'MinSeparation',1400);
Vbnds = time(range(TF_Vol));
plot(time(range),newVolAutoData(range),'-',time(range(TF_Vol)),newVolAutoData(range(TF_Vol)),'*')

%% Plot PPl, PAb, and Pdi with breath bounds 

ColorPl = hex2rgb('#7262C1');
ColorAb = hex2rgb('#F15FA6');
ColorDi = hex2rgb('#E87C3B');
ColorVol = hex2rgb('#F6C433');
ColorFlow = hex2rgb('4B7CB7');


switch FileNum
    case 31 %dev
        ColMod = 1;
    case 32 %spont
        ColMod = 1;
    case 35 %mech
        ColMod = 1;
end

f1 = figure('Position',[00 335 220 280]);
subplot(5,1,1)
hold on

Vbnds0=Vbnds-time(range(1));%normalize time
time0 = time - time(range(1));%normalize time 
% 
% y1=min(norm_filtPPl(range)-0.1);
% y2=max(norm_filtPPl(range)+0.1);
y1 = -2;
y2 = 1.5;
[v,f] = createBreathPatchGuidelines(Vbnds0,y1,y2,time0(range(end)));
patch('Faces',f,'Vertices',v,'FaceColor',[0.95 0.95 0.95],'LineStyle','none')

plot(time0(range),norm_filtPPl(range),'Color',ColorPl.*ColMod)
xlim([time0(range(1)) time0(range(end))])
xticks([0 5 10])
ylabel('P_P_l [cmH_2O]')
ylim([y1 y2])
set(gca, 'Layer', 'top')
hold off

subplot(4,1,2)
hold on
% y1=min(norm_filtPAb(range)-0.1);
% y2=max(norm_filtPAb(range)+0.1);
y1 = -0.2;
y2 = 1.5;
[v,f] = createBreathPatchGuidelines(Vbnds0,y1,y2,time0(range(end)));
patch('Faces',f,'Vertices',v,'FaceColor',[0.95 0.95 0.95],'LineStyle','none')
plot(time0(range),norm_filtPAb(range),'Color',ColorAb.*ColMod)
ylabel('P_A_b [cmH_2O]')
xlim([time0(range(1)) time0(range(end))])
xticks([0 5 10])
ylim([y1 y2])
set(gca, 'Layer', 'top')
hold off

subplot(4,1,3)
hold on
norm_filtPDi = norm_filtPAb(range)-norm_filtPPl(range);
% y1=min(norm_filtPDi-0.1);
% y2=max(norm_filtPDi+0.1);
y1 = -1.5;
y2 = 3;
[v,f] = createBreathPatchGuidelines(Vbnds0,y1,y2,time0(range(end)));
patch('Faces',f,'Vertices',v,'FaceColor',[0.95 0.95 0.95],'LineStyle','none')
plot(time0(range),norm_filtPDi,'Color',ColorDi.*ColMod)
ylabel('P_D_i [cmH_2O]')
xlim([time0(range(1)) time0(range(end))])
xticks([0 5 10])
ylim([y1 y2])
set(gca, 'Layer', 'top')
hold off

subplot(4,1,4)
hold on
% y1=min(FlowData(range)-0.05);
% y2=max(FlowData(range)+0.05);
y1 = -1;
y2 = 0.6;
[v,f] = createBreathPatchGuidelines(Vbnds0,y1,y2,time0(range(end)));
patch('Faces',f,'Vertices',v,'FaceColor',[0.95 0.95 0.95],'LineStyle','none')
plot(time0(range),FlowData(range),'Color',ColorFlow)
xlim([time0(range(1)) time0(range(end))])
xticks([0 5 10])
xlabel('Time [s]')
ylabel('Flow [L/s]')
ylim([y1 y2])
set(gca, 'Layer', 'top')
hold off

% subplot(5,1,5)
% hold on
% y1=0;
% y2=max(newVolAutoData(range)*1000+10);
% [v,f] = createBreathPatchGuidelines(Vbnds0,y1,y2,time0(range(end)));
% patch('Faces',f,'Vertices',v,'FaceColor',[0.95 0.95 0.95],'LineStyle','none')
% plot(time0(range),newVolAutoData(range).*1000,'Color',ColorVol)
% xlim([time0(range(1)) time0(range(end))])
% xlabel('Time [s]')
% ylabel('Tidal Volume [mL]')
% ylim([y1 y2])
% set(gca, 'Layer', 'top')
% hold off

%% Save
pathWithFolderName =  strcat(pwd,'\Figures For Paper\');
Condition = 'ExampleTransDiPressure';
switch FileNum
    case 31 %dev
        Case = 'Device';
    case 32 %spont
        Case = 'Spont';
    case 35 %mech
        Case = 'MechV';
end
figCondition = strcat(Case,Condition);
figPrefix = strcat('FN',num2str(FileNum));

figName = strcat(figPrefix,figCondition);
figFileName = strcat(pathWithFolderName,figName);

exportgraphics(f1,strcat(figFileName,'.eps'),'ContentType','vector') %will save figure f as a .png
savefig(f1,strcat(figFileName,'.fig')) %will save figure f as a .fig

%% Save data to Excel file (NBME request)
writematrix(time0(range)', 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig7',figPrefix), 'Range','A2');
writematrix(norm_filtPPl(range)', 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig7',figPrefix), 'Range','B2');
writematrix(norm_filtPAb(range)', 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig7',figPrefix), 'Range','C2');
writematrix(norm_filtPDi', 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig7',figPrefix), 'Range','D2');
writematrix(FlowData(range)', 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig7',figPrefix), 'Range','E2');

end