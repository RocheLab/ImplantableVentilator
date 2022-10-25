%generateFig_Alignment_Examples


close all 
clear all



% for FileNum = [9 10 36 37 44 45 46 47 69 55]
FileNum = 1;%which file in the list to load
%1: Fig6b December 21st 2020
%3: Fig6a


close all

load('listAlignmentExamples.mat') %read this file to look at names of files and the associated notes
load('fileBlockInfo_C1.mat')
% % Prefix for figures to be saved as
% figPrefix = 'CoherenceFig_';


%% Load data
filepath = strcat('Alignment Example Segment Files/',list.name{FileNum});
load(filepath);
experiment = ExpN;
save('targetBlock.mat','experiment','Block')
[date,x,y,num] = selectBlock(Block,experiment,fileBlockInfo);
[blocktimes,com,comtext,data,dataend,datastart,firstsampleoffset,rangemax,rangemin,samplerate,tickrate,titles,unittext,unittextmap,filename] = loadPowerLabFile(date,x,y);
[time,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,VolAbsData,DateTime] = blockAnalysis(data,num,datastart,dataend,blocktimes,titles);
[time,VolAutoData,VolAbsData,newVolAbsData,Correction] = calibrateVolumeRemoveIntegralNoise(time,VolAutoData,VolAbsData);
[newVolAutoData] = spirometryNormalization(time,newVolAbsData);
%Trim data into the selected portion for the file we are loading
[time,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,VolAbsData,newVolAbsData,newVolAutoData] = trimData_exp(timeStart,2.5,time,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,VolAbsData,newVolAbsData,newVolAutoData);
%%
%%%%%%%%% Plotting Channels %%%%%%%%

%Set value to 1 to turn that channel plotting ON, set to 0 to turn OFF
EKGON = 0;
SpO2ON = 0;
PArtON = 0;
CapnoON = 0;
FlowON = 1;
PActON = 1;
PPlON = 0;
PAbON = 0;
VolAutoON = 1; 
VolAbsON = 0;
minVentON = 0;
fftPActON = 0;
fftFlowON = 0;
pHON=0;
pCO2ON=0;
pO2ON=0;
sO2ON=0;

Title = '';
TminVent=[];
PminVent=[];
xFreqFlow=[];
yMagFlow=[];
xFreqPAct=[];
yMagPAct=[];
ABGTimes=[];
pHvalues=[];
pCO2values=[];
pO2values=[];
sO2values=[];




time=(time-time(1));%in s


%%%%%%%% find max and mins


f=figure('Position',[100 100 200 200]);
% hold on
[dataSAVE] = plotSelectChannelsMS(EKGON,SpO2ON,PArtON,CapnoON, FlowON,PActON,PPlON, PAbON,VolAutoON,VolAbsON,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,newVolAutoData,newVolAbsData, time, Title)
% plotSelectChannelsShort(EKGON,SpO2ON,PArtON,CapnoON, FlowON,PActON,PPlON, PAbON,VolAutoON,VolAbsON,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,newVolAbsData, timeStart,timeDur, Title)

% plotSelectChannels(Title,EKGON,SpO2ON,PArtON,CapnoON, FlowON,PActON,PPlON, PAbON,VolAutoON,VolAbsON,minVentON,fftPActON,fftFlowON, pHON, pCO2ON,pO2ON,sO2ON, EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,newVolAbsData,TminVent,PminVent,xFreqFlow,yMagFlow,xFreqPAct,yMagPAct, ABGTimes,pHvalues,pCO2values,pO2values,sO2values,timeStart,timeDur)
% end
% 
% %%%% plot max and mins
% subplot(3,1,1)
% plot(time(Pmin_ind).*1000,PActData(Pmin_ind),'o')
% subplot(3,1,2)
% plot(time(F_I)*1000,F_M,'o')
% subplot(3,1,3)
% plot(time(V_I)*1000,V_M,'o',time(Vmin_ind).*1000,newVolAutoData(Vmin_ind),'o')
% hold off

%%
%% Save screencapture


pathWithFolderName =  strcat(pwd,'\Figures For Paper\');
Condition = 'SeveredPhrenic_AlignmentExample';
figCondition = strcat('',Condition);
figPrefix = strcat('FN',num2str(FileNum));

figName = strcat(figPrefix,figCondition);
figFileName = strcat(pathWithFolderName,figName);

savefig(f,strcat(figFileName,'.fig')) %will save figure f as a .fig
exportgraphics(f,strcat(figFileName,'.eps'),'ContentType','vector') %will save figure f as a .png

%% Save data to Excel file (NBME request)
writematrix(dataSAVE{1,1}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig6',num2str(FileNum)), 'Range','A2');
writematrix(dataSAVE{1,2}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig6',num2str(FileNum)), 'Range','B2');
writematrix(dataSAVE{1,3}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig6',num2str(FileNum)), 'Range','C2');
writematrix(dataSAVE{1,4}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig6',num2str(FileNum)), 'Range','D2');
writematrix(dataSAVE{2,1}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig6',num2str(FileNum)), 'Range','E2');
writematrix(dataSAVE{2,2}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig6',num2str(FileNum)), 'Range','F2');
writematrix(dataSAVE{2,3}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig6',num2str(FileNum)), 'Range','G2');
writematrix(dataSAVE{2,4}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig6',num2str(FileNum)), 'Range','H2');
writematrix(dataSAVE{3,1}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig6',num2str(FileNum)), 'Range','I2');
writematrix(dataSAVE{3,2}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig6',num2str(FileNum)), 'Range','J2');
writematrix(dataSAVE{3,3}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig6',num2str(FileNum)), 'Range','K2');
writematrix(dataSAVE{3,4}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig6',num2str(FileNum)), 'Range','L2');
writematrix(dataSAVE{3,5}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig6',num2str(FileNum)), 'Range','M2');
writematrix(dataSAVE{3,6}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig6',num2str(FileNum)), 'Range','N2');



