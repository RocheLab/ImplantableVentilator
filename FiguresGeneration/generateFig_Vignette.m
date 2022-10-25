%generateFig_Vignette

close all 
clear all

FileNum = 1;


load('listVignettes.mat') %read this file to look at names of files and the associated notes
load('fileBlockInfo_C1.mat')
% % Prefix for figures to be saved as
% figPrefix = 'CoherenceFig_';


%% Load data
filepath = strcat('Vignette Segment Selection Files/',list.name{FileNum});
load(filepath);
experiment = ExpN;
save('targetBlock.mat','experiment','Block')
[date,x,y,num] = selectBlock(Block,experiment,fileBlockInfo);
[blocktimes,com,comtext,data,dataend,datastart,firstsampleoffset,rangemax,rangemin,samplerate,tickrate,titles,unittext,unittextmap,filename] = loadPowerLabFile(date,x,y);
[time,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,VolAbsData,DateTime] = blockAnalysis(data,num,datastart,dataend,blocktimes,titles);
[time,VolAutoData,VolAbsData,newVolAbsData,Correction] = calibrateVolumeRemoveIntegralNoise(time,VolAutoData,VolAbsData);
[newVolAutoData] = spirometryNormalization(time,newVolAbsData);
%Trim data into the selected portion for the file we are loading
[time,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,VolAbsData,newVolAbsData,newVolAutoData] = trimData_exp(timeStart,timeDur,time,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,VolAbsData,newVolAbsData,newVolAutoData);
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

%%
time=(time-time(1));%in s
f=figure('Position',[00 50 400 260]);
plotSelectChannelsShort(EKGON,SpO2ON,PArtON,CapnoON, FlowON,PActON,PPlON, PAbON,VolAutoON,VolAbsON,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,newVolAutoData,newVolAbsData, time, Title)
% plotSelectChannelsShort(EKGON,SpO2ON,PArtON,CapnoON, FlowON,PActON,PPlON, PAbON,VolAutoON,VolAbsON,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,newVolAbsData, timeStart,timeDur, Title)

% plotSelectChannels(Title,EKGON,SpO2ON,PArtON,CapnoON, FlowON,PActON,PPlON, PAbON,VolAutoON,VolAbsON,minVentON,fftPActON,fftFlowON, pHON, pCO2ON,pO2ON,sO2ON, EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,newVolAbsData,TminVent,PminVent,xFreqFlow,yMagFlow,xFreqPAct,yMagPAct, ABGTimes,pHvalues,pCO2values,pO2values,sO2values,timeStart,timeDur)
% end


%%
%% Save screencapture


pathWithFolderName =  strcat(pwd,'\Figures For Paper\');
Condition = 'Vignette';
figCondition = strcat('',Condition);
figPrefix = strcat('FN',num2str(FileNum));

figName = strcat(figPrefix,figCondition);
figFileName = strcat(pathWithFolderName,figName);

savefig(f,strcat(figFileName,'.fig')) %will save figure f as a .fig
exportgraphics(f,strcat(figFileName,'.eps'),'ContentType','vector') %will save figure f as a .png

%% Save data to Excel file (NBME request)



writematrix(time', 'SourceData_nBME-21-2902.xlsx', 'Sheet', 'Fig3a', 'Range','A2');
writematrix(PActData', 'SourceData_nBME-21-2902.xlsx', 'Sheet', 'Fig3a', 'Range','B2');
writematrix(FlowData', 'SourceData_nBME-21-2902.xlsx', 'Sheet', 'Fig3a', 'Range','C2');
writematrix(VolAutoData', 'SourceData_nBME-21-2902.xlsx', 'Sheet', 'Fig3a', 'Range','D2');
