%generateFig_ABG_revisions
%code for creating the figure showing ABG data for the Apr 2022 revisions
%modified from analyzeWaveformsAndABGs


clear all
close all


% for FileNum = [5:10 48:51 54] %loop through all of the ones we've already set new ABGs for
for FileNum = [5 6 50] % modif
 %FileNum = 5;

load('listFullSegments.mat') %read this file to look at names of files and the associated notes
load('fileBlockInfo_C1.mat')


newAlignABGON = 0; %does the file need the ABG times aligned (1), or have they already been updated?
% prompt = 'Which file # do you want to load?';
% FileNum=input(prompt);
% FileNum =5;


filepath = strcat('Full Segments Selection Files/',list.name{FileNum});
load(filepath);
experiment = ExpN;
save('targetBlock.mat','experiment','Block')
[date,x,y,num] = selectBlock(Block,experiment,fileBlockInfo);
[blocktimes,com,comtext,data,dataend,datastart,firstsampleoffset,rangemax,rangemin,samplerate,tickrate,titles,unittext,unittextmap,filename] = loadPowerLabFile(date,x,y);
[time,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,VolAbsData,PDiData,DateTime] = blockAnalysis(data,num,datastart,dataend,blocktimes,titles);
[time,VolAutoData,VolAbsData,newVolAbsData,Correction] = calibrateVolumeRemoveIntegralNoise(time,VolAutoData,VolAbsData);
[newVolAutoData] = spirometryNormalization(time,newVolAbsData);
blockStart = datetime(blocktimes(num),'ConvertFrom','datenum');
%trim data
% if FileNum==5
% [time,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,VolAbsData,newVolAbsData,newVolAutoData,DateTime] = trimData_exp(timeStart,timeDur-31,time,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,VolAbsData,newVolAbsData,newVolAutoData,DateTime);
% end
%% Switches for plotting different datasets
EKGON = 0;
SpO2ON = 0; %pleth waveform
PArtON = 0;
CapnoON = 0;
FlowON = 1;
PActON = 1;
PPlON = 0;
PAbON = 0;
VolAutoON = 1;
VolAbsON = 0;
minVentON = 0;
fftPActMagON=0;
fftPActPhaseON=0; 
fftFlowMagON=0;
fftFlowPhaseON=0;
coherenceON=0;
fftPActON=0;
fftFlowON=0;

pHON = 1;
pCO2ON = 1;
pO2ON = 0;
sO2ON = 0; %discrete points of data from ABG
plotsON = 1; %for analyzing minvent and abg code

% if minVentON == 1
%     %Calculate minute ventilation over time
%     [TminVent,PminVent,NminVent] = minVentOverTime(FlowData,timeStart,timeDur,20,0,0, 0, 0, 0,0.001);
%     TminVent = DateTime + seconds(TminVent);
% else
%     TminVent=[];
%     PminVent=[];
%     NminVent=[];
% end
    TminVent=[];
    PminVent=[];
    NminVent=[];



%% Trim data loaded into GUI to selected portion
% %Resave blockdata.mat for input into the GUI
% save('blockdata.mat','time','EKGData','SpO2Data','PArtData','CapnoData','FlowData','PActData','PPlData','PAbData','VolAutoData','VolAbsData','DateTime')
% 
% 
% %run GUI to select a smaller wave
% app = viewAndSelectData; %open selection app
% waitfor(app); %wait for app to close
% load('call.mat')
% load(fileOutputName)
% disp(strcat('Notes about this selection: ',notes))

% ABGfilename = [num2str(date) 'ABGData.mat'];
% load(ABGfilename)

if newAlignABGON == 1
    [ABGTimes,pHvalues,pCO2values,pO2values,sO2values] = alignAndPlotABGData(date,DateTime, timeStart, timeDur, time,PArtData);
elseif newAlignABGON == 0
    [ABGTimes,pHvalues,pCO2values,pO2values,sO2values,priorPoint] = callABGData(date,DateTime, timeStart, timeDur);
end



%% Label apnea and maintence sections

apneaVSmaintenence(date,time,PActData,timeStart,ABGTimes)

%% Plotting select channels





if fftPActMagON +fftPActPhaseON+fftFlowMagON+fftFlowPhaseON > 0 
    [xFreqFlow,yDFTFlow,xFreqPAct,yDFTPAct] = fftOfPowerLabData(FlowData,PActData,timeStart,timeDur);
else
    xFreqFlow = [];
    yDFTFlow = [];
    xFreqPAct = [];
    yDFTPAct = [];
    yMagFlow = [];
    yMagPAct = [];
end


Title = list.name{FileNum};


[f, dataSAVE] = plotSelectChannelsABG(Title,EKGON,SpO2ON,PArtON,CapnoON, FlowON,PActON,PPlON, PAbON,VolAutoON,VolAbsON,minVentON,fftPActON,fftFlowON, pHON, pCO2ON,pO2ON,sO2ON, EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,newVolAutoData.*1000,newVolAbsData,TminVent,PminVent,xFreqFlow,yMagFlow,xFreqPAct,yMagPAct, ABGTimes,pHvalues,pCO2values,pO2values,sO2values,priorPoint,timeStart,timeDur,blockStart);

pathWithFolderName =  strcat(pwd,'\Figures For Paper\');
Condition = 'TestABG_0419';
figCondition = strcat('',Condition);
figPrefix = strcat('FN',num2str(FileNum));

figName = strcat(figPrefix,figCondition);
figFileName = strcat(pathWithFolderName,figName);
exportgraphics(f,strcat(figFileName,'.eps'),'ContentType','vector') %will save figure f as a .png
%save('dataSAVE.mat','dataSAVE')

%% Save data to Excel file (NBME request)
writematrix(dataSAVE{1,1}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig5',num2str(FileNum)), 'Range','A2');
writematrix(dataSAVE{1,2}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig5',num2str(FileNum)), 'Range','B2');
writematrix(dataSAVE{2,1}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig5',num2str(FileNum)), 'Range','C2');
writematrix(dataSAVE{2,2}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig5',num2str(FileNum)), 'Range','D2');
writematrix(dataSAVE{3,1}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig5',num2str(FileNum)), 'Range','E2');
writematrix(dataSAVE{3,2}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig5',num2str(FileNum)), 'Range','F2');
writematrix(dataSAVE{4,1}, 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig5',num2str(FileNum)), 'Range','G2');
writematrix(dataSAVE{4,2}, 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig5',num2str(FileNum)), 'Range','H2');
writematrix(dataSAVE{5,1}, 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig5',num2str(FileNum)), 'Range','I2');
writematrix(dataSAVE{5,2}, 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig5',num2str(FileNum)), 'Range','J2');

%% Minute Ventilation and ABG analysis
if  pCO2ON == 1

T = 180;
[baseline,SDbaseline] = calcBaseline(date);

% analyzeMinVentAndABG_newmethod(date,time,newVolAutoData, T,ABGTimes,pHvalues,pCO2values,pO2values,sO2values, baseline, SDbaseline,PActData,DateTime,plotsON) 
    
    
end
%%
 clear all %modif


end %modif

