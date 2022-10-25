%generateFig_CalcWorkPV

close all 
clear all

PVdata = cell(1,6); % {dev P, dev V, spont P, spont V, mech P, mech V}

for FileNum = [2 5];

load('listPVloopfromvignettesample.mat')



% load('listVignettes.mat') %read this file to look at names of files and the associated notes
load('fileBlockInfo_C1.mat')
% % Prefix for figures to be saved as
% figPrefix = 'CoherenceFig_';


%% Load data
filepath = strcat('PV Loop Selection Files/',list.name{FileNum});
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
if FileNum ==2 
cutoff = 20*1000;

range1=[1:cutoff];
range2=[cutoff:length(time)];

ColorDev = hex2rgb('#4c9da9');
ColorMech = hex2rgb('#97cad3');
ColorSpont = hex2rgb('#245767');

f1=figure('Position',[100 100 420 300]);
hold on
plot(norm_filtPPl(range1),newVolAutoData(range1).*1000,'Color',ColorDev)
plot(norm_filtPPl(range2),newVolAutoData(range2).*1000,'Color',ColorSpont)
xlabel('\Delta Pleural Pressure [cmH_2O]')
ylabel('Tidal Volume [mL]')

PVdata{1} = norm_filtPPl(range1); %device loops FN2, P
PVdata{2} = newVolAutoData(range1).*1000; %device loops FN2, V
PVdata{3} = norm_filtPPl(range2); %spont loops FN2, P
PVdata{4} = newVolAutoData(range2).*1000; %spont loops FN2, V


end
% 
% %dev vs mech FN5
if FileNum == 5
cutoff = 20*1000;

range1=[1:cutoff];
range2=[cutoff:length(time)];
% 
% ColorDev = hex2rgb('#7262C1');
% ColorMech = ColorDev.*1.3;
% ColorSpont = ColorDev.*0.6;

% f1=figure('Position',[100 100 300 300]);
hold on
% plot(norm_filtPPl(range1),newVolAutoData(range1).*1000,'Color',ColorDev)
plot(norm_filtPPl(range2),newVolAutoData(range2).*1000,'Color',ColorMech)
xlabel('\Delta Pleural Pressure [cmH_2O]')
ylabel('Tidal Volume [mL]')

PVdata{5} = norm_filtPPl(range2); %mech loops FN2, P
PVdata{6} = newVolAutoData(range2).*1000; %mech loops FN2, V


end
end


%% Work calculation

%regression through mech vent data to generate line for passive chest wall
%compliance

regb1 = 0.001*1000;
% regb2 = 3*1000;
regb2 = 5*1000;
regbounds = find(PVdata{6}(regb1:regb2)<350);

% regbounds = (regb1:regb2);
% figure
X = [ones(length(PVdata{6}(regbounds)),1) PVdata{6}(regbounds)']; %using V as x and P as y bc of integration later
% % X = PVdata{6}(regbounds)';
% b = X\PVdata{5}(regbounds)';
% yCalc = X*b;

% Xsamp = [1 300; 1 420]; %volumes
% Ysamp = Xsamp*b; %pressures

% plot(PVdata{5},PVdata{6},'.', yCalc,PVdata{6}(regbounds),'--k')



mdl = fitlm(PVdata{6}(regbounds)',PVdata{5}(regbounds)');
b=mdl.Coefficients.Estimate;
% b(1)=-4; %artificiall move intercept
yCalc = X*b;

% 
%% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%calculate average work for dev
%find the outer bounds delineated by max and min volumes
[DevMax_pks, DevMax_locs] = findpeaks(PVdata{2},'MinPeakHeight',40,'MinPeakDistance',1400); 
[DevMin_TF] = islocalmin(PVdata{2},'MinProminence',40,'MinSeparation',1400);
DevMin_locs = find(DevMin_TF);
% figure
% plot(PVdata{1},PVdata{2},'-',PVdata{1}(DevMax_locs),PVdata{2}(DevMax_locs),'b*',PVdata{1}(DevMin_locs),PVdata{2}(DevMin_locs),'k*');

bounds = [];%initialize
for i = 1:length(DevMin_locs)
    maxind = find(DevMax_locs>DevMin_locs(i));
    if isempty(maxind) %no more max inds remaining
%         continue
    else 
    nextmax = DevMax_locs(maxind(1));
    bounds(1:2,i) = [DevMin_locs(i);nextmax];
    end
end

% integrate

for dn = 1:size(bounds,2)
%subtract out data from regression
Vsamp = [ones(bounds(2,dn)-bounds(1,dn)+1,1) PVdata{2}(bounds(1,dn):bounds(2,dn))'];
Psamp = Vsamp*b;
% plot(PVdata{5},PVdata{6},'.', yCalc,PVdata{6},'--k',Psamp,Vsamp(:,2),'r*')
P_int = Psamp' - PVdata{1}(bounds(1,dn):bounds(2,dn));

% plot(PVdata{2}(bounds(1,1):bounds(2,1)),P_int)
WorkDev_cmH2OmL = trapz(PVdata{2}(bounds(1,dn):bounds(2,dn)),P_int);%in cmH20*mL
WorkDev_J(dn) = WorkDev_cmH2OmL*(98.0665/(10^6));
WorkDev_JpL(dn) = WorkDev_J(dn)/(PVdata{2}(bounds(2,dn))/1000);%J/L
end

sumWorkDev_J = sum(WorkDev_J(1:end-1));
sumTimeD =  0.001*bounds(1,length(bounds))-0.001*bounds(1,1);

WorkDev_Jpmin = sumWorkDev_J*60/sumTimeD; %in Jpmin
avgWorkDev_J = mean(WorkDev_J);
avgWorkDev_JpL = mean(WorkDev_JpL);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%calculate average work for spont
%find the outer bounds delineated by max and min volumes
[SponMax_pks, SponMax_locs] = findpeaks(PVdata{4},'MinPeakHeight',40,'MinPeakDistance',1400); 
[SponMin_TF] = islocalmin(PVdata{4},'MinProminence',40,'MinSeparation',1400);
SponMin_locs = find(SponMin_TF);
% figure
% plot(PVdata{3},PVdata{4},'-',PVdata{3}(SponMax_locs),PVdata{4}(SponMax_locs),'b*',PVdata{3}(SponMin_locs),PVdata{4}(SponMin_locs),'k*');

Sbounds = [];%initialize
for i = 1:length(SponMin_locs)
    Smaxind = find(SponMax_locs>SponMin_locs(i));
    if isempty(Smaxind) %no more max inds remaining
%         continue
    else
        Snextmax = SponMax_locs(Smaxind(1));
    Sbounds(1:2,i) = [SponMin_locs(i);Snextmax];
    end
end

% integrate

for sn = 1:size(Sbounds,2)
%subtract out data from regression
Vsamp = [ones(Sbounds(2,sn)-Sbounds(1,sn)+1,1) PVdata{4}(Sbounds(1,sn):Sbounds(2,sn))'];
Psamp = Vsamp*b;
% plot(PVdata{5},PVdata{6},'.', yCalc,PVdata{6},'--k',Psamp,Vsamp(:,2),'r*')
P_int = Psamp' - PVdata{3}(Sbounds(1,sn):Sbounds(2,sn));

% plot(PVdata{2}(bounds(1,1):bounds(2,1)),P_int)
WorkSpon_cmH2OmL = trapz(PVdata{4}(Sbounds(1,sn):Sbounds(2,sn)),P_int);%in cmH20*mL
WorkSpon_J(sn) = WorkSpon_cmH2OmL*(98.0665/(10^6));
WorkSpon_JpL(sn) = WorkSpon_J(sn)/(PVdata{4}(Sbounds(2,sn))/1000);%J/L
end

sumWorkSpon_J = sum(WorkSpon_J(1:end-1));
sumTimeS =  0.001*Sbounds(1,length(Sbounds))-0.001*Sbounds(1,1);

WorkSpon_Jpmin = sumWorkSpon_J*60/sumTimeS; %in Jpmin
avgWorkSpon_J = mean(WorkSpon_J);
avgWorkSpon_JpL = mean(WorkSpon_JpL);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%add features to generated figure


%%select a representative loop to plot the integration
rep_dev = 3;
VvertD = [ones(2,1) PVdata{2}(bounds(1:2,rep_dev))'];
PvertD = VvertD*b;
xvertsD = [PvertD(1) PVdata{1}(bounds(1,rep_dev):bounds(2,rep_dev)) PvertD(2)];%pressures
yvertsD = [VvertD(1,2) PVdata{2}(bounds(1,rep_dev):bounds(2,rep_dev)) VvertD(2,2)];%volumes

rep_spon = 2;
VvertS = [ones(2,1) PVdata{4}(bounds(1:2,rep_spon))'];
PvertS = VvertS*b;
xvertsS = [PvertS(1) PVdata{3}(bounds(1,rep_spon):bounds(2,rep_spon)) PvertS(2)];%pressures
yvertsS = [VvertS(1,2) PVdata{4}(bounds(1,rep_spon):bounds(2,rep_spon)) VvertS(2,2)];%volumes



plot(yCalc,PVdata{6}(regbounds),'--k') %linear regression for passive chest wall compliance
p = patch(xvertsD,yvertsD,ColorDev,'EdgeColor','none');
p.EdgeColor = [0 0 0];
p.LineStyle = '- -';
p = patch(xvertsS,yvertsS,ColorSpont,'EdgeColor','none');
p.EdgeColor = [0 0 0];
p.LineStyle = '- -';
alpha(0.3)














%% Save screencapture


pathWithFolderName =  strcat(pwd,'\Figures For Paper\');
Condition = 'PVLoop_WorkCalculation';
figCondition = strcat('',Condition);
%figPrefix = strcat('FN',num2str(FileNum));

figName = strcat(figCondition);
figFileName = strcat(pathWithFolderName,figName);

savefig(f1,strcat(figFileName,'.fig')) %will save figure f as a .fig
exportgraphics(f1,strcat(figFileName,'.eps'),'ContentType','vector') %will save figure f as a .png

%% Save data to Excel file (NBME request)

writematrix(PVdata{1}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', 'Fig7g', 'Range','A2');
writematrix(PVdata{2}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', 'Fig7g', 'Range','B2');
writematrix(PVdata{3}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', 'Fig7g', 'Range','C2');
writematrix(PVdata{4}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', 'Fig7g', 'Range','D2');
writematrix(PVdata{5}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', 'Fig7g', 'Range','E2');
writematrix(PVdata{6}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', 'Fig7g', 'Range','F2');