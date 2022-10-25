%generateFig_VTxPSI_Revision

clear all
close all

f1=figure('Position',[100 100 250 250]);


load('listRevisionCharacterization.mat')

x = [];
y = [];
err = [];
%%%%%%%%%%%%%%%%%
fileset = 2;
%%%%%%%%%%%%%%%%%
switch fileset
    case 1
        experiment = 12;
        load('Block_19_USR_Pressure_Trial.mat')
        Block = 19;
        set = [42 44 46 48 50]; %20s segment
        % set = [43 45 47 49 51]; %10 s segments % for Block 19
    case 2 %%%%%%%%%%%%%%%%%%%%%%%%used this one%%%
        experiment = 12;
        load('Block_16_USL_Shapes_Pressure_Trial.mat')
        Block = 16;
        set = [31 38 36 34 32]; %20 s segments %for Block 16 %spont, 5 10 15 20
        % set = [29 31 33 35 37 ]; %10 s segments % for Block 16
    case 3
        experiment = 10;
        load('Block_34_2022_03_15_Pressurization.mat')
        Block = 34;
        set = [9 3 5 7]; %20 s %spont, 10 15 20
%         set = [8 6 4 10]; %10 s
    case 4 %do not use
        experiment = 14;
        load('Block_22_2022_05_05_Pressurizations.mat')
        Block = 22;
        set = [61 59 57 55 63]; %20 s
%         set = [62 60 58 56 64]; %10 s
    case 5 %do not use
        experiment = 14;
        load('Block_23_2022_05_05_Pressurizations.mat')
        Block = 23;
        set = [73 71 69 67 75];%20s
        set = [74 72 70 68 76];%20s
end


l=1;
for FileNum = set 
%% Load data
filepath = strcat('Revision Characterization Files/',list.name{FileNum});
load(filepath);

[time,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,VolAbsData,DateTime] = blockAnalysis(data,1,datastart,dataend,blocktimes,titles);
[time,VolAutoData,VolAbsData,newVolAbsData,Correction] = calibrateVolumeRemoveIntegralNoise(time,VolAutoData,VolAbsData);
[newVolAutoData] = spirometryNormalization(time,newVolAbsData);
%Trim data into the selected portion for the file we are loading
[time,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,VolAbsData,newVolAbsData,newVolAutoData] = trimData_exp(timeStart,timeDur,time,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,VolAbsData,newVolAbsData,newVolAutoData);
%% Vt analysis
%pull tidal volumes
% figure
% findpeaks(newVolAutoData,'MinPeakDistance',600,'MinPeakHeight',0.008,'MinPeakProminence',0.01);
[Volpks,Vollocs]=findpeaks(newVolAutoData,'MinPeakDistance',600,'MinPeakHeight',0.008,'MinPeakProminence',0.01);

%peak pressure 
% figure
% findpeaks(PActData,'MinPeakDistance',600,'MinPeakHeight',4.5,'MinPeakProminence',0.01);
[Ppks,Plocs]=findpeaks(PActData,'MinPeakDistance',600,'MinPeakHeight',4.5,'MinPeakProminence',0.01);
 
if isempty(Plocs)
    Ppks = zeros(size(Volpks));
    Plocs = zeros(size(Vollocs));
end

%trim any unmatched values
tol = 200; %tolerance in ms
while length(Plocs)~=length(Vollocs)
    if length(Plocs)>length(Vollocs)
        if Plocs(1)-Vollocs(1)>tol || Plocs(1)-Vollocs(1)<tol 
            Plocs = Plocs(2:end);
            Ppks = Ppks(2:end);
        elseif Plocs(end)-Vollocs(end)>tol || Plocs(end)-Vollocs(end)<tol 
            Plocs = Plocs(1:end-1);
            Ppks = Ppks(1:end-1);
        end 
    elseif length(Plocs)<length(Vollocs)
        if Plocs(1)-Vollocs(1)>tol || Plocs(1)-Vollocs(1)<tol 
            Vollocs = Vollocs(2:end);
            Volpks = Volpks(2:end);
        elseif Plocs(end)-Vollocs(end)>tol || Plocs(end)-Vollocs(end)<tol 
            Vollocs = Vollocs(1:end-1);
            Volpks = Volpks(1:end-1);
        end 
    end
end%Vollocs and Plocs should be equal in length now

PVdata{l} = [Ppks;Volpks];
l=l+1;%loop increase
%means and std devs
x = [x mean(Ppks)]; 
y = [y mean(Volpks)];
err = [err std(Volpks.*1000)];

hold on
plot(Ppks,Volpks.*1000,'.','Color',[0.6 0.6 0.6])
xlim([0 21])
xticks([0 5 10 15 20])
xlabel('Actuator Pressure [psi]')
switch fileset
    case {1,2}
ylim([0 100])
    case {3,4,5}
ylim([0 170])
end
ylabel('Tidal Volume [mL]')

end

switch fileset
    case {1}
        [h(1),p(1)] = ttest2(PVdata{3}(2,:), PVdata{4}(2,:));%15,10
        [h(2),p(2)] = ttest2(PVdata{3}(2,:), PVdata{2}(2,:));%15,20
        [h(3),p(3)] = ttest2(PVdata{2}(2,:), PVdata{4}(2,:));%20,10
    case {2}%20 s segments %for Block 16 %spont, 5 10 15 20
        [h(1),p(1)] = ttest2(PVdata{3}(2,:), PVdata{4}(2,:));%15,10
        [h(2),p(2)] = ttest2(PVdata{4}(2,:), PVdata{5}(2,:));%15,20
        [h(3),p(3)] = ttest2(PVdata{3}(2,:), PVdata{5}(2,:));%20,10
        %[31 38 36 34 32]; 
    case {3}%20 s %spont, 10 15 20
        [h(1),p(1)] = ttest2(PVdata{2}(2,:), PVdata{3}(2,:));%10 15
        [h(2),p(2)] = ttest2(PVdata{3}(2,:), PVdata{4}(2,:));%15 20
        [h(3),p(3)] = ttest2(PVdata{2}(2,:), PVdata{4}(2,:));%20 10
end
sigstar({[10 15],[15 20],[10 20]},p)
errorbar(x,y*1000,err,'-_','MarkerSize',8,'MarkerEdgeColor','k','CapSize',4,'Color',[0.4 0.4 0.4])


hold off
%% Save screencapture


pathWithFolderName =  strcat(pwd,'\Figures For Paper\');
Condition = 'TidalVolxPSI_Revision_0506';
figCondition = strcat('',Condition);
figPrefix = strcat('Exp',num2str(experiment),'Block',num2str(Block));

figName = strcat(figPrefix,figCondition);
figFileName = strcat(pathWithFolderName,figName);

% savefig(f1,strcat(figFileName,'.fig')) %will save figure f as a .fig
exportgraphics(f1,strcat(figFileName,'.eps'),'ContentType','vector') %will save figure f as a .png
exportgraphics(f1,strcat(figFileName,'.png'))
