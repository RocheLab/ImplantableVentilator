%generateFig_TidalVol_BarChart_Revision
clear all
close all


% load('Block_16_USL_Shapes_Pressure_Trial.mat')
% load('Block_19_USR_Pressure_Trial.mat')
% 
load('listRevisionCharacterization.mat')

l=1;
x=[];
y=[];
err=[];

        experiment = 12;
        load('Block_14_USR_Shapes_Trial.mat')
        Block = 14;
        set = [11 14 17 15 ]; %curve,sq,tri,spont
        
for FileNum = set; 


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


end

% set = [11 14 17 15 ]; %spont,curve,sq,tri,spont


%% Stats

% set = [11 14 17 15 ]; %spont,curve,sq,tri,spont
        [h(1),p(1)] = ttest2(PVdata{1}(2,:), PVdata{2}(2,:));%
        [h(2),p(2)] = ttest2(PVdata{2}(2,:), PVdata{3}(2,:));%
        [h(3),p(3)] = ttest2(PVdata{3}(2,:), PVdata{4}(2,:));%
        [h(4),p(4)] = ttest2(PVdata{2}(2,:), PVdata{4}(2,:));%
        save('Stats_ED_Fig1_q','p')
        %[31 38 36 34 32]; 

%% Plot
%         Color5 = hex2rgb('#C9DEDB');
%         Color10 = hex2rgb('#83ABA4');
%         Color15 = hex2rgb('#428C80');
        ColorCurve = hex2rgb('#146759');
        ColorTri = hex2rgb('#DDCC77');
        ColorSq = hex2rgb('#CC6677');
        ColorSpont = [0.4 0.4 0.4];
%         ColorMech = [0.8 0.8 0.8];

f1=figure('Position',[100 100 150 250]);
        hold on
% X = categorical({'Spontaneous', 'Curved','Square','Triangle'});
% X = reordercats(X,{'Spontaneous', 'Curved','Square','Triangle'});
X = [1:4];
b = bar(X,y.*1000,'BarWidth',0.45);
b.FaceColor = 'flat';
b.CData(1,:) = ColorSpont;
b.CData(2,:) = ColorCurve;
b.CData(3,:) = ColorSq;
b.CData(4,:) = ColorTri;
ylabel('Tidal Volume [mL]')
xticks([1:4])
xticklabels({'Spontaneous', 'Curved','Square','Triangle'})
    for ii = 1:4
        ydots = PVdata{1,ii}(2,:)*1000; %1 for Ppl, 2 for Pab, 3 for Pdi
        xdots = X(ii).*ones(size(ydots))+0.8.*(rand(size(ydots))-0.5);
        plot(xdots,ydots,'.','Color',[0.7 0.7 0.7])
        %Fig7a{ii,2} = ydots;
    end
errorbar([1:4],y*1000,err,'k_','CapSize',4)

sigstar({[1 2],[2 3],[3 4],[2 4]},p)
% Save
pathWithFolderName =  strcat(pwd,'\Figures For Paper\');
Condition = 'TidalVol_WaveformShapes_BarChart';

%%

figName1 = strcat(Condition);
figFileName1 = strcat(pathWithFolderName,figName1);
exportgraphics(f1,strcat(figFileName1,'.eps'),'ContentType','vector') %will save figure f
