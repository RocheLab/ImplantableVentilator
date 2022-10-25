%generateFig_Alignment_PeakInspandVol_REV


% this is modified from runPeakInspFlowandVol
close all 
clear all



% for FileNum = [9 10 36 37 44 45 46 47 69 55]
FileNum = 50;%which file in the list to load
%45: 6d and 6f and 50: 6c and 6e


removedoubleactuationON = 0;


close all

load('listFullSegments.mat') %read this file to look at names of files and the associated notes
load('fileBlockInfo_C1.mat')
% % Prefix for figures to be saved as
% figPrefix = 'CoherenceFig_';


%% Load data
filepath = strcat('Full Segments Selection Files/',list.name{FileNum});
load(filepath);
experiment = ExpN;
save('targetBlock.mat','experiment','Block')
[date,x,y,num] = selectBlock(Block,experiment,fileBlockInfo);
[blocktimes,com,comtext,data,dataend,datastart,firstsampleoffset,rangemax,rangemin,samplerate,tickrate,titles,unittext,unittextmap,filename] = loadPowerLabFile(date,x,y);
[time,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,VolAbsData,DateTime] = blockAnalysis(data,num,datastart,dataend,blocktimes,titles);
%Trim data into the selected portion for the file we are loading
[time,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,VolAbsData] = trimData(timeStart,timeDur,time,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,VolAbsData);

Condition = list.name{FileNum}; %used for filenaming


%% Use the volume minima to determine the limits of each breath

%lin interp to normalize vol abs data

% figure('Position',[100 50 700 300]);
v_interp = interp1([time(100) time(end-100)],[VolAbsData(100) VolAbsData(end-100)],time);
Vol_interp = VolAbsData-v_interp;
Vminind = islocalmin(Vol_interp,'FlatSelection','last','MinSeparation',600,'MinProminence',0.01);
Vminbounds = find(Vminind);
% hold on
% plot(time,Vol_interp,time(Vminind),Vol_interp(Vminind),'r*')
% xlabel('Time[s]')
[Volpks,Vollocs]=findpeaks(Vol_interp,'MinPeakDistance',0.600,'MinPeakProminence',0.01);
% hold off

PActminind = islocalmin(PActData,'FlatSelection','last','MinSeparation',600,'MinProminence',2);
PActbounds = find(PActminind);
% figure('Position',[100 430 700 300]);
% hold on
% plot(time,PActData,time(PActminind),PActData(PActminind),'r*')
% xlabel('Time[s]')
% findpeaks(PActData,'MinPeakDistance',0.600,'MinPeakProminence',2);
% hold off

%finds pact peak which may not actually be the point we are most interested
%in
% [PActpks,PActlocs]=findpeaks(PActData,'MinPeakDistance',0.600,'MinPeakProminence',2);
%find alt pact peak end 
%(find peak drop in pressure, which corresponds to |10% drop in pressure)
%%%%%% alternate code %%%%%
dPAct = diff(PActData);
PActlocinds = islocalmin(dPAct,'MinSeparation',600,'MinProminence',0.3);
PActlocs = find(PActlocinds);
% plot(time(1:end-1),dPAct,'-',time(PActlocinds),dPAct(PActlocinds),'r*');
removedoubleactuationON = 0;
if removedoubleactuationON == 1
PActlocsdist = PActlocs(2:end)-PActlocs(1:end-1);
PAdind = find(PActlocsdist<900);

    for ii = 1:length(PAdind)
    remove = PAdind(end+1-ii); %go in reverse order
    PActlocs = [PActlocs(1:remove) PActlocs(remove+2:end)];
    end
plot(time(1:end-1),dPAct,'-',time(PActlocinds),dPAct(PActlocinds),'r*',time(1)+0.001*PActlocs,-0.3*ones(size(PActlocs)),'x');
    
end





% figure
% findpeaks(FlowData,time,'MinPeakDistance',0.6,'MinPeakHeight',0.01,'MinPeakProminence',0.1);
[Flowpks,Flowlocs]=findpeaks(FlowData,'MinPeakDistance',0.6,'MinPeakHeight',0.01,'MinPeakProminence',0.1);


% length(find(Vminind))
% length(Volpks)
% length(find(PActminind))
% length(PActpks)
% length(Flowpks)



[Output] = analyzeSnipPFV(FlowData,VolAutoData,Vollocs,PActlocs,Flowlocs,Vminbounds,PActbounds);
%%%                     1       2       3       4      5    6       7     8      9      10   11  12  13      14    15         16     17     
%%% Output content is [bound1 bound2 Fpkloc Fpkval Ppkloc Pminloc Vpkloc Vpkval Vminloc F_P P_V F_V F_Pmin Pmin_V Pmin_Vmin  P_Vmin F_Vmin ];


[Time_Exc] = timeExclude(50,50,Output);


%%%%%%%%%                       Var. Set #                       
% Output(:,11) % P_V %          1
% Output(:,13) % F_Pmin %       2
% Output(:,14) % Pmin_V %       3
% Output(:,15) % Pmin_Vmin %    4  
% Output(:,16) % P_Vmin %       5
% Output(:,17) % F_Vmin %       6
%%%%%
lim1_1 = -800;
lim1_2 = 1000;
lim2_1 = 0;
lim2_2 = 600;
lim3_1 = -1500;
lim3_2 = 000;
lim4_1 = -400;
lim4_2 = 600;
lim5_1 = 300;
lim5_2 = 1500;
lim6_1 = 0; 
lim6_2 = 1500;
[Output_Exc,num_Exc] = exclusionCriteria(Time_Exc,lim1_1,lim1_2,lim2_1,lim2_2,lim3_1,lim3_2,lim4_1,lim4_2,lim5_1,lim5_2,lim6_1,lim6_2);
%%%%%%%%%%%

%% multivariate regression
y1 = Output_Exc(:,8).*1000;
x1 = Output_Exc(:,15);
x2 = -Output_Exc(:,14);


X = [ones(size(x1)) x1 x2 x1.*x2];
[b_vol,bint_v,r_v,rint_v,stats_v] = regress(y1,X);

ff1 = figure
scatter3(x1,x2,y1,'filled')
hold on
x1fit = min(x1):100:max(x1);
x2fit = min(x2):10:max(x2);
[X1FIT_vol,X2FIT_vol] = meshgrid(x1fit,x2fit);
YFIT = b_vol(1) + b_vol(2)*X1FIT_vol + b_vol(3)*X2FIT_vol + b_vol(4)*X1FIT_vol.*X2FIT_vol;
mesh(X1FIT_vol,X2FIT_vol,YFIT)
xlabel('P_0-V_0')
ylabel('P_0-V_p_k')
zlabel('Volume [mL]')
view(50,10)
hold off

y2 = Output_Exc(:,4);
[b_flow,bint_f,r_f,rint_f,stats_f] = regress(y2,X);

ff2 = figure
scatter3(x1,x2,y2,'filled')
hold on
x1fit = min(x1):100:max(x1);
x2fit = min(x2):10:max(x2);
[X1FIT_flow,X2FIT_flow] = meshgrid(x1fit,x2fit);
YFIT = b_flow(1) + b_flow(2)*X1FIT_flow + b_flow(3)*X2FIT_flow + b_flow(4)*X1FIT_flow.*X2FIT_flow;
mesh(X1FIT_flow,X2FIT_flow,YFIT)
xlabel('P_0-V_0')
ylabel('P_0-V_p_k')
zlabel('Flow [L/s]')
view(50,10)
hold off

% [Output_Exc] = timeExclude(50,50,Output);
% [Output_Exc,num_Exc] = exclusionCriteria(Output);



% %hypothesis, F_P (Output_Exc(:,10)) should always be around -420 std 11
% %because it is an inherent property of the system?
% Mean_F_P= mean(Output_Exc(:,10),'omitnan')
% StD_F_P= std(Output_Exc(:,10),'omitnan')



%% figure generation for peak insp flow
close all

sz=8;

f1=figure('Position',[00 50 200 500]);
% subplot(2,4,1)
% scatter(Output_Exc(:,11),Output_Exc(:,4),[],Output_Exc(:,1))
% title({'Effect of peak-peak distance between','volume and PAct on peak insp flow'})
% xlabel('Time between PAct peak-Vol peak [ms]')
% ylabel('Peak Insp Flow [L/s]')
% % if max(Output_Exc(:,11))>1000 && min(Output_Exc(:,11))<-1000
% %     xlim([-1000,1000])
% % elseif max(Output_Exc(:,11))>1000
% %     xlim([-inf,1000])
% % elseif min(Output_Exc(:,11))<-1000
% %     xlim([-1000,inf])
% % end
% 
% % 
% % figure('Position',[500 430 500 300]);
% subplot(2,4,4)
% scatter(Output_Exc(:,10),Output_Exc(:,4),[],Output_Exc(:,1))
% title({'Effect of peak-peak distance between','flow and pressure on peak insp flow'})
% xlabel('Time between Flow peak-PAct peak [ms]')
% ylabel('Peak Insp Flow [L/s]')
% % if max(Output_Exc(:,16))>1000 && min(Output_Exc(:,16))<-1000
% %     xlim([-1000,1000])
% % elseif max(Output_Exc(:,16))>1000
% %     xlim([-inf,1000])
% % elseif min(Output_Exc(:,16))<-1000
% %     xlim([-1000,inf])
% % end
% 
% % figure('Position',[1000 430 500 300]);
% subplot(2,4,3)
% scatter(Output_Exc(:,12),Output_Exc(:,4),[],Output_Exc(:,1))
% title({'Effect of peak-peak distance between','flow and volume on peak insp flow'})
% xlabel('Time between Flow peak-Vol peak [ms]')
% ylabel('Peak Insp Flow [L/s]')
% % if max(Output_Exc(:,17))>1000 && min(Output_Exc(:,17))<-1000
% %     xlim([-1000,1000])
% % elseif max(Output_Exc(:,17))>1000
% %     xlim([-inf,1000])
% % elseif min(Output_Exc(:,17))<-1000
% %     xlim([-1000,inf])
% % end
% 
% 
% 
% % figure('Position',[500 50 500 300]);
% subplot(2,4,8)
% scatter(Output_Exc(:,13),Output_Exc(:,4),[],Output_Exc(:,1))
% title({'Effect of distance between peak flow','and actuation start on peak insp flow'})
% xlabel('Time between Flow peak-PAct start [ms]')
% ylabel('Peak Insp Flow [L/s]')
% % if max(Output_Exc(:,13))>1000 && min(Output_Exc(:,13))<-1000
% %     xlim([-1000,1000])
% % elseif max(Output_Exc(:,13))>1000
% %     xlim([-inf,1000])
% % elseif min(Output_Exc(:,13))<-1000
% %     xlim([-1000,inf])
% % end
% 



% % figure('Position',[0000 430 500 300]);
subplot(2,1,2)
% scatter(Output_Exc(:,15),Output_Exc(:,8),sz,Output_Exc(:,1))
scatter(Output_Exc(:,15),Output_Exc(:,8).*1000,sz,hex2rgb('#D5A427'),'filled')
% title({'Effect of distance between start of breath and' ,'actuation start on volume'})
xlabel('P_0-V_0 [ms]')
ylabel('Tidal Volume [mL]')
ylim([50 220])
% % if max(Output_Exc(:,15))>1000 && min(Output_Exc(:,15))<-1000
% %     xlim([-1000,1000])
% % elseif max(Output_Exc(:,15))>1000
% %     xlim([-inf,1000])
% % elseif min(Output_Exc(:,15))<-1000
% %     xlim([-1000,inf])
% % end

% % figure('Position',[0000 430 500 300]);
subplot(2,1,1)
% scatter(Output_Exc(:,15),Output_Exc(:,4),sz,Output_Exc(:,1))
scatter(Output_Exc(:,15),Output_Exc(:,4),sz,hex2rgb('#273866'),'filled')
% title({'Effect of distance between start of breath and' ,'actuation start on peak insp flow'})
xlabel('P_0-V_0 [ms]')
ylim([0.22 0.41])
ylabel('Peak Inspiratory Flow [L/s]')
% % if max(Output_Exc(:,15))>1000 && min(Output_Exc(:,15))<-1000
% %     xlim([-1000,1000])
% % elseif max(Output_Exc(:,15))>1000
% %     xlim([-inf,1000])
% % elseif min(Output_Exc(:,15))<-1000
% %     xlim([-1000,inf])
% % end
% 
% % figure('Position',[500 430 500 300]);
% subplot(2,4,7)
% scatter(Output_Exc(:,16),Output_Exc(:,4),[],Output_Exc(:,1))
% title({'Effect of distance between beginning of','breath and peak actuation pressure on peak insp flow'})
% xlabel('Time between PAct peak-Vol start [ms]')
% ylabel('Peak Insp Flow [L/s]')
% % if max(Output_Exc(:,16))>1000 && min(Output_Exc(:,16))<-1000
% %     xlim([-1000,1000])
% % elseif max(Output_Exc(:,16))>1000
% %     xlim([-inf,1000])
% % elseif min(Output_Exc(:,16))<-1000
% %     xlim([-1000,inf])
% % end
% % 
% % figure('Position',[1000 430 500 300]);
% subplot(2,4,5)
% scatter(Output_Exc(:,17),Output_Exc(:,4),[],Output_Exc(:,1))
% title({'Effect of distance between beginning','of breath and peak flow on peak insp flow'})
% xlabel('Time between Flow peak-Vol start [ms]')
% ylabel('Peak Insp Flow [L/s]')
% % if max(Output_Exc(:,17))>1000 && min(Output_Exc(:,17))<-1000
% %     xlim([-1000,1000])
% % elseif max(Output_Exc(:,17))>1000
% %     xlim([-inf,1000])
% % elseif min(Output_Exc(:,17))<-1000
% %     xlim([-1000,inf])
% % end
% 
% 
% %%%%%%
% 








%%%% Figure Plotting code for various combinations

% % figure('Position',[500 430 500 300]);
% % scatter(Output(:,11),Output(:,8),[],Output(:,1))
% % title('Effect of peak-peak distance between volume and PAct on tidal volume')
% % xlabel('Time between peak volume and peak actuation pressure [ms]')
% % ylabel('Tidal Volume [L]')

f3 = figure('Position',[1000 50 500 500]);
% scatter(Output_Exc(:,4),Output_Exc(:,8),sz,Output_Exc(:,1))
scatter(Output_Exc(:,4),Output_Exc(:,8).*1000,sz)
title('Effect of peak insp flow on tidal volume')
xlabel({'Peak Inspiratory Flow [L/s]'})
ylabel('Tidal Volume [mL]')

% 
% % figure('Position',[00 430 500 300]);
% % scatter(Output(:,10),Output(:,8),[],Output(:,1))
% % title('Effect of peak-peak distance between flow and PAct on tidal volume')
% % xlabel('Time between peak flow and peak actuation pressure [ms]')
% % ylabel('Tidal Volume [L]')
% % 
% % figure('Position',[00 50 500 300]);
% % scatter(Output(:,10),Output(:,4),[],Output(:,1))
% % title('Effect of peak-peak distance between flow and PAct on peak insp flow')
% % xlabel('Time between peak flow and peak actuation pressure [ms]')
% % ylabel('Peak Insp Flow [L/s]')

%% 
f2=figure('Position',[00 50 200 500]);
% subplot(2,4,1)
% scatter(Output_Exc(:,11),Output_Exc(:,8),[],Output_Exc(:,1))
% title({'Effect of peak-peak distance between','volume and PAct on volume'})
% xlabel('Time between PAct peak-Vol peak [ms]')
% ylabel('Tidal Volume [L]')
% % if max(Output_Exc(:,11))>1000 && min(Output_Exc(:,11))<-1000
% %     xlim([-1000,1000])
% % elseif max(Output_Exc(:,11))>1000
% %     xlim([-inf,1000])
% % elseif min(Output_Exc(:,11))<-1000
% %     xlim([-1000,inf])
% % end

% 
% % figure('Position',[500 430 500 300]);
% subplot(2,4,4)
% scatter(Output_Exc(:,10),Output_Exc(:,8),[],Output_Exc(:,1))
% title({'Effect of peak-peak distance between','flow and pressure on volume'})
% xlabel('Time between Flow peak-PAct peak [ms]')
% ylabel('Tidal Volume [L]')
% % if max(Output_Exc(:,16))>1000 && min(Output_Exc(:,16))<-1000
% %     xlim([-1000,1000])
% % elseif max(Output_Exc(:,16))>1000
% %     xlim([-inf,1000])
% % elseif min(Output_Exc(:,16))<-1000
% %     xlim([-1000,inf])
% % end

% % figure('Position',[1000 430 500 300]);
% subplot(2,4,3)
% scatter(Output_Exc(:,12),Output_Exc(:,8),[],Output_Exc(:,1))
% title({'Effect of peak-peak distance between','flow and volume on volume'})
% xlabel('Time between Flow peak-Vol peak [ms]')
% ylabel('Tidal Volume [L]')
% % if max(Output_Exc(:,17))>1000 && min(Output_Exc(:,17))<-1000
% %     xlim([-1000,1000])
% % elseif max(Output_Exc(:,17))>1000
% %     xlim([-inf,1000])
% % elseif min(Output_Exc(:,17))<-1000
% %     xlim([-1000,inf])
% % end



% % figure('Position',[500 50 500 300]);
% subplot(2,4,8)
% scatter(Output_Exc(:,13),Output_Exc(:,8),[],Output_Exc(:,1))
% title({'Effect of distance between peak flow','and actuation start on volume'})
% xlabel('Time between Flow peak-PAct start [ms]')
% ylabel('Tidal Volume [L]')
% % if max(Output_Exc(:,13))>1000 && min(Output_Exc(:,13))<-1000
% %     xlim([-1000,1000])
% % elseif max(Output_Exc(:,13))>1000
% %     xlim([-inf,1000])
% % elseif min(Output_Exc(:,13))<-1000
% %     xlim([-1000,inf])
% % end

% % figure('Position',[1000 50 500 300]);
subplot(2,1,2)
% scatter(Output_Exc(:,14),Output_Exc(:,8),sz,Output_Exc(:,1))
scatter(-Output_Exc(:,14),Output_Exc(:,8).*1000,sz,hex2rgb('#D5A427'),'filled')
% title({'Effect of distance between peak volume','and actuation start on volume'})
xlabel('V_p_k-P_0 [ms]')
ylabel('Tidal Volume [mL]')
ylim([50 220])
% % if max(Output_Exc(:,14))>1000 && min(Output_Exc(:,14))<-1000
% %     xlim([-1000,1000])
% % elseif max(Output_Exc(:,14))>1000
% %     xlim([-inf,1000])
% % elseif min(Output_Exc(:,14))<-1000
% %     xlim([-1000,inf])
% % end

% % figure('Position',[1000 50 500 300]);
subplot(2,1,1)
% scatter(Output_Exc(:,14),Output_Exc(:,4),sz,Output_Exc(:,1))
scatter(-Output_Exc(:,14),Output_Exc(:,4),sz,hex2rgb('#273866'),'filled')
% title({'Effect of distance between peak volume','and actuation start on peak insp flow'})
xlabel('V_p_k-P_0 [ms]')
ylim([0.22 0.41])
ylabel('Peak Inspiratory Flow [L/s]')
% % if max(Output_Exc(:,14))>1000 && min(Output_Exc(:,14))<-1000
% %     xlim([-1000,1000])
% % elseif max(Output_Exc(:,14))>1000
% %     xlim([-inf,1000])
% % elseif min(Output_Exc(:,14))<-1000
% %     xlim([-1000,inf])
% % end

% % figure('Position',[500 430 500 300]);
% subplot(2,4,5)
% scatter(Output_Exc(:,16),Output_Exc(:,8),[],Output_Exc(:,1))
% title({'Effect of distance between beginning of','breath and peak actuation pressure on volume'})
% xlabel('Time between PAct peak-Vol start [ms]')
% ylabel('Tidal Volume [L]')
% % if max(Output_Exc(:,16))>1000 && min(Output_Exc(:,16))<-1000
% %     xlim([-1000,1000])
% % elseif max(Output_Exc(:,16))>1000
% %     xlim([-inf,1000])
% % elseif min(Output_Exc(:,16))<-1000
% %     xlim([-1000,inf])
% % end
% 
% % figure('Position',[1000 430 500 300]);
% subplot(2,4,7)
% scatter(Output_Exc(:,17),Output_Exc(:,8),[],Output_Exc(:,1))
% title({'Effect of distance between beginning','of breath and peak flow on volume'})
% xlabel('Time between Flow peak-Vol start [ms]')
% ylabel('Tidal Volume [L]')
% % if max(Output_Exc(:,17))>1000 && min(Output_Exc(:,17))<-1000
% %     xlim([-1000,1000])
% % elseif max(Output_Exc(:,17))>1000
% %     xlim([-inf,1000])
% % elseif min(Output_Exc(:,17))<-1000
% %     xlim([-1000,inf])
% % end




%% Save screencapture


pathWithFolderName =  strcat(pwd,'\Figures For Paper\');
Condition = 'SeveredPhrenic';
figCondition = strcat('SCALED_GroupedbyVIndex_NoTimeColor_',Condition);
figPrefix = strcat('FN',num2str(FileNum));

figName = strcat(figPrefix,figCondition);
figFileName = strcat(pathWithFolderName,figName);

savefig(f1,strcat(figFileName,'VolStart.fig')) %will save figure f as a .fig
savefig(f2,strcat(figFileName,'VolPeak.fig')) %will save figure f as a .fig
savefig(f3,strcat(figFileName,'PeakInspFlow.fig')) %will save figure f as a .fig
exportgraphics(f1,strcat(figFileName,'VolStart.eps'),'ContentType','vector')
exportgraphics(f2,strcat(figFileName,'VolPeak.eps'),'ContentType','vector') %will save figure f as a .png
exportgraphics(f3,strcat(figFileName,'FlowVolRelationship.eps'),'ContentType','vector') %will save figure f as a .png
% imageData = screencapture(0,  [0,50,1500,600]);  
% imwrite(imageData,figFileName);  % save the captured image to file

%%
%% Save data to Excel file (NBME request)

if FileNum == 50
writematrix(-Output_Exc(:,14), 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig6c'), 'Range','A2');
writematrix(Output_Exc(:,4), 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig6c'), 'Range','B2');
 writematrix(-Output_Exc(:,14), 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig6e'), 'Range','A2');
 writematrix(Output_Exc(:,8), 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig6e'), 'Range','B2');

else 
writematrix(-Output_Exc(:,14), 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig6d'), 'Range','A2');
writematrix(Output_Exc(:,4), 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig6d'), 'Range','B2');
writematrix(-Output_Exc(:,14), 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig6f'), 'Range','A2');
writematrix(Output_Exc(:,8), 'SourceData_nBME-21-2902.xlsx', 'Sheet', strcat('Fig6f'), 'Range','B2');
end

% clear variables
% end