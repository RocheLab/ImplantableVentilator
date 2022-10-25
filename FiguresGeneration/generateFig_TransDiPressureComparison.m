%generateFig_TransDiPressureComparison


%compareTransDiPressures

%data pulled from test_transdiaphragmOverTime
clear all
close all

datadotsON=1;
significanceON =1;
sigbars = 3;
% load('ImmChangePDidata.mat')
% load('listImmediateChangeSegments.mat') %read this file to look at names of files and the associated notes
% load('deltaComparisonVFdata') 

load('listPressureSelection.mat')
load('PressureSelectionsPDidata.mat')


%transDiPressuredata is [(date) (tPl) (delPPl) (tAb) (delAb) (tDi) (delPDi)]


% indep = [1 6 13 19 21 23 27]; %array of independent file numbers

PPl = [];
PAb = [];
PDi = [];
n_Un = 0;
n_Dev = 0;
n_Mech = 0;
n_Sev = 0;
n_Ap = 0;
Unsup = [];
Dev = [];
Mech = [];
Severed = [];
nonSev = [];
Ap = [];
nonAp = [];


% filerange = [1:2 4:10 12:13 14:29 ];
% filerange = [1 2 5 6 12 13 19 20 21 22 23 24 27 28];
% filerange = [1:74];
% filerange = [1:9 29:74];
% filerange = [1:8  33 34 36 37 51:54 69:72]; %same trial segments as
% % DeltaBeforeandAfter
% filerange = [1:41 61:173];

% filerange=[7:10,29:32,76,77,83,84,122:125,163:166]; %old list matched to new list numbers
% filerange = [9 10 12 31 32 35 76 77 79 124 125 127 165 166 168];%noapnea
% filerange = [ 31 32 35 76 77 79 124 125 127 165 166 168];%noapnea %Exclude exp 2
filerange = [ 31 32 35 76 77 79 124 125 127 165 166 168];%noapnea %Exclude exp 2

N = 1;
for i = filerange

%transDiPressuredata is [(date) (tPl) (delPPl) (tAb) (delAb) (tDi) (delPDi)]
    if i == 8
        
        Pl_ind=find(transDiPressuredata{i,2}>0*1000);%find the first index after 0 s
        Ab_ind=find(transDiPressuredata{i,4}>0*1000);%find the first index after 0 s
        Di_ind=find(transDiPressuredata{i,6}>0*1000);%find the first index after 0 s
    else
        Pl_ind=find(transDiPressuredata{i,2}>0*1000);%find the first index after 0 s
        Ab_ind=find(transDiPressuredata{i,4}>0*1000);%find the first index after 0 s
        Di_ind=find(transDiPressuredata{i,6}>0*1000);%find the first index after 0 s
    end
Pstats(i,1) = mean(transDiPressuredata{i,3}(Pl_ind(1):end)); %PPl mean
Pstats(i,2) = std(transDiPressuredata{i,3}(Pl_ind(1):end));%PPl std
Pstats(i,3) = mean(transDiPressuredata{i,5}(Ab_ind(1):end)); %PAb mean
Pstats(i,4) = std(transDiPressuredata{i,5}(Ab_ind(1):end)); %PAb std
Pstats(i,5) = mean(transDiPressuredata{i,7}(Di_ind(1):end)); %PDi mean
Pstats(i,6) = std(transDiPressuredata{i,7}(Di_ind(1):end)); %PDi std
Pstats(i,7) = str2num(list.name{i}(5));%Experiment number

% Pstats = [(PPl mean) (PPl std) (PAb mean) (PAb std) (PDi mean) (PDi std)]

Pdata{N,1} = transDiPressuredata{i,3}(Pl_ind(1):end); %Ppl data
Pdata{N,2} = transDiPressuredata{i,5}(Ab_ind(1):end); %PAb data
Pdata{N,3} = transDiPressuredata{i,7}(Di_ind(1):end); %Pdi data
Pdata{N,4} = i;%FileNumber

%Pdata = {[Ppl] [PAb] [PDi] Exp N}

%PPl/PAb/PDi = [(data) (file number) (0=Unsup,1=Sup) (0=Severed,1=Intact)
%(0=not Apnea,1=Apnea) ] 
if contains(list.Var2{i},'Apn')
    Apnea = 1;
    n_Ap = n_Ap+1;
    Ap = [Ap i];
% elseif contains(list.Var2{i},'Mech') && contains(list.Var2{i},'Bef') %plots baseline mech vent into apnea side
%     Apnea = 1;
%     n_Ap = n_Ap+1;
%     Ap = [Ap i];
else
    Apnea = 0;
    nonAp = [nonAp i];
end


    if contains(list.Var2{i},'Un') || contains(list.Var2{i},'un') %true if unsupported
        
        if isempty(transDiPressuredata{i,8}) %Intact
            PPl = [PPl ; transDiPressuredata{i,3}(Pl_ind(1):end)' i.*ones(size(transDiPressuredata{i,3}(Pl_ind(1):end)'))...
                zeros(size(transDiPressuredata{i,3}(Pl_ind(1):end)')) ones(size(transDiPressuredata{i,3}(Pl_ind(1):end)')) ...
                Apnea.*ones(size(transDiPressuredata{i,3}(Pl_ind(1):end)'))];
            PAb = [PAb ; transDiPressuredata{i,5}(Ab_ind(1):end)' i.*ones(size(transDiPressuredata{i,5}(Ab_ind(1):end)'))...
                zeros(size(transDiPressuredata{i,5}(Ab_ind(1):end)')) ones(size(transDiPressuredata{i,5}(Ab_ind(1):end)')) ...
                Apnea.*ones(size(transDiPressuredata{i,5}(Ab_ind(1):end)'))];
            PDi = [PDi ; transDiPressuredata{i,7}(Di_ind(1):end)' i.*ones(size(transDiPressuredata{i,7}(Di_ind(1):end)'))...
                zeros(size(transDiPressuredata{i,7}(Di_ind(1):end)')) ones(size(transDiPressuredata{i,7}(Di_ind(1):end)')) ...
                Apnea.*ones(size(transDiPressuredata{i,7}(Di_ind(1):end)'))];
            nonSev = [nonSev i];
        elseif contains(transDiPressuredata(i,8),'Severed') %Severed
            PPl = [PPl ; transDiPressuredata{i,3}(Pl_ind(1):end)' i.*ones(size(transDiPressuredata{i,3}(Pl_ind(1):end)'))...
                zeros(size(transDiPressuredata{i,3}(Pl_ind(1):end)')) zeros(size(transDiPressuredata{i,3}(Pl_ind(1):end)')) ...
                Apnea.*ones(size(transDiPressuredata{i,3}(Pl_ind(1):end)'))];
            PAb = [PAb ; transDiPressuredata{i,5}(Ab_ind(1):end)' i.*ones(size(transDiPressuredata{i,5}(Ab_ind(1):end)'))...
                zeros(size(transDiPressuredata{i,5}(Ab_ind(1):end)')) zeros(size(transDiPressuredata{i,5}(Ab_ind(1):end)')) ...
                Apnea.*ones(size(transDiPressuredata{i,5}(Ab_ind(1):end)'))];
            PDi = [PDi ; transDiPressuredata{i,7}(Di_ind(1):end)' i.*ones(size(transDiPressuredata{i,7}(Di_ind(1):end)'))...
                zeros(size(transDiPressuredata{i,7}(Di_ind(1):end)')) zeros(size(transDiPressuredata{i,7}(Di_ind(1):end)')) ...
                Apnea.*ones(size(transDiPressuredata{i,7}(Di_ind(1):end)'))];
            n_Sev = n_Sev+1;
            Severed = [Severed i];
        end
        n_Un=n_Un+1;  
        Unsup = [Unsup i];
    elseif contains(list.Var2{i},'Mech')
        if isempty(transDiPressuredata{i,8}) %Intact
            PPl = [PPl ; transDiPressuredata{i,3}(Pl_ind(1):end)' i.*ones(size(transDiPressuredata{i,3}(Pl_ind(1):end)'))...
                ones(size(transDiPressuredata{i,3}(Pl_ind(1):end)')) ones(size(transDiPressuredata{i,3}(Pl_ind(1):end)')) Apnea.*ones(size(transDiPressuredata{i,3}(Pl_ind(1):end)'))];
            PAb = [PAb ; transDiPressuredata{i,5}(Ab_ind(1):end)' i.*ones(size(transDiPressuredata{i,5}(Ab_ind(1):end)'))...
                ones(size(transDiPressuredata{i,5}(Ab_ind(1):end)')) ones(size(transDiPressuredata{i,5}(Ab_ind(1):end)')) Apnea.*ones(size(transDiPressuredata{i,5}(Ab_ind(1):end)'))];
            PDi = [PDi ; transDiPressuredata{i,7}(Di_ind(1):end)' i.*ones(size(transDiPressuredata{i,7}(Di_ind(1):end)'))...
                ones(size(transDiPressuredata{i,7}(Di_ind(1):end)')) ones(size(transDiPressuredata{i,7}(Di_ind(1):end)')) Apnea.*ones(size(transDiPressuredata{i,7}(Di_ind(1):end)'))];
            nonSev = [nonSev i];
        elseif contains(transDiPressuredata{i,8},'Severed') %Severed
            PPl = [PPl ; transDiPressuredata{i,3}(Pl_ind(1):end)' i.*ones(size(transDiPressuredata{i,3}(Pl_ind(1):end)'))...
                ones(size(transDiPressuredata{i,3}(Pl_ind(1):end)')) zeros(size(transDiPressuredata{i,3}(Pl_ind(1):end)')) Apnea.*ones(size(transDiPressuredata{i,3}(Pl_ind(1):end)'))];
            PAb = [PAb ; transDiPressuredata{i,5}(Ab_ind(1):end)' i.*ones(size(transDiPressuredata{i,5}(Ab_ind(1):end)'))...
                ones(size(transDiPressuredata{i,5}(Ab_ind(1):end)')) zeros(size(transDiPressuredata{i,5}(Ab_ind(1):end)')) Apnea.*ones(size(transDiPressuredata{i,5}(Ab_ind(1):end)'))];
            PDi = [PDi ; transDiPressuredata{i,7}(Di_ind(1):end)' i.*ones(size(transDiPressuredata{i,7}(Di_ind(1):end)'))...
                ones(size(transDiPressuredata{i,7}(Di_ind(1):end)')) zeros(size(transDiPressuredata{i,7}(Di_ind(1):end)')) Apnea.*ones(size(transDiPressuredata{i,7}(Di_ind(1):end)'))];
 
            n_Sev = n_Sev+1;
            Severed = [Severed i];
        end
        n_Mech=n_Mech+1;  
        Mech = [Mech i];
    else %true if device supported
        if isempty(transDiPressuredata{i,8}) %Intact
            PPl = [PPl ; transDiPressuredata{i,3}(Pl_ind(1):end)' i.*ones(size(transDiPressuredata{i,3}(Pl_ind(1):end)'))...
                ones(size(transDiPressuredata{i,3}(Pl_ind(1):end)')) ones(size(transDiPressuredata{i,3}(Pl_ind(1):end)')) Apnea.*ones(size(transDiPressuredata{i,3}(Pl_ind(1):end)'))];
            PAb = [PAb ; transDiPressuredata{i,5}(Ab_ind(1):end)' i.*ones(size(transDiPressuredata{i,5}(Ab_ind(1):end)'))...
                ones(size(transDiPressuredata{i,5}(Ab_ind(1):end)')) ones(size(transDiPressuredata{i,5}(Ab_ind(1):end)')) Apnea.*ones(size(transDiPressuredata{i,5}(Ab_ind(1):end)'))];
            PDi = [PDi ; transDiPressuredata{i,7}(Di_ind(1):end)' i.*ones(size(transDiPressuredata{i,7}(Di_ind(1):end)'))...
                ones(size(transDiPressuredata{i,7}(Di_ind(1):end)')) ones(size(transDiPressuredata{i,7}(Di_ind(1):end)')) Apnea.*ones(size(transDiPressuredata{i,7}(Di_ind(1):end)'))];
            nonSev = [nonSev i];
        elseif contains(transDiPressuredata{i,8},'Severed') %Severed
            PPl = [PPl ; transDiPressuredata{i,3}(Pl_ind(1):end)' i.*ones(size(transDiPressuredata{i,3}(Pl_ind(1):end)'))...
                ones(size(transDiPressuredata{i,3}(Pl_ind(1):end)')) zeros(size(transDiPressuredata{i,3}(Pl_ind(1):end)')) Apnea.*ones(size(transDiPressuredata{i,3}(Pl_ind(1):end)'))];
            PAb = [PAb ; transDiPressuredata{i,5}(Ab_ind(1):end)' i.*ones(size(transDiPressuredata{i,5}(Ab_ind(1):end)'))...
                ones(size(transDiPressuredata{i,5}(Ab_ind(1):end)')) zeros(size(transDiPressuredata{i,5}(Ab_ind(1):end)')) Apnea.*ones(size(transDiPressuredata{i,5}(Ab_ind(1):end)'))];
            PDi = [PDi ; transDiPressuredata{i,7}(Di_ind(1):end)' i.*ones(size(transDiPressuredata{i,7}(Di_ind(1):end)'))...
                ones(size(transDiPressuredata{i,7}(Di_ind(1):end)')) zeros(size(transDiPressuredata{i,7}(Di_ind(1):end)')) Apnea.*ones(size(transDiPressuredata{i,7}(Di_ind(1):end)'))];
 
            n_Sev = n_Sev+1;
            Severed = [Severed i];
        end
        n_Dev=n_Dev+1;  
        Dev = [Dev i];
    end
 
    
    N = N+1;
end



%% Preliminary Plotting

%PPl/PAb/PDi = [(data) (file number) (0=Unsup,1=Sup) (0=Severed,1=Intact) (0=not Apnea,1=Apnea)]

UnsupIndPPl = find(PPl(:,3)==0);
UnsupIndPAb = find(PAb(:,3)==0);
UnsupIndPDi = find(PDi(:,3)==0);
SupIndPPl = find(PPl(:,3)==1);
SupIndPAb = find(PAb(:,3)==1);
SupIndPDi = find(PDi(:,3)==1);

SeverIndPPl = find(PPl(:,4)==0);
SeverIndPAb = find(PAb(:,4)==0);
SeverIndPDi = find(PDi(:,4)==0);
IntactIndPPl = find(PPl(:,4)==1);
IntactIndPAb = find(PAb(:,4)==1);
IntactIndPDi = find(PDi(:,4)==1);

ColorPl = hex2rgb('#7262C1');
ColorAb = hex2rgb('#F15FA6');
ColorDi = hex2rgb('#E87C3B');
ColorGuidePl = ones(n_Un,1).*ColorPl;
ColorGuideAb = ones(n_Mech,1).*ColorAb;

% Pstats = [(PPl mean) (PPl std) (PAb mean) (PAb std) (PDi mean) (PDi std)]

figure
hold on
b1 = bar(filerange,Pstats(filerange,1));
e1=errorbar(filerange,Pstats(filerange,1),Pstats(filerange,2),'.');
e1.Color=ColorPl.*0.8;
% e1M.Color=ColorMech.*0.5;
% e1D.Color=ColorDev.*0.5;
% plot(4.5.*ones(1,2),[0 200],'k',9.5.*ones(1,2),[0 200],'k',28.5.*ones(1,2),...
%     [0 200],'k',38.5.*ones(1,2),[0 200],'k',58.5.*ones(1,2),[0 200],'k') %lines between animals
plot(Ap,250.*ones(size(Ap)),'*')%indicate apnea points
hold off
b1.FaceColor = 'flat';
b1.CData(Unsup,:) = ColorGuidePl;
b1.CData(Mech,:) = ColorGuideAb;
b1.CData(Severed,:) = b1.CData(Severed,:).*0.6;
xlabel('File Number')
ylabel('Change in Pleural Pressure [cmH_2O]')

% figure
% b2 = boxchart(vol(:,2),vol(:,1),'GroupByColor',vol(:,3));
% xlabel('File Number')
% ylabel('Tidal Volume [mL]')
% 
% figure
% swarmchart(vol(UnsupIndV,2),vol(UnsupIndV,1),10,[0.5 0 0.5],'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
% hold on
% swarmchart(vol(SupIndV,2),vol(SupIndV,1),10,[0 0.5 0.5],'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5 );
% hold off
% xlabel('File Number')
% ylabel('Tidal Volume [mL]')

figure
hold on
b2 = bar(filerange,Pstats(filerange,3));
errorbar(filerange,Pstats(filerange,3),Pstats(filerange,4),'.')
% plot(4.5.*ones(1,2),[0 0.7],'k',9.5.*ones(1,2),[0 0.7],'k',28.5.*ones(1,2),[0 0.7],'k',38.5.*ones(1,2),[0 0.7],'k',58.5.*ones(1,2),[0 0.7],'k')
% plot(Ap,0.75.*ones(size(Ap)),'*')%indicate apnea points
hold off
b2.FaceColor = 'flat';
b2.CData(Unsup,:) = ones(length(Unsup),1)*ColorPl;
b2.CData(Mech,:) = ColorGuideAb;
b2.CData(Severed,:) = b2.CData(Severed,:).*0.6;
xlabel('File Number')
ylabel('Change in Abdominal Pressure [cmH_2O]')
% 
% figure
% b5 = boxchart(flow(:,2),flow(:,1),'GroupByColor',flow(:,3));
% xlabel('File Number')
% ylabel('Peak Inspiratory Flow [L/s]')
% 
% figure
% swarmchart(flow(UnsupIndF,2),flow(UnsupIndF,1),10,[0.5 0 0.5],'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
% hold on
% swarmchart(flow(SupIndF,2),flow(SupIndF,1),10,[0 0.5 0.5],'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
% hold off
% xlabel('File Number')
% ylabel('Peak Inspiratory Flow [L/s]')



%% Plot into 4 figures
range1 = intersect(Ap,nonSev);
range2 = intersect(nonAp,nonSev);
% range3 = intersect(Ap,Severed);
% range4 = intersect(nonAp,Severed);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
files = Pstats(:,7);
  order = [0,0,1,0,0,4,3,2]; %from generateFig_DeltaBeforeandAfter when excluding exp 2 and 5
% order = [0,6,1,0,3,5,4];%from generateFig_DeltaBeforeandAfter
  %Exp     1 2 3 4 5 6 7 8
% order = [0 X X 0 X X X X];%order (will dictate which experiment is 1st, 2nd, 3rd,... on the graph)
% const1 = ones(1,length(range1));
% const1(1:2:end-1) = 2; 
% xind1 = 3.*order(filerange(range1))-const1;
const2 = ones(1,length(range2));
const2(1:3:end-2) = 2; %dev
% const2(2:3:end-1) = 1; %spont
const2(3:3:end) = 3; %mech
xind2 = 4.*order(files(range2))-const2;

Mind = (1:3:length(range2)-2); 
Dind = (2:3:length(range2)-1);
Uind = (3:3:length(range2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
f3=figure('Position',[00 635 210 180]);
% subplot(2,2,1)
% hold on
% b3 = bar(range1,Pstats(range1,1));
% errorbar(range1,Pstats(range1,1),Pstats(range1,2),'.')
% plot(20.5.*ones(1,2),[0 -2],'k',41.5.*ones(1,2),[0 -2],'k',60.5.*ones(1,2),...
%     [0 -2],'k',87.5.*ones(1,2),[0 -2],'k',137.5.*ones(1,2),[0 -2],'k') %lines between animals
% b3.FaceColor = 'flat';
% b3.CData(ismember(range1,Unsup),:) = ones(length(intersect(Unsup,range1)),1).*ColorUn;
% b3.CData(ismember(range1,Mech),:) = ones(length(intersect(Mech,range1)),1).*ColorMech;
% ylabel('Change in Pleural Pressure [cmH_2O]')
% hold off


% subplot(2,2,2)
hold on
b4 = bar(xind2,Pstats(range2,1));
%%% code to add dots onto plots
if datadotsON ==1;
    for ii = 1:size(Pdata,1)
        ydots = Pdata{ii,1}; %1 for Ppl, 2 for Pab, 3 for Pdi
        xdots = xind2(ii).*ones(size(ydots))+0.8.*(rand(size(ydots))-0.5);
        plot(xdots,ydots,'.','Color',[0.7 0.7 0.7])
        Fig7a{ii,2} = ydots;

        %         scatter(xdots,ydots,[],'k.','MarkerAlpha',0.2)
    end
end
PstatsexportationFig7a = Pstats(range2,1);
SD7a = Pstats(range2,4);
for i=1:size(xind2,2)
Fig7a{i,1} = xind2(i);
Fig7a{i,3} = PstatsexportationFig7a(i);
Fig7a{i,4} = mean(Fig7a{i,2}); % check if column 3 and 4 are equal
Fig7a{i,5} = SD7a(i);

end
Fig7asorted = sortrows(Fig7a,1);


%%%
%%% code to add on significance bars
% groups={[Mind(1) Dind(1)] [Dind(1) Uind(1)] [Mind(2) Dind(2)] [Dind(2) Uind(2)]...
%     [Mind(3) Dind(3)] [Dind(3) Uind(3)] [Mind(4) Dind(4)] [Dind(4) Uind(4)]};
if significanceON == 1
    file_list = cell2mat(Pdata(:,4));
    p=[];
    for m = 1:4
        DevIndex = find(file_list==Dev(m));
        UnsupIndex = find(file_list==Unsup(m));
        MechIndex = find(file_list==Mech(m));
        if sigbars == 2
       %for only 2 sets of sig bars, Mech v Dev, and Unsup v Dev
        p(2*m-1)=ranksum(Pdata{MechIndex,1},Pdata{DevIndex,1});% %1 for Ppl, 2 for Pab, 3 for Pdi 
        p(2*m)=ranksum(Pdata{UnsupIndex,1},Pdata{DevIndex,1});% %1 for Ppl, 2 for Pab, 3 for Pdi
        groups{2*m-1,1}=[xind2(MechIndex) xind2(DevIndex)];
        groups{2*m,1}=[xind2(UnsupIndex) xind2(DevIndex)];
        elseif sigbars == 3
       %for 3 sets of sig bars, Mech v Dev, and Unsup v Dev
        p(3*m-2)=ranksum(Pdata{MechIndex,1},Pdata{UnsupIndex,1});%
        p(3*m-1)=ranksum(Pdata{MechIndex,1},Pdata{DevIndex,1});% 
        p(3*m)=ranksum(Pdata{UnsupIndex,1},Pdata{DevIndex,1});% 
        groups{3*m-2,1}=[xind2(MechIndex) xind2(UnsupIndex)];
        groups{3*m-1,1}=[xind2(MechIndex) xind2(DevIndex)];
        groups{3*m,1}=[xind2(UnsupIndex) xind2(DevIndex)];
     
        end

    end
end
sigstar(groups,p);
save('Stats_Fig7a','groups','p');
%%%
errorbar(xind2,Pstats(range2,1),Pstats(range2,4),'.k')
% e4U = errorbar(xind2(Uind),Pstats(range2(Uind),1),Pstats(range2(Uind),2),'.');
% e4U.Color = ColorPl*0.7;
% e4M = errorbar(xind2(Mind),Pstats(range2(Mind),1),Pstats(range2(Mind),2),'.');
% e4M.Color = ColorAb*0.7;
% e4D = errorbar(xind2(Dind),Pstats(range2(Dind),1),Pstats(range2(Dind),2),'.');
% e4D.Color = ColorDi*0.7;
% plot(20.5.*ones(1,2),[0 -2],'k',41.5.*ones(1,2),[0 -2],'k',60.5.*ones(1,2),...
%     [0 -2],'k',87.5.*ones(1,2),[0 -2],'k',137.5.*ones(1,2),[0 -2],'k') %lines between animals
b4.FaceColor = 'flat';
b4.CData = ones(length(b4.CData),1).*ColorPl;
b4.CData(Uind,:) = ones(length(intersect(Unsup,range2)),1).*ColorPl*0.6;
b4.CData(Mind,:) = ones(length(intersect(Mech,range2)),1).*ColorPl*1.8;

ylabel('\Delta Pressure_P_l [cmH_2O]')
% set(gca,'xtick',[])
xticks([0:4:20])
xticklabels({})
hold off

% subplot(2,2,3)
% hold on
% b5 = bar(range3,Pstats(range3,1));
% errorbar(range3,Pstats(range3,1),Pstats(range3,2),'.')
% plot(20.5.*ones(1,2),[0 -2],'k',41.5.*ones(1,2),[0 -2],'k',60.5.*ones(1,2),...
%     [0 -2],'k',87.5.*ones(1,2),[0 -2],'k',137.5.*ones(1,2),[0 -2],'k') %lines between animals
% b5.FaceColor = 'flat';
% b5.CData(ismember(range3,Unsup),:) = ones(length(intersect(Unsup,range3)),1).*ColorUn;
% b5.CData(ismember(range3,Mech),:) = ones(length(intersect(Mech,range3)),1).*ColorMech;
% ylabel('Change in Pleural Pressure [cmH_2O]')
% hold off
% 
% subplot(2,2,4)
% hold on
% b6 = bar(range4,Pstats(range4,1));
% errorbar(range4,Pstats(range4,1),Pstats(range4,2),'.')
% plot(20.5.*ones(1,2),[0 -2],'k',41.5.*ones(1,2),[0 -2],'k',60.5.*ones(1,2),...
%     [0 -2],'k',87.5.*ones(1,2),[0 -2],'k',137.5.*ones(1,2),[0 -2],'k') %lines between animals
% b6.FaceColor = 'flat';
% b6.CData(ismember(range4,Unsup),:) = ones(length(intersect(Unsup,range4)),1).*ColorUn;
% b6.CData(ismember(range4,Mech),:) = ones(length(intersect(Mech,range4)),1).*ColorMech;
% ylabel('Change in Pleural Pressure [cmH_2O]')
% hold off


f4=figure('Position',[00 335 210 180]);
% subplot(2,2,1)
% hold on
% b7 = bar(range1,Pstats(range1,3));
% errorbar(range1,Pstats(range1,3),Pstats(range1,4),'.')
% plot(20.5.*ones(1,2),[0 2],'k',41.5.*ones(1,2),[0 2],'k',60.5.*ones(1,2),...
%     [0 2],'k',87.5.*ones(1,2),[0 2],'k',137.5.*ones(1,2),[0 2],'k') %lines between animals
% b7.FaceColor = 'flat';
% b7.CData(ismember(range1,Unsup),:) = ones(length(intersect(Unsup,range1)),1).*ColorUn;
% b7.CData(ismember(range1,Mech),:) = ones(length(intersect(Mech,range1)),1).*ColorMech;
% ylabel('Change in Abdominal Pressure [cmH_2O]')
% hold off

% subplot(2,2,2)
hold on
b8 = bar(xind2,Pstats(range2,3));
if datadotsON ==1;
    for ii = 1:size(Pdata,1)
        ydots = Pdata{ii,2}; %1 for Ppl, 2 for Pab, 3 for Pdi
        xdots = xind2(ii).*ones(size(ydots))+0.8.*(rand(size(ydots))-0.5);
%         p = plot(xdots,ydots,'o','MarkerSize',2,'Color',[0.7 0.7 0.7])
        plot(xdots,ydots,'.','Color',[0.7 0.7 0.7])
        Fig7b{ii,2} = ydots;

        %         scatter(xdots,ydots,[],'k.','MarkerAlpha',0.2)
    end
end

PstatsexportationFig7b = Pstats(range2,3);
SD7b = Pstats(range2,4);
for i=1:size(xind2,2)
Fig7b{i,1} = xind2(i);
Fig7b{i,3} = PstatsexportationFig7b(i);
Fig7b{i,4} = mean(Fig7b{i,2}); % check if column 3 and 4 are equal
Fig7b{i,5} = SD7b(i);
end
Fig7bsorted = sortrows(Fig7b,1);
%%% code to add on significance bars
% groups={[Mind(1) Dind(1)] [Dind(1) Uind(1)] [Mind(2) Dind(2)] [Dind(2) Uind(2)]...
%     [Mind(3) Dind(3)] [Dind(3) Uind(3)] [Mind(4) Dind(4)] [Dind(4) Uind(4)]};
if significanceON == 1
    file_list = cell2mat(Pdata(:,4));
    p=[];
    for m = 1:4
        DevIndex = find(file_list==Dev(m));
        UnsupIndex = find(file_list==Unsup(m));
        MechIndex = find(file_list==Mech(m));

        if sigbars == 2
       %for only 2 sets of sig bars, Mech v Dev, and Unsup v Dev
        p(2*m-1)=ranksum(Pdata{MechIndex,2},Pdata{DevIndex,2});% %1 for Ppl, 2 for Pab, 3 for Pdi 
        p(2*m)=ranksum(Pdata{UnsupIndex,2},Pdata{DevIndex,2});% %1 for Ppl, 2 for Pab, 3 for Pdi
        groups{2*m-1,1}=[xind2(MechIndex) xind2(DevIndex)];
        groups{2*m,1}=[xind2(UnsupIndex) xind2(DevIndex)];
        elseif sigbars == 3
       %for 3 sets of sig bars, Mech v Dev, and Unsup v Dev
        p(3*m-2)=ranksum(Pdata{MechIndex,2},Pdata{UnsupIndex,2});%
        p(3*m-1)=ranksum(Pdata{MechIndex,2},Pdata{DevIndex,2});% 
        p(3*m)=ranksum(Pdata{UnsupIndex,2},Pdata{DevIndex,2});% 
        groups{3*m-2,1}=[xind2(MechIndex) xind2(UnsupIndex)];
        groups{3*m-1,1}=[xind2(MechIndex) xind2(DevIndex)];
        groups{3*m,1}=[xind2(UnsupIndex) xind2(DevIndex)];
     
        end
    end

save('Stats_Fig7b','groups','p');
sigstar(groups,p);    
end

%%%
errorbar(xind2,Pstats(range2,3),Pstats(range2,4),'.k')

% e8U = errorbar(xind2(Uind),Pstats(range2(Uind),3),Pstats(range2(Uind),4),'.');
% e8U.Color = ColorPl*0.7;
% e8M = errorbar(xind2(Mind),Pstats(range2(Mind),3),Pstats(range2(Mind),4),'.');
% e8M.Color = ColorAb*0.7;
% e8D = errorbar(xind2(Dind),Pstats(range2(Dind),3),Pstats(range2(Dind),4),'.');
% e8D.Color = ColorDi*0.7;

% plot(20.5.*ones(1,2),[0 2],'k',41.5.*ones(1,2),[0 2],'k',60.5.*ones(1,2),...
%     [0 2],'k',87.5.*ones(1,2),[0 2],'k',137.5.*ones(1,2),[0 2],'k') %lines between animals
b8.FaceColor = 'flat';
b8.CData = ones(length(b8.CData),1).*ColorAb;
b8.CData(Uind,:) = ones(length(intersect(Unsup,range2)),1).*ColorAb*0.6;
b8.CData(Mind,:) = ones(length(intersect(Mech,range2)),1).*ColorAb*1.8;
ylabel('\Delta Pressure_A_b [cmH_2O]')
% set(gca,'xtick',[])
xticks([0:4:20])
xticklabels({})
hold off

% 
% subplot(2,2,3)
% hold on
% b9 = bar(range3,Pstats(range3,3));
% errorbar(range3,Pstats(range3,3),Pstats(range3,4),'.')
% plot(20.5.*ones(1,2),[0 2],'k',41.5.*ones(1,2),[0 2],'k',60.5.*ones(1,2),...
%     [0 2],'k',87.5.*ones(1,2),[0 2],'k',137.5.*ones(1,2),[0 2],'k') %lines between animals
% b9.FaceColor = 'flat';
% b9.CData(ismember(range3,Unsup),:) = ones(length(intersect(Unsup,range3)),1).*ColorUn;
% b9.CData(ismember(range3,Mech),:) = ones(length(intersect(Mech,range3)),1).*ColorMech;
% ylabel('Change in Abdominal Pressure [cmH_2O]')
% hold off
% 
% subplot(2,2,4)
% hold on
% b10 = bar(range4,Pstats(range4,3));
% errorbar(range4,Pstats(range4,3),Pstats(range4,4),'.')
% plot(20.5.*ones(1,2),[0 2],'k',41.5.*ones(1,2),[0 2],'k',60.5.*ones(1,2),...
%     [0 2],'k',87.5.*ones(1,2),[0 2],'k',137.5.*ones(1,2),[0 2],'k') %lines between animals
% b10.FaceColor = 'flat';
% b10.CData(ismember(range4,Unsup),:) = ones(length(intersect(Unsup,range4)),1).*ColorUn;
% b10.CData(ismember(range4,Mech),:) = ones(length(intersect(Mech,range4)),1).*ColorMech;
% ylabel('Change in Abdominal Pressure [cmH_2O]')
% hold off




f5=figure('Position',[00 35 210 180]);
% subplot(2,2,1)
% hold on
% b11 = bar(range1,Pstats(range1,5));
% errorbar(range1,Pstats(range1,5),Pstats(range1,6),'.')
% plot(20.5.*ones(1,2),[0 2],'k',41.5.*ones(1,2),[0 2],'k',60.5.*ones(1,2),...
%     [0 2],'k',87.5.*ones(1,2),[0 2],'k',137.5.*ones(1,2),[0 2],'k') %lines between animals
% b11.FaceColor = 'flat';
% b11.CData(ismember(range1,Unsup),:) = ones(length(intersect(Unsup,range1)),1).*ColorUn;
% b11.CData(ismember(range1,Mech),:) = ones(length(intersect(Mech,range1)),1).*ColorMech;
% ylabel('Change in TransDiaphragm Pressure [cmH_2O]')
% hold off

% subplot(2,2,2)
hold on
b12 = bar(xind2,Pstats(range2,5));
if datadotsON ==1;
    for ii = 1:size(Pdata,1)
        ydots = Pdata{ii,3}; %1 for Ppl, 2 for Pab, 3 for Pdi
        xdots = xind2(ii).*ones(size(ydots))+0.8.*(rand(size(ydots))-0.5);
        plot(xdots,ydots,'.','Color',[0.7 0.7 0.7])
        Fig7c{ii,2} = ydots;
        %         scatter(xdots,ydots,[],'k.','MarkerAlpha',0.2)
    end
end

Pstatsexportation = Pstats(range2,5);
SD7c = Pstats(range2,6);
for i=1:size(xind2,2)
Fig7c{i,1} = xind2(i);
Fig7c{i,3} = Pstatsexportation(i);
Fig7c{i,4} = mean(Fig7c{i,2}); % check if column 3 and 4 are equal
Fig7c{i,5} = SD7c(i);
end
Fig7csorted = sortrows(Fig7c,1);



%%% code to add on significance bars
% groups={[Mind(1) Dind(1)] [Dind(1) Uind(1)] [Mind(2) Dind(2)] [Dind(2) Uind(2)]...
%     [Mind(3) Dind(3)] [Dind(3) Uind(3)] [Mind(4) Dind(4)] [Dind(4) Uind(4)]};
if significanceON == 1
    file_list = cell2mat(Pdata(:,4));
    p=[];
    for m = 1:4
        DevIndex = find(file_list==Dev(m));
        UnsupIndex = find(file_list==Unsup(m));
        MechIndex = find(file_list==Mech(m));
        if sigbars == 2
       %for only 2 sets of sig bars, Mech v Dev, and Unsup v Dev
        p(2*m-1)=ranksum(Pdata{MechIndex,3},Pdata{DevIndex,3});% %1 for Ppl, 2 for Pab, 3 for Pdi 
        p(2*m)=ranksum(Pdata{UnsupIndex,3},Pdata{DevIndex,3});% %1 for Ppl, 2 for Pab, 3 for Pdi
        groups{2*m-1,1}=[xind2(MechIndex) xind2(DevIndex)];
        groups{2*m,1}=[xind2(UnsupIndex) xind2(DevIndex)];
        elseif sigbars == 3
       %for 3 sets of sig bars, Mech v Dev, and Unsup v Dev
        p(3*m-2)=ranksum(Pdata{MechIndex,3},Pdata{UnsupIndex,3});%
        p(3*m-1)=ranksum(Pdata{MechIndex,3},Pdata{DevIndex,3});% 
        p(3*m)=ranksum(Pdata{UnsupIndex,3},Pdata{DevIndex,3});% 
        groups{3*m-2,1}=[xind2(MechIndex) xind2(UnsupIndex)];
        groups{3*m-1,1}=[xind2(MechIndex) xind2(DevIndex)];
        groups{3*m,1}=[xind2(UnsupIndex) xind2(DevIndex)];
     
        end
    end
end
save('Stats_Fig7c','groups','p');
sigstar(groups,p);
%%%
errorbar(xind2,Pstats(range2,5),Pstats(range2,6),'k.')
% e12U = errorbar(xind2(Uind),Pstats(range2(Uind),5),Pstats(range2(Uind),6),'.');
% e12U.Color = ColorPl*0.7;
% e12M = errorbar(xind2(Mind),Pstats(range2(Mind),5),Pstats(range2(Mind),6),'.');
% e12M.Color = ColorAb*0.7;
% e12D = errorbar(xind2(Dind),Pstats(range2(Dind),5),Pstats(range2(Dind),6),'.');
% e12D.Color = ColorDi*0.7;
% plot(20.5.*ones(1,2),[0 2],'k',41.5.*ones(1,2),[0 2],'k',60.5.*ones(1,2),...
%     [0 2],'k',87.5.*ones(1,2),[0 2],'k',137.5.*ones(1,2),[0 2],'k') %lines between animals
b12.FaceColor = 'flat';
b12.CData = ones(length(b12.CData),1).*ColorDi; %Device
b12.CData(Uind,:) = ones(length(intersect(Unsup,range2)),1).*ColorDi*0.6; %spont
b12.CData(Mind,:) = ones(length(intersect(Mech,range2)),1).*ColorDi*1.8; %mech
ylabel('\Delta Pressure_D_i [cmH_2O]')
% set(gca,'xtick',[])
xticks([0:4:20])
xticklabels({})
hold off

% 
% subplot(2,2,3)
% hold on
% b13 = bar(range3,Pstats(range3,5));
% errorbar(range3,Pstats(range3,5),Pstats(range3,6),'.')
% plot(20.5.*ones(1,2),[0 2],'k',41.5.*ones(1,2),[0 2],'k',60.5.*ones(1,2),...
%     [0 2],'k',87.5.*ones(1,2),[0 2],'k',137.5.*ones(1,2),[0 2],'k') %lines between animals
% b13.FaceColor = 'flat';
% b13.CData(ismember(range3,Unsup),:) = ones(length(intersect(Unsup,range3)),1).*ColorUn;
% b13.CData(ismember(range3,Mech),:) = ones(length(intersect(Mech,range3)),1).*ColorMech;
% ylabel('Change in TransDiaphragm Pressure [cmH_2O]')
% hold off
% 
% subplot(2,2,4)
% hold on
% b14 = bar(range4,Pstats(range4,5));
% errorbar(range4,Pstats(range4,5),Pstats(range4,6),'.')
% plot(20.5.*ones(1,2),[0 2],'k',41.5.*ones(1,2),[0 2],'k',60.5.*ones(1,2),...
%     [0 2],'k',87.5.*ones(1,2),[0 2],'k',137.5.*ones(1,2),[0 2],'k') %lines between animals
% b14.FaceColor = 'flat';
% b14.CData(ismember(range4,Unsup),:) = ones(length(intersect(Unsup,range4)),1).*ColorUn;
% b14.CData(ismember(range4,Mech),:) = ones(length(intersect(Mech,range4)),1).*ColorMech;
% ylabel('Change in TransDiaphragm Pressure [cmH_2O]')
% hold off
% 

%% Save data to Excel file (NBME request)
Letters = ['A'; 'B'; 'C'; 'D'; 'E'; 'F'; 'G'; 'H'; 'I'; 'J'; 'K'; 'L'];

for i=1:size(Letters,1)
writematrix(Fig7asorted{i,2}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', 'Fig7a', 'Range', strcat(Letters(i), num2str(2)));
writematrix(Fig7bsorted{i,2}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', 'Fig7b', 'Range', strcat(Letters(i), num2str(2)));
writematrix(Fig7csorted{i,2}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', 'Fig7c', 'Range', strcat(Letters(i), num2str(2)));
end

%% Save
pathWithFolderName =  strcat(pwd,'\Figures For Paper\');
Condition = 'ExcExp2BarChartTransDiPressure_3sigbars';
figCondition = strcat('',Condition);
% figPrefix = strcat('FN',num2str(FileNum));
if datadotsON ==1;
    prefix = 'with_data_pts'
else
    prefix =[]
end

figName = strcat(prefix,figCondition);
figFileName = strcat(pathWithFolderName,figName);

savefig(f3,strcat(figFileName,'PPl.fig')) %will save figure f as a .fig
exportgraphics(f3,strcat(figFileName,'PPl.eps'),'ContentType','vector') %will save figure f as a .png
savefig(f4,strcat(figFileName,'PAb.fig')) %will save figure f as a .fig
exportgraphics(f4,strcat(figFileName,'PAb.eps'),'ContentType','vector') %will save figure f as a .png
savefig(f5,strcat(figFileName,'PDi.fig')) %will save figure f as a .fig
exportgraphics(f5,strcat(figFileName,'PDi.eps'),'ContentType','vector') %will save figure f as a .png


