%generateFig_DeltaBeforeandAfter_1subplot

%modified from compareDeltaBeforeandAfter
clear all
close all


load('deltaComparisonVFdata_rev.mat')% from pullVolFlowDataforDeltaComparison
load('listImmediateChangeSegments.mat') %read this file to look at names of files and the associated notes

significanceON = 1;
datadotsON = 1;
% volflowdata = cell(74,6);
 % volflow data is {date Volpks Vollocs Flowpks Flowlocs 'Intact/Severed'}
 
 
 
% indep = [1 6 13 19 21 23 27]; %array of independent file numbers

vol = [];
flow = [];
n_Un = 0;
n_Sev = 0;
n_Ap = 0;
Unsup = [];
Severed = [];
nonSev = [];
Ap = [];
nonAp = [];
Sup = [];

% filelist = [1:8 25:28 33 34 36 37 51:54 69:72];
filelist = [5:8 25:28 33 34 36 37 51:54 69:72];%Exc exp 2
% filerange = [1:74];
N=1;
for i = filelist

    if i == 8
        Vind=find(volflowdata{i,3}>0*1000);%find the first index after 0 s
        Find=find(volflowdata{i,5}>0*1000);%find the first index after 0 s
    else
        Vind=find(volflowdata{i,3}>0*1000);%find the first index after 0 s
        Find=find(volflowdata{i,5}>0*1000);%find the first index after 0 s
    
    end
VFstats(i,1) = mean(volflowdata{i,2}(Vind(1):end)); %vol mean
VFstats(i,2) = std(volflowdata{i,2}(Vind(1):end));%vol std
VFstats(i,3) = mean(volflowdata{i,4}(Vind(1):end)); %flow mean
VFstats(i,4) = std(volflowdata{i,4}(Vind(1):end)); %flow std
VFstats(i,5) = str2num(list.name{i}(5));%Experiment number

VFdata{N,1} = volflowdata{i,2}(Vind(1):end); %vol
VFdata{N,2} = volflowdata{i,4}(Vind(1):end); %flow
VFdata{N,3} = i;%filenumber

%VFdata = {[vol] [flow] [filenum]}

%%%%%% calculate minute ventilation %%%%%%

totvent = sum(volflowdata{i,2}(Vind(1)+1:end)); %sum 1 fewer total vol to account for the time calculated for tottime
tottime = volflowdata{i,3}(end)-volflowdata{i,3}(Vind(1));
minvent = totvent*60/(tottime*0.001); %L/min

VFstats(i,6) = minvent;%min vent

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%vol/flow = [(data) (file number) (0=Unsup,1=Sup) (0=Severed,1=Intact) (0=not Apnea,1=Apnea)]
if contains(list.Var2(i),'Apn')
    Apnea = 1;
    n_Ap = n_Ap+1;
    Ap = [Ap i];
else
    Apnea = 0;
    nonAp = [nonAp i];
end
    if contains(list.Var2(i),'Un') || contains(list.Var2(i),'un') %true if unsupported
        
        if isempty(volflowdata{i,6}) %Intact
            vol = [vol ; volflowdata{i,2}(Vind(1):end)' i.*ones(size(volflowdata{i,2}(Vind(1):end)'))...
                zeros(size(volflowdata{i,2}(Vind(1):end)')) ones(size(volflowdata{i,2}(Vind(1):end)')) Apnea.*ones(size(volflowdata{i,2}(Vind(1):end)'))];
            flow = [flow ; volflowdata{i,4}(Vind(1):end)' i.*ones(size(volflowdata{i,4}(Vind(1):end)'))...
                zeros(size(volflowdata{i,4}(Vind(1):end)')) ones(size(volflowdata{i,4}(Vind(1):end)')) Apnea.*ones(size(volflowdata{i,4}(Vind(1):end)'))];
            nonSev = [nonSev i];
        elseif contains(volflowdata(i,6),'Severed') %Severed
            vol = [vol ; volflowdata{i,2}(Vind(1):end)' i.*ones(size(volflowdata{i,2}(Vind(1):end)'))...
                zeros(size(volflowdata{i,2}(Vind(1):end)')) zeros(size(volflowdata{i,2}(Vind(1):end)')) Apnea.*ones(size(volflowdata{i,2}(Vind(1):end)'))];
            flow = [flow ; volflowdata{i,4}(Vind(1):end)' i.*ones(size(volflowdata{i,4}(Vind(1):end)'))...
                zeros(size(volflowdata{i,4}(Vind(1):end)')) zeros(size(volflowdata{i,4}(Vind(1):end)')) Apnea.*ones(size(volflowdata{i,4}(Vind(1):end)'))];
            n_Sev = n_Sev+1;
            Severed = [Severed i];
        end
        n_Un=n_Un+1;  
        Unsup = [Unsup i];
    else %true if device supported
        if isempty(volflowdata{i,6}) %Intact
            vol = [vol ; volflowdata{i,2}(Vind(1):end)' i.*ones(size(volflowdata{i,2}(Vind(1):end)'))...
                ones(size(volflowdata{i,2}(Vind(1):end)')) ones(size(volflowdata{i,2}(Vind(1):end)')) Apnea.*ones(size(volflowdata{i,2}(Vind(1):end)'))];
            flow = [flow ; volflowdata{i,4}(Vind(1):end)' i.*ones(size(volflowdata{i,4}(Vind(1):end)'))...
                ones(size(volflowdata{i,4}(Vind(1):end)')) ones(size(volflowdata{i,4}(Vind(1):end)')) Apnea.*ones(size(volflowdata{i,4}(Vind(1):end)'))];
            nonSev = [nonSev i];
        elseif contains(volflowdata(i,6),'Severed') %Severed
            vol = [vol ; volflowdata{i,2}(Vind(1):end)' i.*ones(size(volflowdata{i,2}(Vind(1):end)'))...
                ones(size(volflowdata{i,2}(Vind(1):end)')) zeros(size(volflowdata{i,2}(Vind(1):end)')) Apnea.*ones(size(volflowdata{i,2}(Vind(1):end)'))];
            flow = [flow ; volflowdata{i,4}(Vind(1):end)' i.*ones(size(volflowdata{i,4}(Vind(1):end)'))...
                ones(size(volflowdata{i,4}(Vind(1):end)')) zeros(size(volflowdata{i,4}(Vind(1):end)')) Apnea.*ones(size(volflowdata{i,4}(Vind(1):end)'))];
            n_Sev = n_Sev+1;
            Severed = [Severed i];
        end
        Sup=[Sup i];
    end
  N=N+1;
end



%% Statistical Testing
%unpaired t test to compare before and after
 % volflow data is {date Volpks Vollocs Flowpks Flowlocs 'Intact/Severed'}
 vollist = cell(8,4);
 flowlist = cell(8,4);
indicesofdata = find(VFstats(:,1));
for iii = 1:length(indicesofdata) %create vectors
    n = indicesofdata(iii);
    if any(n==Ap) %if this value is an apnea trial 
        if any(n==Unsup)%this value is unsupported before device is turned on 
            vollist{VFstats(n,5),1} = volflowdata{n,2};
            flowlist{VFstats(n,5),1} = volflowdata{n,4};
        elseif any(n==Sup)%this value is supported after device is on
            vollist{VFstats(n,5),2} = volflowdata{n,2};
            flowlist{VFstats(n,5),2} = volflowdata{n,4};
        end
    elseif any(n==nonAp) %if this value is the end of a trial 
        if any(n==Unsup)%this value is unsupported before device is turned on 
            vollist{VFstats(n,5),3} = volflowdata{n,2};
            flowlist{VFstats(n,5),3} = volflowdata{n,4};
        else%this value is supported after device is on
            vollist{VFstats(n,5),4} = volflowdata{n,2};
            flowlist{VFstats(n,5),4} = volflowdata{n,4};
        end
    end
end
for E = [3,5,6,7,8]
[hOFFtoON_V(E),pOFFtoON_V(E),ciOFFtoON_V(E,:),stats] = ttest2(vollist{E,2}, vollist{E,1});
[hONtoOFF_V(E),pONtoOFF_V(E),ciONtoOFF_V(E,:),stats] = ttest2(vollist{E,4}, vollist{E,3});
[hOFFtoON_F(E),pOFFtoON_F(E),ciOFFtoON_F(E,:),stats] = ttest2(flowlist{E,2}, flowlist{E,1});
[hONtoOFF_F(E),pONtoOFF_F(E),ciONtoOFF_F(E,:),stats] = ttest2(flowlist{E,4}, flowlist{E,3});
% d = computeCohen_d(x1, x2, 'independent');
d(E) = computeCohen_d(vollist{E,2}, vollist{E,1}, 'independent');
end

% %paired T test
% %organize data into vectors of (before) and (after)
% vollist = zeros(8,4); %[dev off dev on, dev on dev off]
% flowlist = [];
% minVlist = [];
% indicesofdata = find(VFstats(:,1));
% for iii = 1:length(indicesofdata) %create vectors
%     n = indicesofdata(iii);
%     if any(n==range1) %if this value is an apnea trial 
% 
%         if any(n==Unsup)%this value is unsupported before device is turned on 
%             vollist(VFstats(n,5),1)= VFstats(n,1);%dev off, before
%             flowlist(VFstats(n,5),1)= VFstats(n,3);
%             minVlist(VFstats(n,5),1)= VFstats(n,6);
%         elseif any(n==Sup)%this value is supported after device is on
%             vollist(VFstats(n,5),2) =VFstats(n,1);%dev on, after
%             flowlist(VFstats(n,5),2)= VFstats(n,3);
%             minVlist(VFstats(n,5),2)= VFstats(n,6);
% 
%         end
%     elseif any(n==range2) %if this value is the end of a trial
% 
%             if any(n==Unsup)%this value is unsupported before device is turned on 
%             vollist(VFstats(n,5),3) =VFstats(n,1);%dev on, before
%             flowlist(VFstats(n,5),3)= VFstats(n,3);
%             minVlist(VFstats(n,5),3)= VFstats(n,6);
%         else%this value is supported after device is on
%             vollist(VFstats(n,5),4) =VFstats(n,1);%dev off, after
%             flowlist(VFstats(n,5),4)= VFstats(n,3);
%             minVlist(VFstats(n,5),4)= VFstats(n,6);
%         end
%     end
% end
% 
% % deltaOFFtoON = [vollist(:,2)-vollist(:,1) flowlist(:,2)-flowlist(:,1) minVlist(:,2)-minVlist(:,1)];
% % deltaONtoOFF = [vollist(:,4)-vollist(:,3) flowlist(:,4)-flowlist(:,3) minVlist(:,4)-minVlist(:,3)];
% 
% % figure
% % histogram(deltaOFFtoON(find(deltaOFFtoON(:,1))),6)
% % figure
% % histogram(deltaOFFtoON(find(deltaOFFtoON(:,2))),6)
% % figure
% % histogram(deltaOFFtoON(find(deltaOFFtoON(:,3))),6)
% 
% set = find(vollist(:,2));
% % [p,h,stats]=signrank(flowlist(set,4), flowlist(set,3))%before after volume
% % [h,p] = ttest(minVlist(set,3), minVlist(set,4));









% %% Preliminary Plotting
% 
% 
% vol(:,1)=vol(:,1).*1000;
% 
% UnsupIndV = find(vol(:,3)==0);
% UnsupIndF = find(flow(:,3)==0);
% SupIndV = find(vol(:,3)==1);
% SupIndF = find(flow(:,3)==1);
% 
% SeverIndV = find(vol(:,4)==0);
% SeverIndF = find(flow(:,4)==0);
% IntactIndV = find(vol(:,4)==1);
% IntactIndF = find(flow(:,4)==1);
% 
% 
% ColorGuide = ones(n_Un,1).*[0.5 0 0.5];
% 
% 
% figure
% hold on
% b1 = bar(filerange,VFstats(filerange,1).*1000);
% errorbar(filerange,VFstats(filerange,1).*1000,VFstats(filerange,2).*1000,'.')
% plot(4.5.*ones(1,2),[0 200],'k',9.5.*ones(1,2),[0 200],'k',28.5.*ones(1,2),...
%     [0 200],'k',38.5.*ones(1,2),[0 200],'k',58.5.*ones(1,2),[0 200],'k') %lines between animals
% plot(Ap,250.*ones(size(Ap)),'*')%indicate apnea points
% hold off
% b1.FaceColor = 'flat';
% b1.CData(Unsup,:) = ColorGuide;
% b1.CData(Severed,:) = b1.CData(Severed,:).*1.6;
% xlabel('File Number')
% ylabel('Tidal Volume [mL]')
% 
% 
% figure
% hold on
% b2 = bar(filerange,VFstats(filerange,3));
% errorbar(filerange,VFstats(filerange,3),VFstats(filerange,4),'.')
% plot(4.5.*ones(1,2),[0 0.7],'k',9.5.*ones(1,2),[0 0.7],'k',28.5.*ones(1,2),[0 0.7],'k',38.5.*ones(1,2),[0 0.7],'k',58.5.*ones(1,2),[0 0.7],'k')
% plot(Ap,0.75.*ones(size(Ap)),'*')%indicate apnea points
% hold off
% b2.FaceColor = 'flat';
% b2.CData(Unsup,:) = ones(length(Unsup),1)*[.5 0 .5];
% b2.CData(Severed,:) = b2.CData(Severed,:).*1.6;
% xlabel('File Number')
% ylabel('Peak Inspiratory Flow [L/s]')



%% Plot into 4 figures
range1 = intersect(Ap,nonSev);
range2 = intersect(nonAp,nonSev);
% range3 = intersect(Ap,Severed);
% range4 = intersect(nonAp,Severed);

% filerange = [1:74];
filerange = VFstats(:,5);

%%%%% establish the order they appear on the x axis


% %create the order matrix
% deltas = VFstats(range1(2:2:end),1)-VFstats(range1(1:2:end-1),1);%sort by steady state delta
% % deltas = VFstats(range2(2:2:end),1)-VFstats(range2(1:2:end-1),1);%sort by apnea delta
% [delsort,delI] = sort(deltas,'descend');

%create the order matrix based on minute ventilation delta
deltas = VFstats(range1(2:2:end),6)-VFstats(range1(1:2:end-1),6);%apnea
% deltas = VFstats(range2(2:2:end),6)-VFstats(range2(1:2:end-1),6);%steady
% state
[delsort,delI] = sort(deltas,'descend');

% Exp = [2 3 5:8];
Exp = [3 5:8];
[ord,I] = sort(Exp(delI));
% order = [0 I(1:2) 0 I(3:end)];
%   %Exp     1 2 3 4 5 6 7 8
% % order = [0 X X 0 X X X X];%order (will dictate which experiment is 1st, 2nd, 3rd,... on the graph)
order = [0 0 I(1) 0 I(2:end)];
  %Exp     1 2 3 4 5 6 7 8
% order = [0 0 X 0 X X X X];%order (will dictate which experiment is 1st, 2nd, 3rd,... on the graph)
const1 = ones(1,length(range1));
const1(1:2:end-1) = 2; 
const1 = const1 +4;
xind1 = 7.*order(filerange(range1))-const1;
const2 = ones(1,length(range2));
const2(1:2:end-1) = 2; 
const2 = const2 + 1;
xind2 = 7.*order(filerange(range2))-const2;


xindtot_ = [reshape(xind1,[2 5])' reshape(xind2,[2 5])'];
xindtot = reshape(xindtot_',[1 20]);
%%
f3=figure('Position',[00 50 350 160]);
% subplot(1,2,1)
hold on
b3 = bar(xind1,VFstats(range1,1).*1000);
if datadotsON ==1;
file_list = cell2mat(VFdata(:,3));
    for ii = 1:length(range1)
%VFdata = {[vol] [flow] [filenum]}
        k = find(file_list==range1(ii));
        ydots = VFdata{k,1}*1000; %1 for vol, 2 for flow
        xdots = xind1(ii).*ones(size(ydots))+0.8.*(rand(size(ydots))-0.5);
%         p = plot(xdots,ydots,'o','MarkerSize',2,'Color',[0.7 0.7 0.7])
        plot(xdots,ydots,'.','Color',[0.7 0.7 0.7])
        %         scatter(xdots,ydots,[],'k.','MarkerAlpha',0.2)
        Fig3dbefore{ii,1} = VFdata{k,1}*1000;
        Fig3dbefore{ii,2} = VFdata{k,3};
    end

end

for i=1:size(xind1,2)
Fig3dbefore{i,3} = xind1(i);
end

%%% code to add on significance bars
if significanceON == 1
    file_list = cell2mat(VFdata(:,3));
    p1 = [];
%     p2 = [];
    for m = 1:5
        UnsupIndex1 = find(file_list==range1(2*m-1));
        SupIndex1 = find(file_list==range1(2*m));
        p1(m)=ranksum(VFdata{UnsupIndex1,1},VFdata{SupIndex1,1});% %1 for vol, 2 for flow,        groups{2*m-1,1}=[xind2(MechIndex) xind2(DevIndex)];
        groups1{m} = [xind1(2*m-1) xind1(2*m)];
%         SupIndex2 = find(file_list==range2(2*m-1));
%         UnsupIndex2 = find(file_list==range2(2*m));
%         p2(m)=ranksum(VFdata{UnsupIndex2,1},VFdata{SupIndex2,1});% %1 for vol, 2 for flow,        groups{2*m-1,1}=[xind2(MechIndex) xind2(DevIndex)];
%         groups2{m} = [xind2(2*m-1) xind2(2*m)];
    end
sigstar(groups1,p1);
% sigstar(groups2,p2);
end
errorbar(xind1,VFstats(range1,1).*1000,VFstats(range1,2).*1000,'.k')
% plot(4.5.*ones(1,2),[0 200],'k',9.5.*ones(1,2),[0 200],'k',28.5.*ones(1,2),...
%     [0 200],'k',38.5.*ones(1,2),[0 200],'k',58.5.*ones(1,2),[0 200],'k') %lines between animals
b3.FaceColor = 'flat';
b3.CData = ones(length(b3.CData),1).*(hex2rgb('#E0AF1A'));
b3.CData(ismember(range1,Unsup),:) = ones(length(intersect(Unsup,range1)),1).*(hex2rgb('#FFE281'));
ylabel('Tidal Volume [mL]')
% set(gca,'xtick',[])
xlabel(join(string(Exp(delI))))
xticks([-0.5:3.5:34.5])
xticklabels({})
hold off

% subplot(1,2,2)
hold on
b4 = bar(xind2,VFstats(range2,1).*1000);
if datadotsON ==1;
file_list = cell2mat(VFdata(:,3));
    for ii = 1:length(range2)
%VFdata = {[vol] [flow] [filenum]}
        k = find(file_list==range2(ii));
        ydots = VFdata{k,1}*1000; %1 for vol, 2 for flow
        xdots = xind2(ii).*ones(size(ydots))+0.8.*(rand(size(ydots))-0.5);
%         p = plot(xdots,ydots,'o','MarkerSize',2,'Color',[0.7 0.7 0.7])
        plot(xdots,ydots,'.','Color',[0.7 0.7 0.7])
        %         scatter(xdots,ydots,[],'k.','MarkerAlpha',0.2)
        Fig3dafter{ii,1} = VFdata{k,1}*1000;
        Fig3dafter{ii,2} = VFdata{k,3};
    end
end

for i=1:size(xind2,2)
Fig3dafter{i,3} = xind2(i);
end

Fig3d = [Fig3dbefore ; Fig3dafter];
Fig3dsorted = sortrows(Fig3d,3);
for i=1:size(Fig3dsorted,1)
Fig3dsorted{i,4} = mean(Fig3dsorted{i,1});
Fig3dsorted{i,5} = std(Fig3dsorted{i,1});

end

if significanceON == 1
    file_list = cell2mat(VFdata(:,3));
%     p1 = [];
    p2 = [];
    for m = 1:5
%         UnsupIndex1 = find(file_list==range1(2*m-1));
%         SupIndex1 = find(file_list==range1(2*m));
%         p1(m)=ranksum(VFdata{UnsupIndex1,1},VFdata{SupIndex1,1});% %1 for vol, 2 for flow,        groups{2*m-1,1}=[xind2(MechIndex) xind2(DevIndex)];
%         groups1{m} = [xind1(2*m-1) xind1(2*m)];
        SupIndex2 = find(file_list==range2(2*m-1));
        UnsupIndex2 = find(file_list==range2(2*m));
        p2(m)=ranksum(VFdata{UnsupIndex2,1},VFdata{SupIndex2,1});% %1 for vol, 2 for flow,        groups{2*m-1,1}=[xind2(MechIndex) xind2(DevIndex)];
        groups2{m} = [xind2(2*m-1) xind2(2*m)];
    end
% sigstar(groups1,p1);
sigstar(groups2,p2);
end
errorbar(xind2,VFstats(range2,1).*1000,VFstats(range2,2).*1000,'.k')
% plot(4.5.*ones(1,2),[0 200],'k',9.5.*ones(1,2),[0 200],'k',28.5.*ones(1,2),...
%     [0 200],'k',38.5.*ones(1,2),[0 200],'k',58.5.*ones(1,2),[0 200],'k') %lines between animals
b4.FaceColor = 'flat';
b4.CData = ones(length(b4.CData),1).*(hex2rgb('#E0AF1A'));
b4.CData(ismember(range2,Unsup),:) = ones(length(intersect(Unsup,range2)),1).*(hex2rgb('#FFE281'));
xticks([-0.5:3.5:34.5])
xlim([0 34])
hold off

save('Stats_Fig3d', 'groups1','p1','groups2','p2');

Fig3dMeanSD = [xind1' VFstats(range1,1).*1000 VFstats(range1,2).*1000;
    xind2' VFstats(range2,1).*1000 VFstats(range2,2).*1000];
Fig3dMeanSDsorted = sortrows(Fig3dMeanSD,1);

f4=figure('Position',[00 335 350 160]);
% subplot(1,2,1)
hold on
b7 = bar(xind1,VFstats(range1,3));
if datadotsON ==1;
file_list = cell2mat(VFdata(:,3));
    for ii = 1:length(range1)
%VFdata = {[vol] [flow] [filenum]}
        k = find(file_list==range1(ii));
        ydots = VFdata{k,2}; %1 for vol, 2 for flow
        xdots = xind1(ii).*ones(size(ydots))+0.8.*(rand(size(ydots))-0.5);
%         p = plot(xdots,ydots,'o','MarkerSize',2,'Color',[0.7 0.7 0.7])
        plot(xdots,ydots,'.','Color',[0.7 0.7 0.7])
        %         scatter(xdots,ydots,[],'k.','MarkerAlpha',0.2)
        Fig3cbefore{ii,1} = VFdata{k,2};
        Fig3cbefore{ii,2} = VFdata{k,3};
    end
end

for i=1:size(xind1,2)
Fig3cbefore{i,3} = xind1(i);
end

if significanceON == 1
    file_list = cell2mat(VFdata(:,3));
    p1 = [];
%     p2 = [];
    for m = 1:5
        UnsupIndex1 = find(file_list==range1(2*m-1));
        SupIndex1 = find(file_list==range1(2*m));
        p1(m)=ranksum(VFdata{UnsupIndex1,1},VFdata{SupIndex1,1});% %1 for vol, 2 for flow,        groups{2*m-1,1}=[xind2(MechIndex) xind2(DevIndex)];
        groups1{m} = [xind1(2*m-1) xind1(2*m)];
%         SupIndex2 = find(file_list==range2(2*m-1));
%         UnsupIndex2 = find(file_list==range2(2*m));
%         p2(m)=ranksum(VFdata{UnsupIndex2,1},VFdata{SupIndex2,1});% %1 for vol, 2 for flow,        groups{2*m-1,1}=[xind2(MechIndex) xind2(DevIndex)];
%         groups2{m} = [xind2(2*m-1) xind2(2*m)];
    end
sigstar(groups1,p1);
% sigstar(groups2,p2);
end
errorbar(xind1,VFstats(range1,3),VFstats(range1,4),'.k')
% plot(4.5.*ones(1,2),[0 0.7],'k',9.5.*ones(1,2),[0 0.7],'k',28.5.*ones(1,2),[0 0.7],'k',38.5.*ones(1,2),[0 0.7],'k',58.5.*ones(1,2),[0 0.7],'k')
b7.FaceColor = 'flat';
% b7.CData(ismember(range1,Unsup),:) = ones(length(intersect(Unsup,range1)),1).*[0.5 0 0.5];
b7.CData = ones(length(b7.CData),1).*(hex2rgb('#32478A'));
b7.CData(ismember(range1,Unsup),:) = ones(length(intersect(Unsup,range1)),1).*(hex2rgb('#6493DB'));
hold off

% subplot(1,2,2)
hold on
b8 = bar(xind2,VFstats(range2,3));
if datadotsON ==1;
file_list = cell2mat(VFdata(:,3));
    for ii = 1:length(range2)
%VFdata = {[vol] [flow] [filenum]}
        k = find(file_list==range2(ii));
        ydots = VFdata{k,2}; %1 for vol, 2 for flow
        xdots = xind2(ii).*ones(size(ydots))+0.8.*(rand(size(ydots))-0.5);
%         p = plot(xdots,ydots,'o','MarkerSize',2,'Color',[0.7 0.7 0.7])
        plot(xdots,ydots,'.','Color',[0.7 0.7 0.7])
        %         scatter(xdots,ydots,[],'k.','MarkerAlpha',0.2)
        Fig3cafter{ii,1} = VFdata{k,2};
        Fig3cafter{ii,2} = VFdata{k,3};
    end
end

for i=1:size(xind2,2)
Fig3cafter{i,3} = xind2(i);
end

Fig3c = [Fig3cbefore ; Fig3cafter];
Fig3csorted = sortrows(Fig3c,3);
for i=1:size(Fig3csorted,1)
Fig3csorted{i,4} = mean(Fig3csorted{i,1});
Fig3csorted{i,5} = std(Fig3csorted{i,1});

end

if significanceON == 1
    file_list = cell2mat(VFdata(:,3));
%     p1 = [];
    p2 = [];
    for m = 1:5
%         UnsupIndex1 = find(file_list==range1(2*m-1));
%         SupIndex1 = find(file_list==range1(2*m));
%         p1(m)=ranksum(VFdata{UnsupIndex1,1},VFdata{SupIndex1,1});% %1 for vol, 2 for flow,        groups{2*m-1,1}=[xind2(MechIndex) xind2(DevIndex)];
%         groups1{m} = [xind1(2*m-1) xind1(2*m)];
        SupIndex2 = find(file_list==range2(2*m-1));
        UnsupIndex2 = find(file_list==range2(2*m));
        p2(m)=ranksum(VFdata{UnsupIndex2,1},VFdata{SupIndex2,1});% %1 for vol, 2 for flow,        groups{2*m-1,1}=[xind2(MechIndex) xind2(DevIndex)];
        groups2{m} = [xind2(2*m-1) xind2(2*m)];
    end
% sigstar(groups1,p1);
sigstar(groups2,p2);
end
errorbar(xind2,VFstats(range2,3),VFstats(range2,4),'.k')
% plot(4.5.*ones(1,2),[0 0.7],'k',9.5.*ones(1,2),[0 0.7],'k',28.5.*ones(1,2),[0 0.7],'k',38.5.*ones(1,2),[0 0.7],'k',58.5.*ones(1,2),[0 0.7],'k')
b8.FaceColor = 'flat';
% b8.CData(ismember(range2,Unsup),:) = ones(length(intersect(Unsup,range2)),1).*[0.5 0 0.5];
b8.CData = ones(length(b8.CData),1).*(hex2rgb('#32478A'));
b8.CData(ismember(range2,Unsup),:) = ones(length(intersect(Unsup,range2)),1).*(hex2rgb('#6493DB'));
ylabel('Peak Inspiratory Flow [L/s]')
% set(gca,'xtick',[])
xticks([-0.5:3.5:34.5])
xticklabels({})
xlabel(join(string(Exp(delI))))
ylim([0 0.65])
xlim([0 34])
hold off

save('Stats_Fig3c', 'groups1','p1','groups2','p2');
Fig3cMeanSD = [xind1' VFstats(range1,3) VFstats(range1,4);
    xind2' VFstats(range2,3) VFstats(range2,4)];
Fig3cMeanSDsorted = sortrows(Fig3cMeanSD,1);

f5=figure('Position',[00 580 350 160]);
% subplot(1,2,1)
hold on
b9 = bar(xind1,VFstats(range1,6));
% errorbar(xind1,VFstats(range1,6),VFstats(range1,4),'.k')
% plot(4.5.*ones(1,2),[0 0.7],'k',9.5.*ones(1,2),[0 0.7],'k',28.5.*ones(1,2),[0 0.7],'k',38.5.*ones(1,2),[0 0.7],'k',58.5.*ones(1,2),[0 0.7],'k')
b9.FaceColor = 'flat';
% b7.CData(ismember(range1,Unsup),:) = ones(length(intersect(Unsup,range1)),1).*[0.5 0 0.5];
b9.CData = ones(length(b9.CData),1).*(hex2rgb('#9b3827'));
b9.CData(ismember(range1,Unsup),:) = ones(length(intersect(Unsup,range1)),1).*(hex2rgb('#cf7c49'));
ylabel('Minute Ventilation [L/min]')
ylim([0 6])
% set(gca,'xtick',[])
xticks([-0.5:3.5:34.5])
xticklabels({})
xlabel(join(string(Exp(delI))))
% ylim([0 0.65])
hold off

% subplot(1,2,2)
hold on
b10 = bar(xind2,VFstats(range2,6));
% errorbar(xind2,VFstats(range2,6),VFstats(range2,4),'.k')
% plot(4.5.*ones(1,2),[0 0.7],'k',9.5.*ones(1,2),[0 0.7],'k',28.5.*ones(1,2),[0 0.7],'k',38.5.*ones(1,2),[0 0.7],'k',58.5.*ones(1,2),[0 0.7],'k')
b10.FaceColor = 'flat';
% b8.CData(ismember(range2,Unsup),:) = ones(length(intersect(Unsup,range2)),1).*[0.5 0 0.5];
b10.CData = ones(length(b10.CData),1).*(hex2rgb('#9b3827'));
b10.CData(ismember(range2,Unsup),:) = ones(length(intersect(Unsup,range2)),1).*(hex2rgb('#cf7c49'));
xlim([0 34])
hold off

Fig3e = [xind1' VFstats(range1,6);
    xind2' VFstats(range2,6)];
Fig3esorted = sortrows(Fig3e,1);

%% Save
pathWithFolderName =  strcat(pwd,'\Figures For Paper\');
Condition = 'Together_ExcExp2_DeltaBarChart_Sorted_Ap_revisions';
figCondition = strcat('',Condition);
% figPrefix = strcat('FN',num2str(FileNum));
if significanceON ==1;
    prefix = 'with_significance';
else
    prefix =[];
end
if datadotsON ==1;
    prefix = strcat('with_datadots',prefix);
else
    prefix =prefix;
end
figName = strcat(prefix,figCondition);
figFileName = strcat(pathWithFolderName,figName);

savefig(f3,strcat(figFileName,'Vol.fig')) %will save figure f as a .fig
exportgraphics(f3,strcat(figFileName,'Vol.eps'),'ContentType','vector') %will save figure f as a .eps
savefig(f4,strcat(figFileName,'Flow.fig')) %will save figure f as a .fig
exportgraphics(f4,strcat(figFileName,'Flow.eps'),'ContentType','vector') %will save figure f as a .eps
savefig(f5,strcat(figFileName,'MinVent.fig')) %will save figure f as a .fig
exportgraphics(f5,strcat(figFileName,'MinVent.eps'),'ContentType','vector') %will save figure f as a .eps


%% Save data to Excel file (NBME request)

Letters = ['A'; 'B'; 'C'; 'D'; 'E'; 'F'; 'G'; 'H'; 'I'; 'J'; 'K'; 'L'; 'M'; 'N'; 'O'; 'P'; 'Q'; 'R'; 'S'; 'T'];

for i=1:size(Letters,1)
writematrix(Fig3csorted{i,1}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', 'Fig3c', 'Range', strcat(Letters(i), num2str(2)));
end

for i=1:size(Letters,1)
writematrix(Fig3dsorted{i,1}', 'SourceData_nBME-21-2902.xlsx', 'Sheet', 'Fig3d', 'Range', strcat(Letters(i), num2str(2)));
end

for i=1:size(Letters,1)
writematrix(Fig3esorted(i,2), 'SourceData_nBME-21-2902.xlsx', 'Sheet', 'Fig3e', 'Range', strcat(Letters(i), num2str(2)));
end

