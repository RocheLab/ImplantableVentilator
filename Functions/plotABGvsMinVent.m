function [FitTable,relFitTable] = plotABGvsMinVent(date,T,ExperimentType)
%plotABGvsMinVent
%will pull data about avg min vent from files named 'dateABGandMinVentData.mat'
%will take the nonzero values and plot them as a function of pH
%   date = vector of dates of the experiments [date date]
%   T = the period of time (in seconds) over which to calculate an average minute
%   ventilation
%   ExperimentType = 'Apnea' or 'Maintenence' or 'Baseline'
%   FitTable = table of linear regression fit lines [slope y-int] for the
%   ABG data x minute ventilation
%   relFitTable = table of linear regression fit lines [slope y-int] for the
%   ABG data normalized to a delta from avg baseline per animal x minute ventilation

%cell of color labels 
% C = {'k','b','r','g','y',[.5 .6 .7],[.8 .2 .6]};
C = {'b','r','g','y',[.5 .6 .7],[.8 .2 .6]};

figure
hold all 

minventTname = sprintf('Avg Min Vent with Period T = %d[s]',T);
title(minventTname)
for d = 1:length(date)
ABGandMinVentfilename = [num2str(date(d)) 'ABGandMinVentData.mat'];
load(ABGandMinVentfilename)

ABGTypeDatafilename = [num2str(date(d)) 'ABGData.mat'];
load(ABGTypeDatafilename)

%find correct column

col = find(strcmp(ABGMinVentTable.Properties.VariableNames,minventTname));
if any(col)==0
    continue
end

%% Regression
%absolute values%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fullind = find(ABGMinVentTable{:,col}) %which ones have a avg min vent calculated

m=1;
ind = [];
for ii = 1:length(fullind)
    switch ExperimentType
    case 'Apnea'
        if modifiedABGDataCompiled{fullind(ii),7}=="Apnea"
            ind(m) = fullind(ii);
            m=m+1;
        end
    case 'Maintenence'
        if modifiedABGDataCompiled{fullind(ii),7}=="Maintenence"
            ind(m) = fullind(ii);
            m=m+1;
        end
    case 'Baseline'
        if modifiedABGDataCompiled{fullind(ii),7}=="Baseline"
            ind(m) = fullind(ii);
            m=m+1;
        end
    end    
        
    end
    

end



MinVentVector = ABGMinVentTable{ind,col};
pHvalues = ABGMinVentTable{ind,2};
pCO2values = ABGMinVentTable{ind,3};
pO2values = ABGMinVentTable{ind,4};
sO2values = ABGMinVentTable{ind,5};
    
[fpH] = polyfit(MinVentVector,pHvalues,1);
[fpCO2] = polyfit(MinVentVector,pCO2values,1);
[fpO2] = polyfit(MinVentVector,pO2values,1);
[fsO2] = polyfit(MinVentVector,sO2values,1);


YfpH = polyval(fpH,MinVentVector);
YfpCO2 = polyval(fpCO2,MinVentVector);
YfpO2 = polyval(fpO2,MinVentVector);
YfsO2 = polyval(fsO2,MinVentVector);

pHRsq = Rsquared(pHvalues,YfpH);
pCO2Rsq = Rsquared(pCO2values,YfpCO2);
pO2Rsq = Rsquared(pO2values,YfpO2);
sO2Rsq = Rsquared(sO2values,YfsO2);

if d==1
    FitTable = table({num2str(date(d))},fpH,pHRsq,fpCO2,pCO2Rsq,fpO2,pO2Rsq,fsO2,sO2Rsq);
else
    newrow = table({num2str(date(d))},fpH,pHRsq,fpCO2,pCO2Rsq,fpO2,pO2Rsq,fsO2,sO2Rsq);
    FitTable = [FitTable ; newrow];
end

%relative values%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
relpHvalues = ABGAdjBaseline{ind,7};
relpCO2values = ABGAdjBaseline{ind,8};
relpO2values = ABGAdjBaseline{ind,9};
relsO2values = ABGAdjBaseline{ind,10};
    
[relfpH] = polyfit(MinVentVector,relpHvalues,1);
[relfpCO2] = polyfit(MinVentVector,relpCO2values,1);
[relfpO2] = polyfit(MinVentVector,relpO2values,1);
[relfsO2] = polyfit(MinVentVector,relsO2values,1);


relYfpH = polyval(relfpH,MinVentVector);
relYfpCO2 = polyval(relfpCO2,MinVentVector);
relYfpO2 = polyval(relfpO2,MinVentVector);
relYfsO2 = polyval(relfsO2,MinVentVector);

relpHRsq = Rsquared(relpHvalues,relYfpH);
relpCO2Rsq = Rsquared(relpCO2values,relYfpCO2);
relpO2Rsq = Rsquared(relpO2values,relYfpO2);
relsO2Rsq = Rsquared(relsO2values,relYfsO2);

if d==1
    relFitTable = table({num2str(date(d))},relfpH,relpHRsq,relfpCO2,relpCO2Rsq,relfpO2,relpO2Rsq,relfsO2,relsO2Rsq);
else
    relnewrow = table({num2str(date(d))},relfpH,relpHRsq,relfpCO2,relpCO2Rsq,relfpO2,relpO2Rsq,relfsO2,relsO2Rsq);
    relFitTable = [relFitTable ; relnewrow];
end



%% Plot
    

    figure(1)
    subplot (2,2,1)
%     plot(MinVentVector,pHvalues,'marker','x','color',C{d},'LineStyle', 'none',V,YfpH,'color')
    plot(MinVentVector,YfpH,'-',MinVentVector,pHvalues,'x','color',C{d})
    xlabel('Average Minute Ventilation [L]')
    ylabel('pH')
    hold on
    
    subplot (2,2,2)
    plot(MinVentVector,YfpCO2,'-',MinVentVector,pCO2values,'x','color',C{d})
    xlabel('Average Minute Ventilation [L]')
    ylabel('pCO2 [mmHg]')
    hold on

    subplot (2,2,3)
    plot(MinVentVector,YfpO2,'-',MinVentVector,pO2values,'x','color',C{d})
    xlabel('Average Minute Ventilation [L]')
    ylabel('pO2 [mmHg]')
    hold on

    subplot (2,2,4)
    plot(MinVentVector,YfsO2,'-',MinVentVector,sO2values,'x','color',C{d})
    xlabel('Average Minute Ventilation [L]')
    ylabel('sO2 [%]')
    hold on
    
    
    figure(2)
    subplot (2,2,1)
%     plot(MinVentVector,pHvalues,'marker','x','color',C{d},'LineStyle', 'none',V,YfpH,'color')
    plot(MinVentVector,relYfpH,'-',MinVentVector,relpHvalues,'x','color',C{d})
    xlabel('Average Minute Ventilation [L]')
    ylabel('pH')
    hold on
    
    subplot (2,2,2)
    plot(MinVentVector,relYfpCO2,'-',MinVentVector,relpCO2values,'x','color',C{d})
    xlabel('Average Minute Ventilation [L]')
    ylabel('pCO2 [mmHg]')
    hold on

    subplot (2,2,3)
    plot(MinVentVector,relYfpO2,'-',MinVentVector,relpO2values,'x','color',C{d})
    xlabel('Average Minute Ventilation [L]')
    ylabel('pO2 [mmHg]')
    hold on

    subplot (2,2,4)
    plot(MinVentVector,relYfsO2,'-',MinVentVector,relsO2values,'x','color',C{d})
    xlabel('Average Minute Ventilation [L]')
    ylabel('sO2 [%]')
    hold on
end


