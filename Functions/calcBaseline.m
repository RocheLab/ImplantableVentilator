function [baseline,SDbaseline] = calcBaseline(date)
%calcBaseline Generates an average and standard deviation of baseline ABGs
%for a given experimental date
%   date = the date of the experiment
%   baseline = vector of baseline values [pH pCO2 pO2 sO2]

%load ABG data that was imported from excel file
ABGfilename = [num2str(date) 'ABGData.mat'];
load(ABGfilename)


ind = find(strcmp(ABGDataCompiled{:,7},"Baseline"));
baseline = [mean(ABGDataCompiled{ind,3}) mean(ABGDataCompiled{ind,4}) mean(ABGDataCompiled{ind,5}) mean(ABGDataCompiled{ind,6})];
SDbaseline = [std(ABGDataCompiled{ind,3}) std(ABGDataCompiled{ind,4}) std(ABGDataCompiled{ind,5}) std(ABGDataCompiled{ind,6})];

end

