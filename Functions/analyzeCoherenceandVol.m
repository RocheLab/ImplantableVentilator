function analyzeCoherenceandVol(wcoh,f,VolAutoData,time_interp,f_interp,coherences_D,vpks,vlocs)
%analyzeCoherenceandVol: 
%   these variables are stored in the file name 
%   'coherenceVars_FN_',num2str(FileNum),'_',Condition,'.mat' from the
%   runWaveletCoherenceTest.m code
%   
%   wcoh = array of coherence values of flow and PAct over time
%   f = frequency values corresponding to wcoh
%   VolAutoData = the volume data corrected with the spirometric data
%   time_interp = the time vector corresponding to f_interp that spans the
%   same time as wcoh
%   f_interp = the interpolated frequency ridge (a dynamic value of RR in
%   Hz over time
%   fStat = the static frequency that is the average RR over a block
%   coherences_D = coherences from the dynamic freq coherence values
%   coherences_S = coherences from the static freq coherence values
%   vpks = peaks of VolAutoData
%   vlocs = peak locations of VolAutoData


figure('Position',[500 430 500 300]);
scatter(coherences_D(vlocs),vpks,[],vlocs)
ylabel('Tidal Volume [L]')
xlabel('Coherence')
title('Dynamic Frequency Coherence vs Tidal Volume')
colormap(gca,'autumn')

% figure('Position',[500 50 500 300]);
% scatter(coherences_S(vlocs),vpks,[],vlocs)
% ylabel('Tidal Volume [L]')
% xlabel('Coherence')
% title('Static Freqeuncy Coherence vs Tidal Volume')
% colormap(gca,'autumn')
% 
% % bucket coherences and evaluate their spread
% figure('Position',[500 50 500 300]);
groupings_D = discretize(coherences_D(vlocs),0:0.1:1);
% boxplot(vpks,groupings_D);
% ylabel('Tidal Volume [L]')
% xlabel('Coherence')
% title('Dynamic Frequency Coherence vs Tidal Volume')

% figure('Position',[1000 50 500 300]);
% groupings_S = discretize(coherences_S(vlocs),0:0.1:1);
% boxplot(vpks,groupings_S);
% ylabel('Tidal Volume [L]')
% xlabel('Coherence')
% title('Static Frequency Coherence vs Tidal Volume')

% plot them together

Coh_out = []; %initialize a vector of all A outliers
vpks_out = []; %initialize a vector of all B outliers
tot_idx = []; %initialize a vector of all the indices of all outliers

%%%%% code for plotting outliers seperately 
for i = 1:10 %find outliers
    vpkind = find(groupings_D == i);%find indices for each bin
    subvpk = vpks(vpkind);
    vpks_out_idx = isoutlier(subvpk,'quartiles');
    if isempty(subvpk(vpks_out_idx))
        continue
    end
    Coh = coherences_D(vlocs);
    
    out_idx = vpkind(vpks_out_idx);
    vpks_out = [vpks_out vpks(out_idx)];
    Coh_out = [Coh_out Coh(out_idx)];
    tot_idx = [tot_idx out_idx];
end


% pull outlier data out of the main scatter plot data
    Coh_wo_out = coherences_D(vlocs);
    vpks_wo_out = vpks;
    Coh_wo_out(tot_idx) = [];
    vpks_wo_out(tot_idx) = [];
    
size(vpks)
length(tot_idx)
size(Coh_wo_out)
size(vpks_wo_out)

%%%%%

figure('Position',[500 50 500 300]);
hold on
b = boxchart(groupings_D*0.1-0.05,vpks,'MarkerStyle','none'); %'none' suppresses outlier plotting
b.BoxWidth = 0.09;
scatter(Coh_wo_out,vpks_wo_out,[],[0 0.4470 0.7410])
scatter(Coh_out,vpks_out,[],[0.8500 0.3250 0.0980],'*')
xlim([0,1])
ylabel('Tidal Volume [L]')
xlabel('Coherence')
title('Dynamic Frequency Coherence vs Tidal Volume')
hold off



end

