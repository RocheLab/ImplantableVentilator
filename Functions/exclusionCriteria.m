function [Output_Exc,num_Exc] = exclusionCriteria(Output,lim1_1,lim1_2,lim2_1,lim2_2,lim3_1,lim3_2,lim4_1,lim4_2,lim5_1,lim5_2,lim6_1,lim6_2)
%
%   Takes the [Output] of analyzeSnipPFV.m and runs an exclusion criteria 
%for points that are outside of the bounds expected of the data points.
%Excluded points are most likely the result of numbers misassigned from
%analyzeSnipPFV. num_Exc is the number of points excluded  by this function
%%% Output content is [bound1 bound2 Fpkloc Fpkval Ppkloc Pminloc Vpkloc Vpkval Vminloc F_P P_V F_V F_Pmin Pmin_V Pmin_Vmin  P_Vmin F_Vmin ];
%%%                     1       2       3       4      5    6       7     8      9      10   11  12  13      14    15         16     17     


%%%%%                       Var. Set #                       
% Output(:,11) % P_V %          1
% Output(:,13) % F_Pmin %       2
% Output(:,14) % Pmin_V %       3
% Output(:,15) % Pmin_Vmin %    4  
% Output(:,16) % P_Vmin %       5
% Output(:,17) % F_Vmin %       6
%%%%%

Output_Exc = Output;
col = size(Output,2);
num_Exc = 0; %initialize



%% 1) PAct peak to Vol peak
% P_V %          1
row_ind1_1 = find(Output_Exc(:,11) < lim1_1); 
Output_Exc(row_ind1_1,:) = NaN(length(row_ind1_1),col); 
num_Exc = num_Exc + length(row_ind1_1);
row_ind1_2 = find(Output_Exc(:,11) > lim1_2); 
Output_Exc(row_ind1_2,:) = NaN(length(row_ind1_2),col); 
num_Exc = num_Exc + length(row_ind1_2);
%% 2) Flow peak to Pmin
% F_Pmin %       2
row_ind2_1 = find(Output_Exc(:,13) < lim2_1); % F_Pmin %       2
Output_Exc(row_ind2_1,:) = NaN(length(row_ind2_1),col); 
num_Exc = num_Exc + length(row_ind2_1);
row_ind2_2 = find(Output_Exc(:,13) > lim2_2); 
Output_Exc(row_ind2_2,:) = NaN(length(row_ind2_2),col); 
num_Exc = num_Exc + length(row_ind2_2);
%% 3) Pmin to Vol peak
% Pmin_V %       3
row_ind3_1 = find(Output_Exc(:,14) < lim3_1); 
Output_Exc(row_ind3_1,:) = NaN(length(row_ind3_1),col); 
num_Exc = num_Exc + length(row_ind3_1);
row_ind3_2 = find(Output_Exc(:,14) > lim3_2); 
Output_Exc(row_ind3_2,:) = NaN(length(row_ind3_2),col); 
num_Exc = num_Exc + length(row_ind3_2);
%% 4) Pmin to Volmin
% Pmin_Vmin %    4 
row_ind4_1 = find(Output_Exc(:,15) < lim4_1); 
Output_Exc(row_ind4_1,:) = NaN(length(row_ind4_1),col); 
num_Exc = num_Exc + length(row_ind4_1);
row_ind4_2 = find(Output_Exc(:,15) > lim4_2); 
Output_Exc(row_ind4_2,:) = NaN(length(row_ind4_2),col); 
num_Exc = num_Exc + length(row_ind4_2);
%% 5) PAct peak to Vmin
% P_Vmin %       5
row_ind5_1 = find(Output_Exc(:,16) < lim5_1); 
Output_Exc(row_ind5_1,:) = NaN(length(row_ind5_1),col); 
num_Exc = num_Exc + length(row_ind5_1);
row_ind5_2 = find(Output_Exc(:,16) > lim5_2); 
Output_Exc(row_ind5_2,:) = NaN(length(row_ind5_2),col); 
num_Exc = num_Exc + length(row_ind5_2);
%% 6) Flow to Volmin
% F_Vmin %       6
row_ind6_1 = find(Output_Exc(:,17) < lim6_1); 
Output_Exc(row_ind6_1,:) = NaN(length(row_ind6_1),col); 
num_Exc = num_Exc + length(row_ind6_1);
row_ind6_2 = find(Output_Exc(:,17) > lim6_2); 
Output_Exc(row_ind6_2,:) = NaN(length(row_ind6_2),col); 
num_Exc = num_Exc + length(row_ind6_2);

end

