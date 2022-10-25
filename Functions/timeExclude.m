function [Time_Exc] = timeExclude(Exc1,Exc2,Output)
%UNTITLED4 Summary of this function goes here
%   Exc1 exclusion from the beginning in number of breaths
%   Exc2 exclusion from the end of the segment in number of breaths
%%% Output content is [bound1 bound2 Fpkloc Fpkval Ppkloc Pminloc Vpkloc Vpkval Vminloc F_P P_V F_V F_Pmin Pmin_V Pmin_Vmin  P_Vmin F_Vmin ];
%%%                     1       2       3       4      5    6       7     8      9      10   11  12  13      14    15         16     17     

Time_Exc = Output;
col = size(Output,2);

rowsfirst = find(~ismissing(Output(:,7)),Exc1,'first');
Time_Exc(rowsfirst,:) = NaN(length(rowsfirst),col); 

rowslast = find(~ismissing(Output(:,7)),Exc2,'last');
Time_Exc(rowslast,:) = NaN(length(rowslast),col); 
end

