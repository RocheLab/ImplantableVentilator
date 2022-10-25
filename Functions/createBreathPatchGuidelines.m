function [v,f] = createBreathPatchGuidelines(Vbnds,y1,y2,tend)
%createBreathPatchGuidelines creates v and f for patch function to
%highlight every other breath bound
%   outputs v,f will be fed into the patch function
%   tend is the last element of the time meatrix that will be plotted
%   y1 and y2 are the y values to plot the patches as

if (-1)^(length(Vbnds)) == 1 %if length is even
    matlength = 2*length(Vbnds);
else %if odd
    matlength = 2*length(Vbnds)+2;
    Vbnds = [Vbnds tend];
end

Ys=[y1;y1;y2;y2];
ycol = repmat(Ys,[matlength*0.25,1]);
Vcol=[];
for n = 1:(matlength/4)
    Vcol = [Vcol; Vbnds(2*n-1);Vbnds(2*n);Vbnds(2*n);Vbnds(2*n-1)];
end
length(Vcol)
length(ycol)
v = [Vcol ycol];
f = reshape(1:matlength,[4,matlength/4])';
end

