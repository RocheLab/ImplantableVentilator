function [blocktimes,com,comtext,data,dataend,datastart,firstsampleoffset,rangemax,rangemin,samplerate,tickrate,titles,unittext,unittextmap,filename] = loadPowerLabFile(date,x,y)
%loadPowerLabFiles loads appropriate
%   date = string of date, ex, 08252020
%   x and y are integers specifying which x of y export file for a specific
%   date
    filename = [date,'_C1 PowerLab Data ',int2str(x),'of',int2str(y),'.mat']; 
    load(filename)
end

