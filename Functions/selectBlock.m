function [date,x,y,num] = selectBlock(Block,experiment,fileBlockInfo)
%selectBlock selects the correct export file to load for a desired block
%from specified experiment day
%   experiment = integer dating which experimental day to load from [1:8]
%   Block = integer indicating the desired block to load
%   date = date of the experiment
%   to pull up the desired block, it will be in export file x of y, and it will be the num-th block within that file  
    date = fileBlockInfo{experiment,1}{1,1};
    experimentBlockRanges = fileBlockInfo{experiment,1}{1,2};
    y = length(experimentBlockRanges);
    for i = 1:y
        if Block >= experimentBlockRanges{i}(1) && Block <= experimentBlockRanges{i}(2)
            x=i; 
            num = Block - experimentBlockRanges{i}(1) + 1; 
            break
        end
    end
end

