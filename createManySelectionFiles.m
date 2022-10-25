%createManySelectionFiles

%loops through *all* the experiment files available so you can select out
%blocks of relevant data

%Which experiments would you like to look at? Set them as array "exp"
%Which blocks would you like to look at? Set them as array "b"

clear all
close all

load('fileBlockInfo_C1.mat') %a cell array of all the dates and block ranges for each experimentsom
for exp = 5 %what experiments to look at? 
    A = fileBlockInfo{exp,1}{1,2}(end);

%     for b = 1:A{1}(end)%total number of blocks per experiment %use this code to look at all the blocks in that experiemnt 
      for  b = [7] %what blocks within each experiment would you like to look at?
        %%%%%%%%% Variables to specify %%%%%%%%
        experiment = exp %which experiment [1:9]
        Block = b %desire block of data to load 
        
        save('targetBlock.mat','experiment','Block')
        %%%%%%%%% Functions that will run %%%%%

        [date,x,y,num] = selectBlock(Block,experiment,fileBlockInfo);
        [blocktimes,com,comtext,data,dataend,datastart,firstsampleoffset,rangemax,rangemin,samplerate,tickrate,titles,unittext,unittextmap,filename] = loadPowerLabFile(date,x,y);
        [time,EKGData,SpO2Data,PArtData,CapnoData,FlowData,PActData,PPlData,PAbData,VolAutoData,VolAbsData,PDiData,DateTime] = blockAnalysis(data,num,datastart,dataend,blocktimes,titles);

            app = viewAndSelectData; %open selection app
            waitfor(app); %wait for app to close

    end
end