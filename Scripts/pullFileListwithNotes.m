%pullFileListwithNotes
%pulls list of selection files, after generating them from
%createManySelectionFiles, extracts notes from the files and amends them to
%the table
clear all 
close all

list = struct2table(dir('Exp *.mat')); 

FileNames = list.name;

for ii=[1:size(FileNames,1)]
    load(FileNames{ii})
    Tnote{ii,1} = (notes);
end

T=cell2table(Tnote);
list = addvars(list,T.Tnote,'Before','folder');
list=list;
save list