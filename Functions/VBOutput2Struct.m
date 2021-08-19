function [VBResults] = VBOutput2Struct(Folder_Name)
%VBOUTPUT2STRUCT Takes in Output Folder Name, and loops through, gathering
% results (AQ = 1.1*EdLN/ESIA) and x values (spans) for each InfCase.
% VBResults, the return variable has parts VBResults.AQ and VBResults.x for
% all the branches. VBResults can be easily plotted (hint, use VBTriPlot,
% with a small preparation code)

% Ensure file list is succinct
File_List = GetFileList(Folder_Name);

% Read in .mat results variables into a single OInfo variable
for i = 1:length(File_List)
    load(['Output/' Folder_Name '/' File_List(i).name])
    OInfo(i) = OutInfo;
end

% Clear OutInfo to avoid confusion (we now use complete OInfo)
clear OutInfo

% Loop through OutInfo
for i = 1:length(OInfo)
    
    ILNames = string({OInfo(i).ILData.Name}');
    
    % Loop through ILNames
    for j = 1:length(ILNames)
        ILSplit(j,:) = strsplit(ILNames(j),'.'); 
    end
    % Delete first column (ILLIb)
    ILSplit(:,1) = [];
    % Remove "S"s from spans
    ILSplit(:,8) = extractAfter(ILSplit(:,8),1);
    % Remove "AGB" and "MAT"
    ILSplit(:,1) = extractAfter(ILSplit(:,1),3);
    
    % Loop through ILNames
    for j = 1:length(ILNames)
        ILJoin(j) = strjoin(ILSplit(j,1:7),'.');
    end
    ILJoin = ILJoin';
    % Gather the names and group them into unique ones...
    [C,ia,ic] = unique(ILJoin);
    % Put into structure... with end table including the BaseData.Traffic
    for k = 1:length(ia)
        VBResults.AQ.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(OInfo(i).BaseData.Traffic{:}) = OInfo(i).AQ(ic == k);
        VBResults.x.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(OInfo(i).BaseData.Traffic{:}) = cellfun(@str2num,ILSplit(ic == k,8));
    end
    clear ILSplit
end

end

