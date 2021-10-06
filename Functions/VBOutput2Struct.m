function [VBResults] = VBOutput2Struct(Folder_Name)
%VBOUTPUT2STRUCT Takes in Output Folder Name, and loops through, gathering
% results (AQ = 1.1*EdLN/ESIA) and x values (spans) for each InfCase.
% VBResults, the return variable has parts VBResults.AQ and VBResults.x for
% all the branches. VBResults can be easily plotted (hint, use VBTriPlot,
% with a small preparation code)

% Ensure file list is succinct
File_List = GetFileList(Folder_Name);

% Load even if not WIM, just in case
load('Sites.mat')
load('SiteGroups.mat')

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
    %ILSplit(:,1) = extractAfter(ILSplit(:,1),3);
    
    % Loop through ILNames
    for j = 1:length(ILNames)
        ILJoin(j) = strjoin(ILSplit(j,1:7),'.');
    end
    ILJoin = ILJoin';
    % Gather the names and group them into unique ones...
    [~,ia,ic] = unique(ILJoin);
    % Put into structure... with end table including the BaseData.Traffic
    for k = 1:length(ia)
        if strcmp(OInfo(i).BaseData.AnalysisType,'Sim')
            VBResults.AQ.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(OInfo(i).BaseData.Traffic{:}) = OInfo(i).AQ(ic == k);
            VBResults.x.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(OInfo(i).BaseData.Traffic{:}) = cellfun(@str2num,ILSplit(ic == k,8));
        elseif strcmp(OInfo(i).BaseData.AnalysisType,'WIM')
            if OInfo(i).BaseData.SITE == 11, Traffic = SiteGroups.('Uni2L');
            elseif OInfo(i).BaseData.SITE == 111, Traffic = SiteGroups.('Uni3L');
            elseif OInfo(i).BaseData.SITE == 12, Traffic = SiteGroups.('Bi2L');
            elseif OInfo(i).BaseData.SITE == 1122, Traffic = SiteGroups.('Bi4L');
            elseif OInfo(i).BaseData.SITE == 110, Traffic = SiteGroups.('LSVAUni2L');
            else
                Traffic = strcat(Sites.SName(Sites.SITE == OInfo(i).BaseData.SITE),num2str(OInfo(i).BaseData.SITE));
            end
            try
                if OInfo(i).SimStop
                    VBResults.AQ.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).StopAll = OInfo(i).AQ.All.Weekly(ic == k);
                    VBResults.AQ.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).StopClassOW = OInfo(i).AQ.ClassOW.Weekly(ic == k);
                    VBResults.AQ.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).StopClass = OInfo(i).AQ.Class.Weekly(ic == k);
                    VBResults.x.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).StopAll = cellfun(@str2num,ILSplit(ic == k,8));
                    VBResults.x.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).StopClassOW = cellfun(@str2num,ILSplit(ic == k,8));
                    VBResults.x.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).StopClass = cellfun(@str2num,ILSplit(ic == k,8));
                    
                else
                    
                    VBResults.AQ.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).All = OInfo(i).AQ.All.Weekly(ic == k);
                    VBResults.AQ.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).ClassOW = OInfo(i).AQ.ClassOW.Weekly(ic == k);
                    VBResults.AQ.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).Class = OInfo(i).AQ.Class.Weekly(ic == k);
                    VBResults.x.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).All = cellfun(@str2num,ILSplit(ic == k,8));
                    VBResults.x.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).ClassOW = cellfun(@str2num,ILSplit(ic == k,8));
                    VBResults.x.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).Class = cellfun(@str2num,ILSplit(ic == k,8));
                end
            catch
                VBResults.AQ.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).All = OInfo(i).AQ.All.Weekly(ic == k);
                VBResults.AQ.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).ClassOW = OInfo(i).AQ.ClassOW.Weekly(ic == k);
                VBResults.AQ.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).Class = OInfo(i).AQ.Class.Weekly(ic == k);
                VBResults.x.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).All = cellfun(@str2num,ILSplit(ic == k,8));
                VBResults.x.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).ClassOW = cellfun(@str2num,ILSplit(ic == k,8));
                VBResults.x.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).Class = cellfun(@str2num,ILSplit(ic == k,8));
                
            end
        end
        
    end
clear ILSplit
end

