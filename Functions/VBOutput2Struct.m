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
            VBResults.LaneTrDistr.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(char(OInfo(i).BaseData.Traffic)) = char(OInfo(i).BaseData.LaneTrDistr);
       
        elseif strcmp(OInfo(i).BaseData.AnalysisType,'WIM')
            %VBResults.LaneTrDistr.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(char(OInfo(i).BaseData.Traffic)) = ' , ';
            if OInfo(i).BaseData.SITE == 11, Traffic = 'Uni2L'; Group = SiteGroups.('Uni2L');
            elseif OInfo(i).BaseData.SITE == 111, Traffic = 'Uni3L'; Group = SiteGroups.('Uni3L');
            elseif OInfo(i).BaseData.SITE == 12, Traffic = 'Bi2L'; Group = SiteGroups.('Bi2L');
                if OInfo(i).BaseData.StopSim, Group(Group == 441) = [];
                end
            elseif OInfo(i).BaseData.SITE == 1122, Traffic = 'Bi4L'; Group = SiteGroups.('Bi4L');
            elseif OInfo(i).BaseData.SITE == 110, Traffic = 'LSVAUni2L'; Group = SiteGroups.('LSVAUni2L');
            else
                Traffic = strcat(Sites.SName(Sites.SITE == OInfo(i).BaseData.SITE),num2str(OInfo(i).BaseData.SITE)); Group = [];
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
            if ~isempty(Group)
                % Cycle through each element in Group and filter Max
                
                try
                    if OInfo(i).SimStop
                    
                        BM = {'Daily', 'Weekly', 'Yearly'};             % j
                        ClassType = {'All', 'ClassOW', 'Class'};        % i
                        DistTypes = {'Lognormal'};
                    
                        for z = 1:length(Group)
                            % Plot BlockMax, find Design Values, Ed, using Beta, rather than 99th percentile
                            for r = ia(k):ia(k)+sum(ic == k)-1
                                for p = 1:length(ClassType)
                                    Class = ClassType{p};
                                    BlockM = BM{2};
                                    % Filter Max for Group element
                                    Maxi = OInfo(i).Max(r).(Class).(BlockM).Max(OInfo(i).Max(r).(Class).(BlockM).SITE == Group(z));
                                
                                
                                    %[~,OutInfo.x_values.(Class).(BlockM)(:,r),OutInfo.y_valuespdf.(Class).(BlockM)(:,r),~] = GetBlockMaxFit(Max(i,'Lognormal',BaseData.Plots(g));
                                    %[ECDF,ECDFRank,PPx,PPy,Fity,OutInfo.LNFitR2] = GetLogNormPPP(Max(r).(Class).(BlockM).Max,false);
                                    [~, AQ, ~] = GetBlockMaxEd(Maxi,BlockM,'Lognormal',OInfo(i).ESIA.Total(r),OInfo(i).ESIA.EQ(:,r),OInfo(i).ESIA.Eq(:,r),0.6,0.5);
                                
                                    Traffic = strcat(Sites.SName(Sites.SITE == Group(z)),num2str(Group(z)));
                                
                                    VBResults.AQ.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).(['Stop' Class])(r-ia(k)+1) = AQ;
                                    VBResults.x.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).(['Stop' Class])(r-ia(k)+1) = cellfun(@str2num,ILSplit(r,8));
                                
                                end
                            end
                        end
                    else
                    
                        BM = {'Daily', 'Weekly', 'Yearly'};             % j
                        ClassType = {'All', 'ClassOW', 'Class'};        % i
                        DistTypes = {'Lognormal'};
                    
                        for z = 1:length(Group)
                            % Plot BlockMax, find Design Values, Ed, using Beta, rather than 99th percentile
                            for r = ia(k):ia(k)+sum(ic == k)-1
                                for p = 1:length(ClassType)
                                    Class = ClassType{p};
                                    BlockM = BM{2};
                                    % Filter Max for Group element
                                    Maxi = OInfo(i).Max(r).(Class).(BlockM).Max(OInfo(i).Max(r).(Class).(BlockM).SITE == Group(z));
                                
                                
                                    %[~,OutInfo.x_values.(Class).(BlockM)(:,r),OutInfo.y_valuespdf.(Class).(BlockM)(:,r),~] = GetBlockMaxFit(Max(i,'Lognormal',BaseData.Plots(g));
                                    %[ECDF,ECDFRank,PPx,PPy,Fity,OutInfo.LNFitR2] = GetLogNormPPP(Max(r).(Class).(BlockM).Max,false);
                                    [~, AQ, ~] = GetBlockMaxEd(Maxi,BlockM,'Lognormal',OInfo(i).ESIA.Total(r),OInfo(i).ESIA.EQ(:,r),OInfo(i).ESIA.Eq(:,r),0.6,0.5);
                                
                                    Traffic = strcat(Sites.SName(Sites.SITE == Group(z)),num2str(Group(z)));
                                
                                    VBResults.AQ.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).(Class)(r-ia(k)+1) = AQ;
                                    VBResults.x.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).(Class)(r-ia(k)+1) = cellfun(@str2num,ILSplit(r,8));
                                
                                end
                            end
                        end
                
                    end
                    
                catch
                    
                    BM = {'Daily', 'Weekly', 'Yearly'};             % j
                    ClassType = {'All', 'ClassOW', 'Class'};        % i
                    DistTypes = {'Lognormal'};
                    
                    for z = 1:length(Group)
                        % Plot BlockMax, find Design Values, Ed, using Beta, rather than 99th percentile
                        for r = ia(k):ia(k)+sum(ic == k)-1
                            for p = 1:length(ClassType)
                                Class = ClassType{p};
                                BlockM = BM{2};
                                % Filter Max for Group element
                                Maxi = OInfo(i).Max(r).(Class).(BlockM).Max(OInfo(i).Max(r).(Class).(BlockM).SITE == Group(z));
                                
                                
                                %[~,OutInfo.x_values.(Class).(BlockM)(:,r),OutInfo.y_valuespdf.(Class).(BlockM)(:,r),~] = GetBlockMaxFit(Max(i,'Lognormal',BaseData.Plots(g));
                                %[ECDF,ECDFRank,PPx,PPy,Fity,OutInfo.LNFitR2] = GetLogNormPPP(Max(r).(Class).(BlockM).Max,false);
                                [~, AQ, ~] = GetBlockMaxEd(Maxi,BlockM,'Lognormal',OInfo(i).ESIA.Total(r),OInfo(i).ESIA.EQ(:,r),OInfo(i).ESIA.Eq(:,r),0.6,0.5);
                                
                                Traffic = strcat(Sites.SName(Sites.SITE == Group(z)),num2str(Group(z)));
                                
                                VBResults.AQ.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).(Class)(r-ia(k)+1) = AQ;
                                VBResults.x.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).(Class)(r-ia(k)+1) = cellfun(@str2num,ILSplit(r,8));
                                
                            end
                        end
                    end
                end
            end
            Group = [];
        end
    end
clear ILSplit ILJoin
end

