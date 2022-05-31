% ProcessOutput ----------------------------------------------------------

clear, clc, close all

% ----- TRAFFIC ANALYSES -----
% VBSim (VBSimInput) OR VBWIM (VBWIMInput)
% Goal is to get OverMax (Sim) MaxEvents (WIM)

% ----- TRAFFIC ANALYSES POST-PROCESSING -----
% Goals of PROCESSOUTPUT
% 0. Remove unnecessary things
% 1. Move from OverMaxT/MaxEvents to BlockMax (Max)             GetBlockMax
% 2. Perform fits of Max & get design values (Eds)              GetFit
% 3. Compare to code values to get partial factors/alphas       VBGetECode
% 4. Organize in easily plottable structures                    GetPlotFormat
%    Replaces Output2Struct

% ----- TRAFFIC ANALYSES PLOTTING -----
% AlphaSummaryPlot(_Q2)
% VBPlotProResults
% VBPlotFatigue

% Create a script that will go into an Output folder, load the contents,
% delete unnecessary things... and re-add stuff we care about.

% Minimum what we require:
%   Name (datestring)
%   BaseData (table)
%   ILData (struct)
%   MaxEvents (table) daily maxes
%       {'DTS', 'SITE', 'MaxLE', 'InfCase', 'm', 'DayRank', 'BrStInd', 'PlatType'}
% OR
%   OverMaxT (table) yearly maxes
%       {'InfCase','SimNum','BatchNum','MaxLE','DLF','BrStInd','MaxDamage'}

% Things we will add:
%   Max (struct)
%   pd, Eds (structs)
%   Ecodes (struct)
%   Alphas (struct)

% This script will work with any output file that contains MaxEvents or
% OverMaxT... will hunt for unneeded things and delete them or archive them
% When MaxEvents doesn't exist, it will work with Max

Folder_Name = 'WIM60tv11';
NewFolder = 'WIM1160tAll';
IncZ = 0; % Line 123-124 modify

% Ensure file list is succinct
File_List = GetFileList(Folder_Name);

% Load even if not WIM, just in case
load('Sites.mat'); LenPrint = []; RamUsed = [];

% Start Progress Bar
u = StartProgBar(length(File_List), 1, 1, 5); tic; st = now;

% Read in .mat results variables into a single OInfo variable
for v = 1:length(File_List)
    load(['Output/' Folder_Name '/' File_List(v).name])
    OInfo(v) = OutInfo;
    OInfo(v).BaseData = BaseDataDefaults(OInfo(v).BaseData);
    % Update progress bar
    user = memory;
    RamUsed = [RamUsed;user.MemUsedMATLAB/(user.MemAvailableAllArrays+user.MemUsedMATLAB)*100];
    LenPrint = VBUpProgBar(st,RamUsed(end),v,LenPrint);
end
clear OutInfo

% Hunt for the following to delete:
% .pd, .Max, .x_values, .y_valuespdf, .EdLN, ESIA, .AQ, .SimStop (use BaseData)

for v = 1:length(OInfo)
    fields = fieldnames(OInfo(v));
    if any(contains(fields,'SimStop'))
        if OInfo(v).SimStop == 1
            OInfo(v).BaseData.StopSim = 1;
        else
            OInfo(v).BaseData.StopSim = 0;
        end
    end
end
if any(contains(fields,'SimStop'))
    OInfo = rmfield(OInfo,'SimStop');
end
if any(contains(fields,'y_valuespdf'))
    OInfo = rmfield(OInfo,'y_valuespdf');
    fields = fieldnames(OInfo(v));
end
if any(contains(fields,'pd'))
    OInfo = rmfield(OInfo,'pd');
end
if any(contains(fields,'x_values'))
    OInfo = rmfield(OInfo,'x_values');
end
if any(contains(fields,'EdLN'))
    OInfo = rmfield(OInfo,'EdLN');
end
if any(contains(fields,'AQ'))
    OInfo = rmfield(OInfo,'AQ');
end
if any(contains(fields,'ESIA'))
    OInfo = rmfield(OInfo,'ESIA');
end
if any(contains(fields,'MaxEvents'))
    if any(matches(fields,'max','IgnoreCase',true))
        OInfo = rmfield(OInfo,'Max');
    end
end
fields = fieldnames(OInfo(v));

% If we have MaxEvents then perform GetBlockMax
% Otherwise, 2 scenarios
% 1. We are doing WIM, but only have Max Available (skip to next step)
%    This is not ideal, but OK (GetBlockMax is better than predecessors
%    that were probably used to create Max... elimantes low values from
%    partial blocks for example.
% 2. We are doing SIM, and therefore only have OverMax

% GetBlockMax and GetFit
BlockMax = {'Daily', 'Weekly', 'Monthly', 'Yearly'};        % j
ClassTypes = {'ClassOW'}; %{'All', 'ClassOW', 'Class'};     % i
%DistTypes = {'All'};                                        % k
DistTypes = {'NormalLM', 'LognormalLM', 'LognormalTF', 'gev', 'gevGumbel'}; % For the 60t analyses
if strcmp(OInfo(1).BaseData.AnalysisType,'WIM')
    % Start Progress Bar
    u = StartProgBar(length(File_List), 1, 2, 5); tic; st = now;
    for v = 1:length(OInfo)
        % GetBlockMax can handle however MaxEvents is structured...
        % no need to loop ILs, BlockMaxs, or ClassTypes
        if any(contains(fields,'MaxEvents'))
            OInfo(v).Max = GetBlockMax(OInfo(v).MaxEvents,ClassTypes,BlockMax);
        else
            BlockMax = fieldnames(OInfo(v).Max(1).All);        % j
        end
        % GetFit, only the other hand we need to give specific Data and BlockM
        % for this we loop between ILs, BlockMax, and ClassTypes
        for r = 1:length(OInfo(v).ILData)
            for j = 1:length(BlockMax)
                BM = BlockMax{j};
                for i = 1:length(ClassTypes)
                    CT = ClassTypes{i};
                    OInfo(v).pd(r).(CT).(BM) = GetFit(OInfo(v).Max(r).(CT).(BM).Max,BM,DistTypes,0,IncZ);
                end
            end
        end
        % Update progress bar
        user = memory;
        RamUsed = [RamUsed;user.MemUsedMATLAB/(user.MemAvailableAllArrays+user.MemUsedMATLAB)*100];
        LenPrint = VBUpProgBar(st,RamUsed(end),v,LenPrint);
    end
end

% For Sim
if strcmp(OInfo(1).BaseData.AnalysisType,'Sim')
    % Start Progress Bar
    u = StartProgBar(length(File_List), 1, 2, 5); tic; st = now;
    BlockMax = OInfo(v).BaseData.Period; BM = BlockMax{1};
    ClassTypes = {'Class'}; CT = ClassTypes{1};
    for v = 1:length(OInfo)
        for r = 1:length(OInfo(v).ILData)
            for j = 1:length(BlockMax)
                for i = 1:length(ClassTypes)
                    OInfo(v).pd(r).(CT).(BM) = GetFit(OInfo(v).OverMax(:,r),BM,DistTypes,0,IncZ);
                end
            end
        end
        % Update progress bar
        user = memory;
        RamUsed = [RamUsed;user.MemUsedMATLAB/(user.MemAvailableAllArrays+user.MemUsedMATLAB)*100];
        LenPrint = VBUpProgBar(st,RamUsed(end),v,LenPrint);
    end
end

Gamma = 1.5;
AQ1 = 0.7; AQ2 = 0.5;
% VBGetECode and GetAlphas
% Start Progress Bar
u = StartProgBar(length(File_List), 1, 3, 5); tic; st = now;
for v = 1:length(OInfo)
    OInfo(v).E = VBGetECode(OInfo(v).ILData,OInfo(v).BaseData.ILRes);
    for r = 1:length(OInfo(v).ILData)
        for j = 1:length(BlockMax)
            BM = BlockMax{j};
            for i = 1:length(ClassTypes)
                CT = ClassTypes{i};
                Ed = OInfo(v).pd(r).(CT).(BM).(OInfo(v).pd(r).(CT).(BM).Best).Ed;
                OInfo(v).AQ.(CT).(BM)(r) = Ed./OInfo(v).E(r).SIA.Total;
                OInfo(v).Aq.(CT).(BM)(r) = ((Ed/1.5)-AQ1*OInfo(v).E(r).SIA.EQ(1)-AQ2*OInfo(v).E(r).SIA.EQ(2))./(sum(OInfo(v).E(r).SIA.Eq));
            end
        end
    end
    % Update progress bar
    user = memory;
    RamUsed = [RamUsed;user.MemUsedMATLAB/(user.MemAvailableAllArrays+user.MemUsedMATLAB)*100];
    LenPrint = VBUpProgBar(st,RamUsed(end),v,LenPrint);
end

% Repeat some of the above for each site... save as OInfo(v).SNameSITE
if strcmp(OInfo(1).BaseData.AnalysisType,'WIM')
    % Start Progress Bar
    u = StartProgBar(length(File_List), 1, 4, 5); tic; st = now;
    BlockMax = {'Weekly'};        % j
    for v = 1:length(OInfo)
        if isempty(strcat(Sites.SName(Sites.SITE == OInfo(v).BaseData.SITE),num2str(OInfo(v).BaseData.SITE)))
            Sitex = VBGetSiteSet(OInfo(v).BaseData.SITE,OInfo(v).BaseData.LightVehs,0,OInfo(v).BaseData.Country);
            for z = 1:length(Sitex)
                Traffic = strcat(Sites.SName(Sites.SITE == Sitex(z)),num2str(Sitex(z)));
                % GetFit for individual SITE
                for r = 1:length(OInfo(v).ILData)
                    for j = 1:length(BlockMax)
                        BM = BlockMax{j};
                        for i = 1:length(ClassTypes)
                            CT = ClassTypes{i};
                            Maxi = OInfo(v).Max(r).(CT).(BM).Max(OInfo(v).Max(r).(CT).(BM).SITE == Sitex(z));
                            OInfo(v).(Traffic).pd(r).(CT).(BM) = GetFit(Maxi,BM,DistTypes,0,IncZ);
                            Ed = OInfo(v).(Traffic).pd(r).(CT).(BM).(OInfo(v).(Traffic).pd(r).(CT).(BM).Best).Ed;
                            OInfo(v).(Traffic).AQ.(CT).(BM)(r) = Ed./OInfo(v).E(r).SIA.Total;
                            OInfo(v).(Traffic).Aq.(CT).(BM)(r) = ((Ed/1.5)-AQ1*OInfo(v).E(r).SIA.EQ(1)-AQ2*OInfo(v).E(r).SIA.EQ(2))./(sum(OInfo(v).E(r).SIA.Eq));
                        end
                    end
                end
            end
        end
        % Update progress bar
        user = memory;
        RamUsed = [RamUsed;user.MemUsedMATLAB/(user.MemAvailableAllArrays+user.MemUsedMATLAB)*100];
        LenPrint = VBUpProgBar(st,RamUsed(end),v,LenPrint);
    end
end

% check if NewFolder folder exist, if not create one
Dir_List = dir('Output/');
Folder_List = {Dir_List.name}';
if sum(strcmp(Folder_List,NewFolder))>=1
else
   mkdir(append('Output/',NewFolder));
end

% Save
for v = 1:length(OInfo)
    OutInfo = OInfo(v);
    %save(['Output' char(OutInfo.BaseData.Folder) '/' OutInfo.Name], 'OutInfo')
    save(['Output'  '/' NewFolder '/' OutInfo.Name], 'OutInfo')
end

% GetPlotFormat
%BlockM = {'Weekly'};
BlockM = BlockMax; %{'Daily', 'Weekly', 'Monthly', 'Yearly'};
ClassTypes = {'All', 'ClassOW', 'Class'}; 

% Make a table
for v = 1:length(BlockM)
VBResults.(BlockM{v}) = array2table(zeros(0,12), 'VariableNames',{'Type','SubType','Width','Layout','Support','Trans','AE','Traffic','Span','All','ClassOW','Class'});
% Delete at the end if they are empty?
VBResults.SS.(BlockM{v}) = array2table(zeros(0,12), 'VariableNames',{'Type','SubType','Width','Layout','Support','Trans','AE','Traffic','Span','All','ClassOW','Class'});
end

% Start Progress Bar
u = StartProgBar(length(File_List), 1, 5, 5); tic; st = now;

% Loop through OutInfo
for v = 1:length(OInfo)
    
    ILNames = string({OInfo(v).ILData.Name}');
    clear ILSplit ILJoin
    
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
        if strcmp(OInfo(v).BaseData.AnalysisType,'Sim')
            Traffic = OInfo(v).BaseData.Traffic{:};
            BlockMax = OInfo(v).BaseData.Period; BlockM = BlockMax{1};
            ClassTypes = {'Class'}; CT = ClassTypes{1};
            VBResults.AQ.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic) = OInfo(v).AQ.(CT).(BlockM)(ic == k);
            VBResults.x.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(OInfo(v).BaseData.Traffic{:}) = cellfun(@str2num,ILSplit(ic == k,8));
            VBResults.LaneTrDistr.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(OInfo(v).BaseData.Traffic{:}) = OInfo(v).BaseData.LaneTrDistr{:};          
        elseif strcmp(OInfo(v).BaseData.AnalysisType,'WIM')
            % Similar to VBGetSiteSet
            Traffic = strcat(Sites.SName(Sites.SITE == OInfo(v).BaseData.SITE),num2str(OInfo(v).BaseData.SITE)); Sitex = [];
            if isempty(Traffic)
                Sitex = VBGetSiteSet(OInfo(v).BaseData.SITE,OInfo(v).BaseData.LightVehs,0,OInfo(v).BaseData.Country);
                Traffic = ConvertLayoutName(OInfo(v).BaseData.SITE);
            end
            for n = 1:length(BlockM)
                for m = 1:length(ClassTypes)
            CT = ClassTypes{m};
            T1 = array2table(ILSplit(ic == k,1:7),'VariableNames',{'Type','SubType','Width','Layout','Support','Trans','AE'});
            try T8 = array2table(repmat(string(Traffic),length(OInfo(v).AQ.(CT).(BlockM{n})(ic == k)),1),'VariableNames',{'Traffic'}); catch end
            T9 = array2table(cellfun(@str2num,ILSplit(ic == k,8)),'VariableNames',{'Span'});
            %try T10 = array2table(OInfo(v).AQ.All.(BlockM{n})(ic == k)','VariableNames',{'All'}); catch end
            %T10.(append('table',int2str(n))) = array2table(OInfo(v).AQ.ClassOW.(BlockM{n})(ic == k)','VariableNames',{CT});
            try T10(:,m) = array2table(OInfo(v).AQ.(CT).(BlockM{n})(ic == k)');
            catch T10(:,m) = array2table(ic(ic == k).*0);
            end
            T10.Properties.VariableNames{m} = CT;
            %try T12 = array2table(OInfo(v).AQ.Class.(BlockM{n})(ic == k)','VariableNames',{'Class'}); catch end
                end
                if OInfo(v).BaseData.StopSim
                    for m = 1:length(ClassTypes)
                    CT = ClassTypes{m};
                    try VBResults.AQ.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).(['Stop' CT]) = OInfo(v).AQ.(CT).(BlockM{n})(ic == k);
                    catch VBResults.AQ.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).(['Stop' CT]) = (ic(ic == k).*0)';
                    end
                    VBResults.x.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).(['Stop' CT]) = cellfun(@str2num,ILSplit(ic == k,8));
                    end
                    VBResults.SS.(BlockM{n}) = [VBResults.(BlockM{n}); T1 T8 T9 T10];
                elseif OInfo(v).BaseData.Plat
                    for m = 1:length(ClassTypes)
                    CT = ClassTypes{m};
                    try VBResults.P.AQ.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).(strcat('Size',num2str(OInfo(v).BaseData.PlatSize))).(strcat('Rate',num2str(10*OInfo(v).BaseData.PlatRate))).(strcat('FolDist',num2str(10*OInfo(v).BaseData.PlatFolDist))).(CT) = OInfo(v).AQ.(CT).(BlockM{n})(ic == k);
                    catch VBResults.P.AQ.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).(strcat('Size',num2str(OInfo(v).BaseData.PlatSize))).(strcat('Rate',num2str(10*OInfo(v).BaseData.PlatRate))).(strcat('FolDist',num2str(10*OInfo(v).BaseData.PlatFolDist))).(CT) = (ic(ic == k).*0)';
                    end
                    VBResults.P.x.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).(strcat('Size',num2str(OInfo(v).BaseData.PlatSize))).(strcat('Rate',num2str(10*OInfo(v).BaseData.PlatRate))).(strcat('FolDist',num2str(10*OInfo(v).BaseData.PlatFolDist))).(CT) = cellfun(@str2num,ILSplit(ic == k,8));
                    end
                else
                    for m = 1:length(ClassTypes)
                    CT = ClassTypes{m};
                    try VBResults.AQ.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).(CT) = OInfo(v).AQ.(CT).(BlockM{n})(ic == k);
                    catch VBResults.AQ.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).(CT) = (ic(ic == k).*0)';
                    end
                    VBResults.x.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).(CT) = cellfun(@str2num,ILSplit(ic == k,8));
                    end
                    VBResults.(BlockM{n}) = [VBResults.(BlockM{n}); T1 T8 T9 T10];
                end
                clear T1 T8 T9 T10
            end
            for z = 1:length(Sitex)
                for p = 1:length(ClassTypes)
                    CT = ClassTypes{p};
                    Traffic = ConvertLayoutName(Sitex(z));
                    if OInfo(v).BaseData.StopSim
                        try VBResults.AQ.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).(['Stop' CT]) = OInfo(v).(Traffic).AQ.(CT).Weekly(ic == k);
                        catch VBResults.AQ.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).(['Stop' CT]) = (ic(ic == k).*0)';
                        end
                        VBResults.x.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).(['Stop' CT]) = cellfun(@str2num,ILSplit(ic == k,8));
                        %VBResults.SS.(BlockM{n}) = [VBResults.(BlockM{n}); T1 T8 T9 T10 T11 T12];
                    else
                        try VBResults.AQ.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).(CT) = OInfo(v).(Traffic).AQ.(CT).Weekly(ic == k);
                        catch VBResults.AQ.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).(CT) = (ic(ic == k).*0)';
                        end
                        VBResults.x.(ILSplit(ia(k),1)).(ILSplit(ia(k),2)).(ILSplit(ia(k),3)).(ILSplit(ia(k),4)).(ILSplit(ia(k),5)).(ILSplit(ia(k),6)).(ILSplit(ia(k),7)).(Traffic).(CT) = cellfun(@str2num,ILSplit(ic == k,8));
                        %VBResults.(BlockM{n}) = [VBResults.(BlockM{n}); T1 T8 T9 T10 T11 T12];
                    end
                end
            end
        end % Sim vs WIM
    end % k  Unique AEs
    % Update progress bar
    user = memory;
    RamUsed = [RamUsed;user.MemUsedMATLAB/(user.MemAvailableAllArrays+user.MemUsedMATLAB)*100];
    LenPrint = VBUpProgBar(st,RamUsed(end),v,LenPrint);
end % v OInfo


%if isempty(VBResults.(BlockM))
%    VBResults.(BlockM) = [];
%end
%if isempty(VBResults.SS.(BlockM))
%    VBResults.SS.(BlockM) = [];
%end

% User Save

% The end result of VBOutput2Struct is a var, WIM/AGB/MAT

    