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

Folder_Name = 'WIMOct5';
NewFolder = 'WIMOct5Redopr';
IncZ = 0; % Line 123-124 modify

% Ensure file list is succinct
File_List = GetFileList(Folder_Name);

% Load even if not WIM, just in case
load('Sites.mat'); LenPrint = []; RamUsed = [];

% Start Progress Bar
u = StartProgBar(length(File_List), 1, 1, 4); tic; st = now;

% Read in .mat results variables into a single OInfo variable
for v = 1:length(File_List)
    load(['Output/' Folder_Name '/' File_List(v).name])
    OInfo(v) = OutInfo;
    OInfo(v).BaseData = BaseDataDefaults(OInfo(v).BaseData);
    % TEMP
    OInfo(v).MaxEvents(OInfo(v).MaxEvents.SITE == 409,:) = [];
    % Update progress bar
    try
        user = memory;
        RamUsed = [RamUsed;user.MemUsedMATLAB/(user.MemAvailableAllArrays+user.MemUsedMATLAB)*100];
        LenPrint = VBUpProgBar(st,RamUsed(end),v,LenPrint);
    catch
        LenPrint = VBUpProgBar(st,1,v,LenPrint);
    end
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
BlockMax = {'Weekly'};        % j
ClassTypes = {'All', 'ClassOW', 'Class'}; %{'ClassOW'}; %{'All', 'ClassOW', 'Class'};     % i
DistTypes = {'All'};                                        % k
%DistTypes = {'NormalLM', 'LognormalLM', 'LognormalTF', 'gev', 'gevGumbel'}; % For the 60t analyses
if strcmp(OInfo(1).BaseData.AnalysisType,'WIM')
    % Start Progress Bar
    u = StartProgBar(length(File_List), 1, 2, 4); tic; st = now;
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
                    OInfo(v).pd(r).(CT).(BM) = GetFit(OInfo(v).Max(r).(CT).(BM).Max,BM,DistTypes,0,0);
                end
            end
        end
        % Update progress bar
        try
            user = memory;
            RamUsed = [RamUsed;user.MemUsedMATLAB/(user.MemAvailableAllArrays+user.MemUsedMATLAB)*100];
            LenPrint = VBUpProgBar(st,RamUsed(end),v,LenPrint);
        catch
            LenPrint = VBUpProgBar(st,1,v,LenPrint);
        end
    end
end

% For Sim
if strcmp(OInfo(1).BaseData.AnalysisType,'Sim')
    % Start Progress Bar
    u = StartProgBar(length(File_List), 1, 2, 4); tic; st = now;
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
        try
            user = memory;
            RamUsed = [RamUsed;user.MemUsedMATLAB/(user.MemAvailableAllArrays+user.MemUsedMATLAB)*100];
            LenPrint = VBUpProgBar(st,RamUsed(end),v,LenPrint);
        catch
            LenPrint = VBUpProgBar(st,1,v,LenPrint);
        end
    end
end

Gamma = 1.5;
AQ1 = 0.7; AQ2 = 0.5;
% VBGetECode and GetAlphas
% Start Progress Bar
u = StartProgBar(length(File_List), 1, 3, 4); tic; st = now;
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
    try
        user = memory;
        RamUsed = [RamUsed;user.MemUsedMATLAB/(user.MemAvailableAllArrays+user.MemUsedMATLAB)*100];
        LenPrint = VBUpProgBar(st,RamUsed(end),v,LenPrint);
    catch
        LenPrint = VBUpProgBar(st,1,v,LenPrint);
    end
end

% Repeat some of the above for each site... save as OInfo(v).SNameSITE
if strcmp(OInfo(1).BaseData.AnalysisType,'WIM')
    % Start Progress Bar
    u = StartProgBar(length(File_List), 1, 4, 4); tic; st = now;
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
                            if ~isempty(Maxi)
                                OInfo(v).(Traffic).pd(r).(CT).(BM) = GetFit(Maxi,BM,DistTypes,0,IncZ);
                                Ed = OInfo(v).(Traffic).pd(r).(CT).(BM).(OInfo(v).(Traffic).pd(r).(CT).(BM).Best).Ed;
                                OInfo(v).(Traffic).AQ.(CT).(BM)(r) = Ed./OInfo(v).E(r).SIA.Total;
                                OInfo(v).(Traffic).Aq.(CT).(BM)(r) = ((Ed/1.5)-AQ1*OInfo(v).E(r).SIA.EQ(1)-AQ2*OInfo(v).E(r).SIA.EQ(2))./(sum(OInfo(v).E(r).SIA.Eq));
                        
                            end
                        end
                    end
                end
            end
        end
        % Update progress bar
        try
            user = memory;
            RamUsed = [RamUsed;user.MemUsedMATLAB/(user.MemAvailableAllArrays+user.MemUsedMATLAB)*100];
            LenPrint = VBUpProgBar(st,RamUsed(end),v,LenPrint);
        catch
            LenPrint = VBUpProgBar(st,1,v,LenPrint);
        end
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
    save(['Output'  '/' NewFolder '/' OutInfo.Name], 'OutInfo','-v7.3')
end

    