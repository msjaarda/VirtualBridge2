% Results2MATStruct
% Takes in the folder name where results are and gives a structure variable
% (often called MAT) which stores the results in a way that is condusive to
% plotting them. You can optionally probe folder to see contents using
% "OutputFolderPeek.m"

clear, clc

% For the revamp, should we do multiple Eq ? No ETotal?
% Be sure to figure out how to do deterministic...
% BaseData.Traffic as column name? Do we really need static...

% After AGB... shall we generalize by including Jammed/Flow or others as
% variables?


% Folder name where results are
Folder_Name = 'AGB2002_VB';

% Structure Name, AGB/MAT/WIM
Struct_Name = 'MAT';
%MAT.(Section).(Config).(Dist).(--PLoc--).(AE).(Loc).(Span/10)

% Ensure file list is succinct
File_List = GetFileList(Folder_Name);
% Read in .mat results variables into a signle OInfo variable
for i = 1:length(File_List)
    load(['Output/' Folder_Name '/' File_List(i).name])
    OInfo(i) = OutInfo;
end

% Clear OutInfo to avoid confusion (we now use complete OInfo)
clear OutInfo

% Name each column of the final matrix/table in the structure
% AGB
%ColumnNames = {'EQ1','EQ2','Eq','GS','GD','CS','CD','DS','DD','DetS','DetD','E'};
% AGBMAT
ColumnNames = {'EQ1','EQ2','Eq','GD','MD','DD','DetD','E'};
% WIM
%ColumnNames = {'EQ1','EQ2','Eq','YClass','YClassOW','YAll','WClass','WClassOW','WAll','E'};

% Add to existing Structure
%load('AGBMATResults.mat')

% Initialize DetFlag
DetFlag = 0;
% Search for the deterministic one within the folder... this will populate
% the Det and the ESIA parts...
for i = 1:length(OInfo)
    if strcmp(OInfo(i).BaseData.AnalysisType,'Det')
        DetFlag = i;
    end
end
if DetFlag == 0
    fprintf('No Deterministic / ESIA Included in the Folder\n')
end

% Step through all output files and put them into Structure Variable
for i = 1:length(OInfo)
    
    if i == DetFlag
        continue
    else
        
        % COLLECT TRAFFIC INFO (UNCHANGED FOR ALL) ---------------------------
        % Get Configuration
        if strcmp(OInfo(i).BaseData.LaneDir{:},'1,2')
            Config = 'Bi'; else Config = 'Mo';
        end
        % Get Distribution
        if strcmp(OInfo(i).BaseData.LaneTrDistr{:},'50,50')
            Dist = 'Split';
        elseif strcmp(OInfo(i).BaseData.LaneTrDistr{:},'96,4')
            Dist = 'Stand';
        elseif strcmp(OInfo(i).BaseData.LaneTrDistr{:},'100,0')
            Dist = 'ExSlow'; else Dist = 'ExFast';
        end
        % Get Location
        %Loc = OInfo(i).BaseData.Traffic; % For the future!
        if OInfo(i).BaseData.TrRate == 0.29
            Loc = 'GD';
        elseif OInfo(i).BaseData.TrRate == 0.14
            Loc = 'MD';
        elseif OInfo(i).BaseData.TrRate == 0.07
            Loc = 'DD'; else Loc = 'DetD';
        end
        
        % COLLECT IL INFO (MULTIPLE) -----------------------------------------
        for k = 1:length(OInfo(i).ILData)
            
            Temp = OInfo(i).ILData(k).Name; Temp2 = strfind(Temp,'S'); Temp3 = strfind(Temp,'.');
            % Get Span from InfName
            Span = str2num(Temp(Temp2(end)+1:end));
            % Get Action Effect from InfName
            AE = Temp(Temp3(end-1)+1:Temp3(end)-1);
            % Get Section... ideally ILData matches Results struct better...
            Section = Temp(Temp3(1)+4:Temp3(2)-1);
            if strcmp(Section,'Twin')
                SubSect = Temp(Temp3(2)+1:Temp3(3)-1);
                if ~strcmp(SubSect,'Standard')
                    Section = [Section SubSect(1:3)];
                end
            elseif strcmp(Section,'Slab')
                SubSect = Temp(Temp3(2)+1:Temp3(3)-1);
                if strcmp(SubSect,'Fixed') || strcmp(SubSect,'Pinned')
                    Section = [Section SubSect(1:3)];
                else
                    Section = [Section 'Semi'];
                end
            end
            
            % Treat each case individually... to diificult to generalize
            if contains(Section,'Box') || contains(Section,'Twin')
                SpanDiv = 10;
                % Try assigning values directly, catch if table is not yet initialized
                try 
                    if height(MAT.(Section).(Config).(Dist).(AE)) < 8
                        MAT.(Section).(Config).(Dist).(AE)(end+1,:) = array2table(NaN(1,length(ColumnNames)),'VariableNames',ColumnNames);
                    end
                catch MAT.(Section).(Config).(Dist).(AE) = array2table(NaN(1,length(ColumnNames)),'VariableNames',ColumnNames); end
                MAT.(Section).(Config).(Dist).(AE).(Loc)(Span/SpanDiv) = OInfo(i).ESIM(k);
                
                Indet = find(string({OInfo(DetFlag).ILData.Name}) == OInfo(i).ILData(k).Name);
                MAT.(Section).(Config).(Dist).(AE).DetD(Span/SpanDiv) = OInfo(DetFlag).ESIM(Indet);
                MAT.(Section).(Config).(Dist).(AE).E(Span/SpanDiv) = OInfo(DetFlag).ESIA.Total(Indet);
                MAT.(Section).(Config).(Dist).(AE).Eq(Span/SpanDiv) = OInfo(DetFlag).ESIA.Eq(Indet);
                MAT.(Section).(Config).(Dist).(AE).EQ1(Span/SpanDiv) = OInfo(DetFlag).ESIA.EQ(1,Indet);
                MAT.(Section).(Config).(Dist).(AE).EQ2(Span/SpanDiv) = OInfo(DetFlag).ESIA.EQ(2,Indet);
                
            elseif contains(Section,'Multi')
                % Define PLoc
                PLoc = Temp(Temp3(end-2)+1:Temp3(end-1)-1);
                SpanDiv = 10;
                SpanDivShift = 10;
                try  
                    if height(MAT.(Section).(Config).(Dist).(PLoc).(AE)) < 2
                        MAT.(Section).(Config).(Dist).(PLoc).(AE)(end+1,:) = array2table(NaN(1,length(ColumnNames)),'VariableNames',ColumnNames);
                    end
                catch MAT.(Section).(Config).(Dist).(PLoc).(AE) = array2table(NaN(1,length(ColumnNames)),'VariableNames',ColumnNames); end
                MAT.(Section).(Config).(Dist).(PLoc).(AE).(Loc)((Span-SpanDivShift)/SpanDiv) = OInfo(i).ESIM(k);
                
                Indet = find(string({OInfo(DetFlag).ILData.Name}) == OInfo(i).ILData(k).Name);
                MAT.(Section).(Config).(Dist).(PLoc).(AE).DetD((Span-SpanDivShift)/SpanDiv) = OInfo(DetFlag).ESIM(Indet);
                MAT.(Section).(Config).(Dist).(PLoc).(AE).E((Span-SpanDivShift)/SpanDiv) = OInfo(DetFlag).ESIA.Total(Indet);
                MAT.(Section).(Config).(Dist).(PLoc).(AE).Eq((Span-SpanDivShift)/SpanDiv) = OInfo(DetFlag).ESIA.Eq(Indet);
                MAT.(Section).(Config).(Dist).(PLoc).(AE).EQ1((Span-SpanDivShift)/SpanDiv) = OInfo(DetFlag).ESIA.EQ(1,Indet);
                MAT.(Section).(Config).(Dist).(PLoc).(AE).EQ2((Span-SpanDivShift)/SpanDiv) = OInfo(DetFlag).ESIA.EQ(2,Indet);
                
            else % Slab
                % Define PLoc
                PLoc = Temp(Temp3(end-2)+1:Temp3(end-1)-1);
                SpanDiv = 5;
                try 
                    if height(MAT.(Section).(Config).(Dist).(PLoc).(AE)) < 6
                        MAT.(Section).(Config).(Dist).(PLoc).(AE)(end+1,:) = array2table(NaN(1,length(ColumnNames)),'VariableNames',ColumnNames);
                    end
                catch MAT.(Section).(Config).(Dist).(PLoc).(AE) = array2table(NaN(1,length(ColumnNames)),'VariableNames',ColumnNames); end
                MAT.(Section).(Config).(Dist).(PLoc).(AE).(Loc)(Span/SpanDiv) = OInfo(i).ESIM(k);
              
                Indet = find(string({OInfo(DetFlag).ILData.Name}) == OInfo(i).ILData(k).Name);
                MAT.(Section).(Config).(Dist).(PLoc).(AE).DetD(Span/SpanDiv) = OInfo(DetFlag).ESIM(Indet);
                
                MAT.(Section).(Config).(Dist).(PLoc).(AE).E(Span/SpanDiv) = OInfo(DetFlag).ESIA.Total(Indet);
                MAT.(Section).(Config).(Dist).(PLoc).(AE).Eq(Span/SpanDiv) = OInfo(DetFlag).ESIA.Eq(Indet);
                MAT.(Section).(Config).(Dist).(PLoc).(AE).EQ1(Span/SpanDiv) = OInfo(DetFlag).ESIA.EQ(1,Indet);
                MAT.(Section).(Config).(Dist).(PLoc).(AE).EQ2(Span/SpanDiv) = OInfo(DetFlag).ESIA.EQ(2,Indet);
            end
        end
    end
end

% User may choose to save the file (they do so manually)