clear, clc, close all, load('Sites.mat')
% Script to analyse individual results from AlphaSummaryPlot

%% INPUT
% Select the Infl Line to inspect
InfLine = 'Box.Stand.Wid12.Uni.Simp.p0.V.S40';
OutputFolder = 'WIMMatt020622Output';
BlockM = 'Weekly';
Class = 'ClassOW';
DistTypes = 'All'; % fitting
NumAnalyses = 5;

%% CODE
Dir_List = dir('Output');
File_List = {Dir_List.name}';

% check if WIM folder exist
if sum(strcmp(File_List,OutputFolder))>=1
else
   error(append('No ',OutputFolder,' Folder!!'));
end

Dir_List = dir(append('Output/',OutputFolder));
File_List = {Dir_List(:).name}';

%cleaning file list
File_List = File_List(~strcmp(File_List,'.')&~strcmp(File_List,'..')&contains(File_List,'.mat'));
File_List = erase(File_List,'.mat');

for i=1:height(File_List)
    
    load(append('Output/',OutputFolder,'/',File_List{i},'.mat'));
    
    if OutInfo.BaseData.StopSim == 0
    ILDataNames = {OutInfo.ILData.Name}';
    if sum(contains(ILDataNames,InfLine))>=1
        File_List = File_List{i};
        InflPosi = find(contains(ILDataNames,InfLine));
        break
    end
    end
end

Max = OutInfo.Max(InflPosi).(Class).(BlockM);
Data = Max.Max;
Max = sortrows(Max,3,'descend');

pd = GetFit(Data,BlockM,DistTypes,1,1);
%pd = GetFit(Data(~isoutlier(Data,'gesd')),BlockM,DistTypes,1,1);

% Find sites that generate Max cases
%[SitesAct,~,c] = unique(Max.SITE(1:NumAnalyses));
%for i=1:max(c)
%    NumSitesAct(i) = sum(c==i);
%end
%NumSitesAct = NumSitesAct';
SitesAct = Max.SITE(1:NumAnalyses);

% Create table
BaseData = table('Size',[height(SitesAct),10],'VariableTypes',["double","string","double","double","double","double","double","string","string","datetime"]);
BaseData.Properties.VariableNames = {'SITE','ILs','ILRes','RunDyn','Parallel','Apercu','NumAnalyses','AnalysisType','ClassType','Date'};
BaseData.SITE = SitesAct;
BaseData.ILs(:) = InfLine;
BaseData.ILRes(:) = OutInfo.BaseData.ILRes;
BaseData.Apercu(:) = 1;
BaseData.NumAnalyses(:) = 1;
BaseData.AnalysisType(:) = 'WIM';
BaseData.ClassType(:) = Class;
BaseData.Date = Max.DTS(1:NumAnalyses);

% Initialize parpool if necessary and initialize progress bar
u = StartProgBar(height(BaseData), 1, 1, 1); st = now;

% Each row of BaseData represents one analysis
%parfor g = 1:height(BaseData)
for g = 1:height(BaseData)
    
    % Update analysis data for current row of BaseData
    [Num,Lane,ILData,~,~] = VBUpdateData(BaseData(g,:));
    
    % Load File
    load(['WIM/',num2str(BaseData.SITE(g)),'.mat']);
    PDs = Stage2Prune(PDs);
    PDs = PDs(PDs.DTS > (BaseData.Date(g)-1),:);
    PDs = PDs(PDs.DTS < (BaseData.Date(g)+1),:);
    
%      % Find and remove duplicates
%      if Sites.Layout(Sites.SITE == BaseData.SITE(g)) == 11
%             % Get Duplicates
%             PDs = FindDup2(PDs,0,0);
%             % Delete Duplicates - from L1
%             PDs(PDs.Dup & PDs.LANE == 1,:) = [];
%      end
    
    % Get Only the ClassType Specified
    try
        if strcmp(BaseData.ClassType(g),'Class')
            PDs = PDs(PDs.CLASS > 0 & (PDs.CLASS > 50 | PDs.CLASS < 40),:);
        elseif strcmp(BaseData.ClassType(g),'ClassOW')
            PDs = PDs(PDs.CLASS > 0,:);
        end
    catch end
    
    % Convert PDC to AllTrAx - Spacesave at MaxLength
    MaxLength = (max(arrayfun(@(x) size(x.v,1),ILData))-1)*BaseData.ILRes(g);
    [PDs, AllTrAx, TrLineUp] = VBWIMtoAllTrAx(PDs,MaxLength,Lane,BaseData.ILRes(g));
    
    % Round TrLineUp first row, move unrounded to fifth row
    TrLineUp(:,5) = TrLineUp(:,1); TrLineUp(:,1) = round(TrLineUp(:,1)/BaseData.ILRes(g));
    % TrLineUp [     1            2         3         4          5     ]
    %           AllTrAxIndex  AxleValue   Truck#    LaneID   Station(m)
    
    % For each influence case
    for t = 1:Num.InfCases
        
        % Reset for each t
        AllTrAxt = AllTrAx;
        TrLineUpt = TrLineUp;
        k = 0;
        
        % For each analysis
        while k < BaseData.NumAnalyses(g) && sum(AllTrAxt,'all') > 0
            
            % Subject Influence Line to Truck Axle Stream
            [MaxLE,DLF,BrStInd,R] = VBGetMaxLE(AllTrAxt,ILData(t).v,BaseData.RunDyn(g));
            
            % Get length of bridge in number of indices
            BrLengthInds = size(ILData(t).v,1);
            
            % Add Padding if necessary
            if BrStInd < 1 || BrStInd + BrLengthInds - 1 > height(AllTrAxt)
                % Add Padding
                PadLen = BrLengthInds -1;
                AllTrAxt = [zeros(PadLen,size(AllTrAxt,2)); AllTrAxt; zeros(PadLen,size(AllTrAxt,2))];
                BrStInd = BrStInd + PadLen;
                % Also need to modify TrLineUp
                TrLineUpt(:,1) = TrLineUpt(:,1) + PadLen; TrLineUpt(:,5) = TrLineUpt(:,5) + BaseData.ILRes(g)*PadLen;
            end
            
            BrEndInds = BrStInd + BrLengthInds-1;
            BrInds = [BrStInd:BrEndInds]';
            AxOnBr = sum(AllTrAxt(BrInds,:),2);
            
            % Now add to k since continue has passed
            k = k+1;

            % Optional Apercu
            if BaseData.Apercu(g) == 1
                ApercuTitle = Lane.Sites.SName + " " + num2str(BaseData.SITE(g)) + " Max " + num2str(k);
                T = VBApercuv2(PDs,ApercuTitle,ILData(t),BrStInd,TrLineUpt,DLF,Lane,BaseData.ILRes(g));
                %exportgraphics(gcf,"Apercu" + "/" + ApercuTitle + ".jpg",'Resolution',600)
            end
            
            % Prepare for next run - Set Axles to zero in AllTrAx (can't delete because indices are locations)
            AllTrAxt(BrInds,:) = 0;
            
        end
    end
    
    % Update progress bar
    UpProgBar(u, st, g, 1, height(BaseData), 1)
    
end

