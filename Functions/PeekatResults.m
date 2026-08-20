function PeekatResults(InfLine,OutputFolder,BlockM,Class,DistTypes,NumAnalyses,Dyna)
% function to analyse individual results from AlphaSummaryPlot

%clear, clc, close all,
load('Sites.mat')
IncZ = 0;

%% INPUT
% Select the Infl Line to inspect
%InfLine = 'Box.Stand.Wid12.Bi.Simp.p0.Mn.S80';
%OutputFolder = 'WIMv17pr';
%BlockM = 'Weekly';
%Class = 'Class';
%DistTypes = 'All'; % fitting
%NumAnalyses = 4;
%Dyna = 0; % Run dynamique, 0) No 1) Yes

%% CODE
Dir_List = dir('Output');
File_List = {Dir_List.name}';

% check if WIM/SIM folder exist
if sum(strcmp(File_List,OutputFolder))>=1
else
   error(append('No ',OutputFolder,' Folder!!'));
end

Dir_List = dir(append('Output/',OutputFolder));
File_List = {Dir_List(:).name}';

%cleaning file list
File_List = File_List(~strcmp(File_List,'.')&~strcmp(File_List,'..')&contains(File_List,'.mat'));
File_List = erase(File_List,'.mat');
File_Listsim = {};
InflPosi = [];

for i=1:height(File_List)
    
    load(append('Output/',OutputFolder,'/',File_List{i},'.mat'));
    
    if OutInfo.BaseData.StopSim == 0
    ILDataNames = {OutInfo.ILData.Name}';
    if strcmp(OutInfo.BaseData.AnalysisType,'WIM')
        if sum(contains(ILDataNames,InfLine))>=1
            File_List = File_List{i};
            InflPosi = find(contains(ILDataNames,InfLine));
            break
        end
    elseif strcmp(OutInfo.BaseData.AnalysisType,'Sim')
        ClassSimCheck = strcmp(char(fieldnames(OutInfo.pd)),Class);
        if sum(contains(ILDataNames,InfLine))>=1 && ClassSimCheck
            File_Listsim{end+1} = File_List{i};
            InflPosi(end+1) = find(contains(ILDataNames,InfLine));
        end
    end
    end
end

if strcmp(OutInfo.BaseData.AnalysisType,'WIM') %Possible un jour checker si on peut avoir plusieurs résultats WIM... Normalement pas un probleme
Max = OutInfo.Max(InflPosi).(Class).(BlockM);
Data = Max.Max;
Max = sortrows(Max,3,'descend');
try EventDuration = OutInfo.BaseData.EventDuration; catch EventDuration = 1; end
try BETATarget = OutInfo.BaseData.Beta; catch BETATarget = 4.2; end
elseif strcmp(OutInfo.BaseData.AnalysisType,'Sim')
    Max = [];
    Data = [];
    for i = 1:length(File_Listsim)
    load(append('Output/',OutputFolder,'/',File_Listsim{i},'.mat'));
    Maxtmp = OutInfo.OverMaxT(OutInfo.OverMaxT.InfCase == InflPosi(i),:);
    Maxtmp.FileName = repmat(string(File_Listsim{i}), height(Maxtmp), 1);
    Maxtmp.SITE = repmat(string(OutInfo.BaseData.Traffic), height(Maxtmp), 1);
    Data(:,i) = Maxtmp.MaxLE;
    DataSite{i} = unique(Maxtmp.SITE);
    try EventDuration(i) = OutInfo.BaseData.EventDuration; catch EventDuration = 1; end
    try BETATarget(i) = OutInfo.BaseData.Beta; catch BETATarget = 4.2; end
    Maxtmp = sortrows(Maxtmp,4,'descend');
    Max = [Max;Maxtmp(1:NumAnalyses,:)];
    end
    Max = sortrows(Max,4,'descend');
    Max = Max(1:NumAnalyses,:);
    Max.DTS = datetime(str2double(regexp(Max.SITE, '\d{4}', 'match', 'once')), 1, 1);
else
    error('Ni SIM ni WIM reconnu');
end

if strcmp(OutInfo.BaseData.AnalysisType,'WIM')
    pd = GetFit(Data,BlockM,DistTypes,1,[IncZ BETATarget EventDuration]);
elseif strcmp(OutInfo.BaseData.AnalysisType,'Sim')
    for i = 1:width(Data)
    disp(DataSite{i});
    pd = GetFit(Data(:,i),BlockM,DistTypes,1,[IncZ BETATarget(i) EventDuration(i)]);
    end
end
%pd = GetFit(Data,BlockM,DistTypes,1,0);
%pd = GetFit(Data(~isoutlier(Data,'gesd')),BlockM,DistTypes,1,1);

% Find sites that generate Max cases
%[SitesAct,~,c] = unique(Max.SITE(1:NumAnalyses));
%for i=1:max(c)
%    NumSitesAct(i) = sum(c==i);
%end

if strcmp(OutInfo.BaseData.AnalysisType,'WIM')

%NumSitesAct = NumSitesAct';
SitesAct = Max.SITE(1:NumAnalyses);

% Create table
BaseData = table('Size',[height(SitesAct),11],'VariableTypes',["double","string","double","double","double","double","double","string","string","datetime","logical"]);
BaseData.Properties.VariableNames = {'SITE','ILs','ILRes','RunDyn','Parallel','Apercu','NumAnalyses','AnalysisType','ClassType','Date','SlowOnly'};
BaseData.SITE = SitesAct;
BaseData.ILs(:) = InfLine;
BaseData.ILRes(:) = OutInfo.BaseData.ILRes;
BaseData.Apercu(:) = 1;
BaseData.NumAnalyses(:) = 1;
BaseData.AnalysisType(:) = "WIM";
BaseData.ClassType(:) = Class;
BaseData.Date = Max.DTS(1:NumAnalyses);
BaseData.RunDyn(:) = Dyna;
BaseData.SlowOnly(:) = false;
RType = {'Uni';'Bi';'PUN'};
BaseData.RType(:) = convertCharsToStrings(RType{[contains(InfLine,'Uni');contains(InfLine,'Bi');contains(InfLine,'PUN')]});

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
    
     % Find and remove duplicates
     if Sites.Layout(Sites.SITE == BaseData.SITE(g)) == 11
            % Get Duplicates
            PDs = FindDup2(PDs,0,0);
            % Delete Duplicates - from L1
            PDs(PDs.Dup & PDs.LANE == 1,:) = [];
     end
    
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
            
            % Correction des positions pour retrouver le cas d'origine : Lucas 31.08.22 
            if MaxLE-Max.Max(g) ~=0
                TrLineUpt = CorrectPosApercu(MaxLE,BrStInd,BrLengthInds,AllTrAxt,TrLineUpt,Max.Max(g),ILData(t).v);
            end
            
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
                %ApercuTitle = Lane.Sites.SName + " " + num2str(BaseData.SITE(g)) + " Max " + num2str(k);
                ApercuTitle = Lane.Sites.SName + " " + num2str(BaseData.SITE(g)) + " Max " + num2str(g);
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

elseif strcmp(OutInfo.BaseData.AnalysisType,'Sim')

    load('ILLib.mat');
    [~,Lane,ILData,~,~] = VBUpdateData(OutInfo.BaseData);
    ILData = ILData([matches({ILData.Name},append('ILLib.',InfLine))]');

    Dir_List2 = dir('Apercu/');
    File_List2 = {Dir_List2(:).name}';

    %cleaning file list
    File_List2 = File_List2(~strcmp(File_List2,'.')&~strcmp(File_List2,'..')&~contains(File_List2,'.mat'));

    if sum(strcmp(File_List2,OutputFolder))>=1
        OutputFolder2 = OutputFolder;
    else
        % check if the file name can be similar to the OutputFolder name
        TF = cellfun(@(x) contains(OutputFolder, x), File_List2);
        if sum(TF)==1
            OutputFolder2 = string(File_List2(TF));
        elseif sum(TF)==0
            error(append('Apercu file not found! Name of the file : ',OutputFolder));
        else
            error('Too many files with similar names inside Apercu...');
        end
    end

    for i = 1:NumAnalyses
    load(append('Apercu/',OutputFolder2,'/','AWIM_',Max.FileName(i),'.mat'));
    Max.SITE(i) = Max.SITE(i) + " : Max " + num2str(i) + " / Simulation n°" + Max.SimNum(i);
    % Because we want to use the VBGetApercu code, we have to adapt the
    % given datas, PD is "pruned" to contain only the ploted data
    PD = PD(PD.SimNum==Max.SimNum(i),:);
    PD = PD(PD.InfCase==Max.InfCase(i),:);
    Max.InfCase(i) = 1; %set infl lane to 1 because we will only study InfLine
    PD.InfCase(:) = 1;
    PD.SimNum(:) = 1;
    [~,~,~] = VBGetApercu(PD,Max(i,:),1,ILData,Dyna,Lane,OutInfo.BaseData.ILRes);
    end

end

end

