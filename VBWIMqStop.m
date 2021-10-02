% ------------------------------------------------------------------------
%                            VBWIMq
% ------------------------------------------------------------------------
% Generate a summary of maximum effects on bridges from real WIM
% This has a sister live script which loads the var and perform analyses.

% AxleStatsBasic >> AxTandem      >> Q1   Investigation
% VBWIMQ1Q2      >> Q1Q2MaxEvents >> Q1Q2 Investigation
% VBWIMq         >> qMaxEvents    >> q    Investigation

% Must re-write a little for the sake of memory... split PDs up into years
% perhaps.

% Initializing commands
clear, clc, tic, format long g, rng('shuffle'), close all;

% Read Input File
BaseData = VBReadInputFile('Input/VBWIMqInputx.xlsx');

% Initialize parpool if necessary and initialize progress bar
if BaseData.Parallel(1) > 0, gcp; clc; end

% Each row of BaseData represents one analysis OR analysis 'set'/'group'
for g = 1:height(BaseData)
    
    % Initialize variables and start row counter
    MaxEvents = []; RamUsed = []; LenPrint = [];
    MaxEventsStop = [];
    
    % Recognize if BaseData.SITE(g) is actually a 'set'
    load('SiteGroups.mat')
    if BaseData.SITE(g) == 11, Sites = SiteGroups.('Uni2L');
    elseif BaseData.SITE(g) == 111, Sites = SiteGroups.('Uni3L');
    elseif BaseData.SITE(g) == 12, Sites = SiteGroups.('Bi2L');
    elseif BaseData.SITE(g) == 1122, Sites = SiteGroups.('Bi4L');
    else Sites = BaseData.SITE(g);
    end
        
    for w = 1:length(Sites)
        
        % Modify BaseData.SITE(g) based on SiteSet
        BaseData.SITE(g) = Sites(w);
        % Update analysis data for current row of BaseData
        [Num,Lane,ILData,~,~,ESIA] = VBUpdateData(BaseData(g,:));
        
        % Get MaxLength for Spacesave
        MaxLength = (max(arrayfun(@(x) size(x.v,1),ILData))-1)*BaseData.ILRes(g);
        
        % Load File
        if BaseData.LSVA(g)
            load(['WIMLSVA/',num2str(BaseData.SITE(g)),'.mat']);
        else
            load(['WIM/',num2str(BaseData.SITE(g)),'.mat']);
            % Stage2Prune (not for LSVA)
            PDs = Stage2Prune(PDs);
        end
        
        % Separate for each year...
        if ismember('Year',BaseData.Properties.VariableNames)
            UYears = BaseData.Year(g);
        else
            UYears = unique(year(PDs.DTS));
        end
        
        % Start Progress Bar
         u = StartProgBar(length(UYears), 1, g, height(BaseData)); tic; st = now;
        
        for r = 1:length(UYears)
            
            % Try to save mem by clearing PDs and re-loading
            [PDsy] = LoadPDYear(['WIMLSVA/',num2str(BaseData.SITE(g)),'.mat'],UYears(r));
            
            % Get TrLineUp, AllTrAx, Starti and Endi in sliced form
            [TrLineUpGr,StartiGr,~,AllTrAxGr,PDsy] = GetSlicedPDs2AllTrAx(PDsy,MaxLength,Lane,BaseData.ILRes(g));
            
%             % Convert PDsy to AllTrAx
%             [PDsy, AllTrAx, TrLineUp] = VBWIMtoAllTrAx(PDsy,MaxLength,Lane,BaseData.ILRes(g));
%             
%             % Make groups out of each unique day
%             PDsy.Group = findgroups(dateshift(PDsy.DTS,'start','day'));
%             
%             % Round TrLineUp first row, move unrounded to fifth row
%             TrLineUp(:,5) = TrLineUp(:,1); TrLineUp(:,1) = round(TrLineUp(:,1)/BaseData.ILRes(g));
%             % Expand TrLineUp to include groups
%             TrLineUp(:,6) = PDsy.Group(TrLineUp(:,3));
%             % TrLineUp [ 1: AllTrAxIndex  2: AxleValue  3: Truck#  4: LaneID  5: Station(m)  6: Group  ]
%             
%             % In order to prevent broadcast variables, and instead have sliced
%             % variables, particularly for TrLineUp, and AllTrAx
%             for z = 1:max(PDsy.Group)
%                 TrLineUpGr{z} = TrLineUp(TrLineUp(:,6) == z,:);
%                 StartiGr{z} = max(1,min(TrLineUpGr{z}(:,1)));
%                 EndiGr{z} = min(max(TrLineUpGr{z}(:,1)),length(AllTrAx));
%                 AllTrAxGr{z} = AllTrAx(StartiGr{z}:EndiGr{z},:);
%             end
            
            % Perform search for maximums for each day
            %parfor z = 1:max(PDsy.Group)
            for z = 1:max(PDsy.Group)
                
                % Initialize
                MaxEvents1 = []; MaxEvents1Stop = [];
                
                % Get TrLineUpSub and AllTrAxSub
                TrLineUpSub = TrLineUpGr{z};
                Starti = StartiGr{z};
                AllTrAxSub = AllTrAxGr{z};
                
                % Don't bother running if the segment is too small
                if length(AllTrAxSub) < 20000/BaseData.ILRes(g), continue, end
                
                % For each InfCase
                for t = 1:Num.InfCases
                    
                    % Get length of bridge in number of indices
                    BrLengthInd = size(ILData(t).v,1);
                    
                    % Reset for each t
                    AllTrAxSub = AllTrAxGr{z};
                    
                    % Eliminate the need for padding or BrStInd index issues
                    AllTrAxSub(1:BrLengthInd,:) = 0; AllTrAxSub(end-BrLengthInd:end,:) = 0;
                    
                    % For each analysis
                    k = 0; % Initialize k
                    while k < BaseData.NumAnalyses(g) && sum(AllTrAxSub,'all') > 0
                        
                        % Subject Influence Line to Truck Axle Stream
                        [MaxLE,DLF,BrStInd,R] = VBGetMaxLE(AllTrAxSub,ILData(t).v,BaseData.RunDyn(g));
                        if MaxLE == 0, k = k+1; continue, end
                        k = k+1; % Add to k
                        
                        % Adjust BrStInd for Starti
                        BrStIndx = BrStInd + Starti -1;
                        % Get BrEndIndx
                        BrEndIndx = BrStIndx + BrLengthInd - 1;
                        % Get Bridge Indices
                        BrInds = [BrStIndx:BrEndIndx]';
                        %AxOnBr = sum(AllTrAxt(StripInds,:),2);
                        
                        % Get Key Info to Save
                        TrNums = TrLineUpSub(TrLineUpSub(:,1) >= min(BrInds) & TrLineUpSub(:,1) <= max(BrInds),3);
                        TrLineUpOnBr = TrLineUpSub(TrLineUpSub(:,1) >= min(BrInds) & TrLineUpSub(:,1) <= max(BrInds),:);
                        [MaxM, MaxI] = max(TrLineUpOnBr(:,2));
                        TrIdMax = TrLineUpOnBr(MaxI,3);
                        
                        TrNumsU = unique(TrNums);
                        
                        NumExtra = 20;
                        TrNumsUE = [max(1,TrNumsU(1)-NumExtra):min(TrNumsU(end)+NumExtra,height(PDsy))]';
                        
                        PDe = PDsy(TrNumsUE,:);
                    
                        % Call VBWIMtoAllTrAx w/ mods... must give stationary point or truck
                        [PDe, AllTrAxStop, TrLineUpStop] = VBWIMtoAllTrAxStop(PDe,MaxLength,Lane,BaseData.ILRes(g),find(TrNumsUE == TrIdMax));
                        
                        % Round TrLineUp first row, move unrounded to fifth row
                        TrLineUpStop(:,5) = TrLineUpStop(:,1); TrLineUpStop(:,1) = round(TrLineUpStop(:,1)/BaseData.ILRes(g));
                        % TrLineUpStop [ 1: AllTrAxIndex  2: AxleValue  3: Truck#  4: LaneID  5: Station(m) ]
                        
                        [MaxLEe,DLFe,BrStInde,Re] = VBGetMaxLE(AllTrAxStop,ILData(t).v,BaseData.RunDyn(g));
  
                        % We use TrNums because they don't depend on Starti shift
                        MaxLETime = PDsy.DTS(TrNums(1));
                        Vehs = PDsy.CLASS(TrNumsU);
                        
                        %T = VBApercu(PDsy,'',ILData(t),BrStIndx,TrLineUp,MaxLE/ESIA.Total(t),1,Lane,BaseData.ILRes(g));
                        %TStop = VBApercu(PDe,'',ILData(t),BrStInde,TrLineUpStop,MaxLEe/ESIA.Total(t),1,Lane,BaseData.ILRes(g));
                        
                        % Only collect detailed info if desired
                        %[L1Veh,L2Veh,L1Spd,L2Spd,L1Load,L2Load,L1Ax,L2Ax] = DetailedVBWIM(PDsy,TrNumsU,Vehs,AllTrAxSub,BrInds,Starti);
                        
                        % Get ClassT (in m form for now)
                        if min(Vehs) == 0
                            m = 1;
                        elseif sum(Vehs > 39 & Vehs < 50) > 0
                            m = 2;
                        else
                            m = 3;
                        end
                        
                        % Save MaxEvents... save Times and Datenums and then convert
                        MaxEvents1 = [MaxEvents1; datenum(MaxLETime), BaseData.SITE(g), MaxLE, t, m, k, BrStIndx];
                        MaxEvents1Stop = [MaxEvents1Stop; datenum(MaxLETime), BaseData.SITE(g), MaxLEe, t, m, k, BrStInde];
                        
                        if m == 3
                            % Bump k up so that analysis doesn't continue!
                            k = 100;
                        end
                        
                        % Prepare for next run - Set Axles to zero in AllTrAx (can't delete because indices are locations)
                        AllTrAxSub(BrInds-(Starti-1),:) = 0;
                        
                    end % k, analyses
                end % t, InfCases
                
                MaxEvents = [MaxEvents; MaxEvents1];
                MaxEventsStop = [MaxEventsStop; MaxEvents1Stop];
            end % z, groups
            
            % Update progress bar
            user = memory;
            RamUsed = [RamUsed;user.MemUsedMATLAB/(user.MemAvailableAllArrays+user.MemUsedMATLAB)*100];
            LenPrint = VBUpProgBar(u,st,g,1,length(UYears),1,RamUsed(end),r,LenPrint);
            
        end % r, years
    end % w, SiteGroups
    
    % Optional Apercu... just for the first InfCase if there are multiple... works best if BaseData.Year is given
    if BaseData.Apercu(g) == 1
        SortedME = sortrows(MaxEvents,3);
        for c = 1:BaseData.NumAnalyses(g)
            ApercuTitle = Lane.Sites.SName + " " + num2str(BaseData.SITE(g)) + " " + num2str(BaseData.Year(g)) + " Max " + num2str(c);
            T = VBApercu(PDsy,ApercuTitle,ILData(1),SortedME(c,7),TrLineUp,SortedME(c,3)/ESIA.Total(1),1,Lane,BaseData.ILRes(g));
            %T = VBApercu(PDsy,ApercuTitle,ILData(t),BrStIndx,TrLineUp,MaxLE/ESIA.Total(t),DLF,Lane,BaseData.ILRes(g));
            %exportgraphics(gcf,"Apercu" + BaseData.Folder(g) + "/" + ApercuTitle + ".jpg",'Resolution',600)
        end
    end
  
    % Convert back to datetime
    MaxEvents = array2table(MaxEvents,'VariableNames',{'Datenum', 'SITE', 'MaxLE', 'InfCase', 'm', 'DayRank', 'BrStInd'});
    MaxEventsStop = array2table(MaxEventsStop,'VariableNames',{'Datenum', 'SITE', 'MaxLE', 'InfCase', 'm', 'DayRank', 'BrStInd'});
    %MaxEvents = array2table(MaxEvents,'VariableNames',{'Datenum', 'SITE', 'MaxLE', 'InfCase', 'm', 'DayRank', 'L1Veh', 'L2Veh', 'L1Load', 'L2Load', 'L1Ax', 'L2Ax', 'L1Sp', 'L2Sp'});
    MaxEvents.DTS = datetime(MaxEvents.Datenum,'ConvertFrom','datenum');
    MaxEventsStop.DTS = datetime(MaxEventsStop.Datenum,'ConvertFrom','datenum');
    MaxEvents.Datenum = [];
    MaxEventsStop.Datenum = [];
    
    % Add Column for All, Class, ClassOW and delete former m
    % CONSIDER DOING ALL THIS AFTER qINVESTINITIAL... delete first part of
    % qINVESTINITIAL...
    MaxEvents.ClassT(MaxEvents.m == 1) = "All";
    MaxEvents.ClassT(MaxEvents.m == 2) = "ClassOW";
    MaxEvents.ClassT(MaxEvents.m == 3) = "Class";
    MaxEvents.m = [];
    MaxEventsStop.ClassT(MaxEventsStop.m == 1) = "All";
    MaxEventsStop.ClassT(MaxEventsStop.m == 2) = "ClassOW";
    MaxEventsStop.ClassT(MaxEventsStop.m == 3) = "Class";
    MaxEventsStop.m = [];
    
    % qInvestInitial Inputs
    %BM = {'Daily', 'Weekly', 'Yearly'};             % j
    BM = {'Daily'}; 
    ClassType = {'All', 'ClassOW', 'Class'};        % i
    DistTypes = {'Lognormal'};
    [Max,~,~,~] = qInvestInitial(BM,ClassType,DistTypes,MaxEvents,ILData);
    
    TName = datestr(now,'mmmdd-yy HHMMSS');
    OutInfo.Name = TName; OutInfo.BaseData = BaseData(g,:);
    OutInfo.ESIA = ESIA; %OutInfo.ESIM = ESIM;
    %OutInfo.OverMax = OverMax; OutInfo.OverMaxT = OverMaxT;
    OutInfo.ILData = ILData;
    OutInfo.SimStop = false;
    
    % Plot BlockMax, find Design Values, Ed, using Beta, rather than 99th percentile
    for r = 1:Num.InfCases
        for i = 1:length(ClassType)
            Class = ClassType{i};
            %BlockM = BM{2};
            BlockM = BM{1};
            [~,OutInfo.x_values.(Class).(BlockM)(:,r),OutInfo.y_valuespdf.(Class).(BlockM)(:,r),~] = GetBlockMaxFit(Max(r).(Class).(BlockM).Max,'Lognormal',BaseData.Plots(g));
            %[ECDF,ECDFRank,PPx,PPy,Fity,OutInfo.LNFitR2] = GetLogNormPPP(Max(r).(Class).(BlockM).Max,false);
            [OutInfo.EdLN.(Class).(BlockM)(r), OutInfo.AQ.(Class).(BlockM)(r), ~] = GetBlockMaxEd(Max(r).(Class).(BlockM).Max,BlockM,'Lognormal',ESIA.Total(r),ESIA.EQ(:,r),ESIA.Eq(:,r),0.6,0.5);
        end
    end
    
    % Create folders where there are none
    CreateFolders(BaseData.Folder{g},BaseData.VWIM(g),BaseData.Apercu(g),BaseData.Save(g))
    
    if BaseData.Save(g) == 1
        save(['Output' BaseData.Folder{g} '/' OutInfo.Name], 'OutInfo')
    end
    
    MaxEventsStop(MaxEventsStop.MaxLE <= 0,:) = [];
    [Max,~,~,~] = qInvestInitial(BM,ClassType,DistTypes,MaxEventsStop,ILData);
    
    TName = datestr(now+1/86400,'mmmdd-yy HHMMSS');
    OutInfo.Name = TName;
    OutInfo.SimStop = true;
    
    % Plot BlockMax, find Design Values, Ed, using Beta, rather than 99th percentile
    for r = 1:Num.InfCases
        for i = 1:length(ClassType)
            Class = ClassType{i};
            %BlockM = BM{2};
            BlockM = BM{1};
            [~,OutInfo.x_values.(Class).(BlockM)(:,r),OutInfo.y_valuespdf.(Class).(BlockM)(:,r),~] = GetBlockMaxFit(Max(r).(Class).(BlockM).Max,'Lognormal',BaseData.Plots(g));
            %[ECDF,ECDFRank,PPx,PPy,Fity,OutInfo.LNFitR2] = GetLogNormPPP(Max(r).(Class).(BlockM).Max,false);
            [OutInfo.EdLN.(Class).(BlockM)(r), OutInfo.AQ.(Class).(BlockM)(r), ~] = GetBlockMaxEd(Max(r).(Class).(BlockM).Max,BlockM,'Lognormal',ESIA.Total(r),ESIA.EQ(:,r),ESIA.Eq(:,r),0.6,0.5);
        end
    end
    
    if BaseData.Save(g) == 1
        save(['Output' BaseData.Folder{g} '/' OutInfo.Name], 'OutInfo')
    end
    
end % g, BaseData

%T = VBApercu(PDsy,'',ILData(t),BrStIndx,TrLineUp,MaxLE/ESIA.Total(t),1,Lane,BaseData.ILRes(g));
%T = VBApercu(PDe,'',ILData(t),BrStInde,TrLineUpStop,MaxLEe/ESIA.Total(t),1,Lane,BaseData.ILRes(g));
