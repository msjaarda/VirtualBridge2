% ------------------------------------------------------------------------
%                            VBWIMq
% ------------------------------------------------------------------------
% Generate a summary of maximum effects on bridges from real WIM
% This has a sister live script which loads the var and perform analyses.

% AxleStatsBasic >> AxTandem      >> Q1   Investigation
% VBWIMQ1Q2      >> Q1Q2MaxEvents >> Q1Q2 Investigation
% VBWIMq         >> qMaxEvents    >> q    Investigation

% Initializing commands
clear, clc, tic, format long g, rng('shuffle'), close all;

% Read Input File
BaseData = VBReadInputFile('Input/VBWIMqInputx.xlsx');

% Initialize parpool if necessary and initialize progress bar
if BaseData.Parallel(1) > 0, gcp; clc; end

% Initialize variables and start row counter
%MaxEvents = []; RamUsed = []; LenPrint = []; MaxEvents1 = [];

% Each row of BaseData represents one analysis OR analysis 'set'/'group'
for g = 1:height(BaseData)
    
    % Initialize variables and start row counter
    MaxEvents = []; RamUsed = []; LenPrint = []; %MaxEvents1 = [];

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
        
        % Load File
        load(['WIM/',num2str(BaseData.SITE(g)),'.mat']);
        
        % Stage2Prune
        PDs = Stage2Prune(PDs);
        
        % Separate for each year...
        if ismember('Year',BaseData.Properties.VariableNames)
            UYears = BaseData.Year(g);
        else
            UYears = unique(year(PDs.DTS));
        end
        
        % Start Progress Bar
        u = StartProgBar(length(UYears), 1, g, height(BaseData)); tic; st = now;
        
        %parfor r = 1:length(UYears)
        for r = 1:length(UYears)
            
            %MaxEvents1 = [];
            
            PDsy = PDs(year(PDs.DTS) == UYears(r),:);
            
            % Convert PDC to AllTrAx - Spacesave at MaxLength
            MaxLength = (max(arrayfun(@(x) size(x.v,1),ILData))-1)*BaseData.ILRes(g);
            [PDsy, AllTrAx, TrLineUp] = VBWIMtoAllTrAx(PDsy,MaxLength,Lane,BaseData.ILRes(g));
            
            % Make groups out of each unique day
            PDsy.Group = findgroups(dateshift(PDsy.DTS,'start','day'));
            
            % Round TrLineUp first row, move unrounded to fifth row
            TrLineUp(:,5) = TrLineUp(:,1); TrLineUp(:,1) = round(TrLineUp(:,1)/BaseData.ILRes(g));
            % Expand TrLineUp to include groups
            TrLineUp(:,6) = PDsy.Group(TrLineUp(:,3));
            % TrLineUp [ 1: AllTrAxIndex  2: AxleValue  3: Truck#  4: LaneID  5: Station(m)  6: Group  ]
            
            % In order to prevent broadcast variables, and instead have sliced
            % variables, particularly for TrLineUp, AllTrAx, and PDsy -- Necessary?
            for z = 1:max(PDsy.Group)
                TrLineUpGr{z} = TrLineUp(TrLineUp(:,6) == z,:);
                StartiGr{z} = max(1,min(TrLineUpGr{z}(:,1)));
                EndiGr{z} = min(max(TrLineUpGr{z}(:,1)),length(AllTrAx));
                AllTrAxGr{z} = AllTrAx(StartiGr{z}:EndiGr{z},:);
                %PDsyDTGr{z} = PDsy.DTS(TrLineUpGr{z}(1):TrLineUpGr{z}(end));
                %PDsyCLGr{z} = PDsy.CLASS(TrLineUpGr{z}(1):TrLineUpGr{z}(end));
            end
            
            % Perform search for maximums for each day
            parfor z = 1:max(PDsy.Group)
            %for z = 1:max(PDsy.Group)
            
                MaxEvents1 = [];
                
                % Get TrLineUpSub and AllTrAxSub
                TrLineUpSub = TrLineUpGr{z};
                Starti = StartiGr{z};
                Endi = EndiGr{z};
                AllTrAxSub = AllTrAxGr{z};
                %PDsyDTSub = PDsyDTGr{z};
                %PDsyCLSub = PDsyCLGr{z};
                
                % Don't bother running if the segment is too small
                if length(AllTrAxSub) < 20000/BaseData.ILRes(g), continue, end
                
                % For each InfCase
                %parfor t = 1:Num.InfCases
                for t = 1:Num.InfCases
                    
                    % Get length of bridge in number of indices
                    BrLengthInd = size(ILData(t).v,1);
                    
                    % Reset for each t
                    %AllTrAxSub = AllTrAx(Starti:Endi,:);
                    AllTrAxSub = AllTrAxGr{z};
                    
                    % Eliminate the need for padding or BrStInd index issues
                    AllTrAxSub(1:BrLengthInd,:) = 0; AllTrAxSub(end-BrLengthInd:end,:) = 0;
                    
                    % For each analysis
                    k = 0; % Initialize k
                    while k < BaseData.NumAnalyses(g) && sum(AllTrAxSub,'all') > 0
                        
                        % Subject Influence Line to Truck Axle Stream
                        [MaxLE,DLF,BrStInd,R] = VBGetMaxLE(AllTrAxSub,ILData(t).v,BaseData.RunDyn(g));
                        if MaxLE == 0, k = k+1, continue, end
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
                        TrNumsU = unique(TrNums);
                        
                        % We use TrNums because they don't depend on Starti shift
                        MaxLETime = PDsy.DTS(TrNums(1));
                        Vehs = PDsy.CLASS(TrNumsU);
                        
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
                        %MaxEvents1 = [MaxEvents1; datenum(MaxLETime), BaseData.SITE(g), MaxLE, t, m, k, L1Veh, L2Veh, L1Load, L2Load, L1Ax, L2Ax, L1Spd, L2Spd];
                        
                        % Prepare for next run - Set Axles to zero in AllTrAx (can't delete because indices are locations)
                        AllTrAxSub(BrInds-(Starti-1),:) = 0;
                        
                    end % k, analyses
                end % t, InfCases
                
                MaxEvents = [MaxEvents; MaxEvents1]; %MaxEvents1 = [];
            end % z, groups
            
            % Update progress bar
            user = memory;
            RamUsed = [RamUsed;user.MemUsedMATLAB/user.MemAvailableAllArrays*100];
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
    %MaxEvents = array2table(MaxEvents,'VariableNames',{'Datenum', 'SITE', 'MaxLE', 'InfCase', 'm', 'DayRank', 'L1Veh', 'L2Veh', 'L1Load', 'L2Load', 'L1Ax', 'L2Ax', 'L1Sp', 'L2Sp'});
    MaxEvents.DTS = datetime(MaxEvents.Datenum,'ConvertFrom','datenum');
    MaxEvents.Datenum = [];
    
    % Add Column for All, Class, ClassOW and delete former m
    MaxEvents.ClassT(MaxEvents.m == 1) = "All";
    MaxEvents.ClassT(MaxEvents.m == 2) = "ClassOW";
    MaxEvents.ClassT(MaxEvents.m == 3) = "Class";
    MaxEvents.m = [];
    
    % qInvestInitial Inputs
    BM = {'Daily', 'Weekly', 'Yearly'};             % j
    ClassType = {'All', 'ClassOW', 'Class'};        % i
    ClassT = {'All', 'Classified+', 'Classified'};
    DistTypes = {'Lognormal'};
    [Max,pd,x_values,y_values] = qInvestInitial(BM,ClassType,ClassT,DistTypes,MaxEvents,ESIA,ILData,BaseData(g,:));
    
    TName = datestr(now,'mmmdd-yy HHMMSS');
    OutInfo.Name = TName; OutInfo.BaseData = BaseData(g,:);
    OutInfo.ESIA = ESIA; %OutInfo.ESIM = ESIM;
    %OutInfo.OverMax = OverMax; OutInfo.OverMaxT = OverMaxT;
    OutInfo.ILData = ILData;
    
    % Plot BlockMax, find Design Values, Ed, using Beta, rather than 99th percentile
    for r = 1:Num.InfCases
        for i = 1:length(ClassType)
            Class = ClassType{i};
            BlockM = BM{2};
            %[pd,OutInfo.x_values(:,r),OutInfo.y_valuespdf(:,r),y_valuescdf] = GetBlockMaxFit(OutInfo.Max(:,r),'Lognormal',BaseData.Plots(g));
            %[ECDF,ECDFRank,PPx,PPy,Fity,OutInfo.LNFitR2] = GetLogNormPPP(Max(r).(Class).(BlockM).Max,false);
            [OutInfo.EdLN.(Class).(BlockM)(r), OutInfo.AQ.(Class).(BlockM)(r), ~] = GetBlockMaxEd(Max(r).(Class).(BlockM).Max,BlockM,'Lognormal',ESIA.Total(r),ESIA.EQ(:,r),ESIA.Eq(:,r),0.6,0.5);
        end
    end
    
    % Create folders where there are none
    CreateFolders(BaseData.Folder{g},BaseData.VWIM(g),BaseData.Apercu(g),BaseData.Save(g))
    
    if BaseData.Save(g) == 1
        save(['Output' BaseData.Folder{g} '/' OutInfo.Name], 'OutInfo')
    end
    
end % g, BaseData

% % Convert back to datetime !!
% MaxEvents = array2table(MaxEvents,'VariableNames',{'Datenum', 'SITE', 'MaxLE', 'InfCase', 'm', 'DayRank', 'BrStInd'});
% %MaxEvents = array2table(MaxEvents,'VariableNames',{'Datenum', 'SITE', 'MaxLE', 'InfCase', 'm', 'DayRank', 'L1Veh', 'L2Veh', 'L1Load', 'L2Load', 'L1Ax', 'L2Ax', 'L1Sp', 'L2Sp'});
% MaxEvents.DTS = datetime(MaxEvents.Datenum,'ConvertFrom','datenum');
% MaxEvents.Datenum = [];
% 
% % Add Column for All, Class, ClassOW and delete former m
% MaxEvents.ClassT(MaxEvents.m == 1) = "All";
% MaxEvents.ClassT(MaxEvents.m == 2) = "ClassOW";
% MaxEvents.ClassT(MaxEvents.m == 3) = "Class";
% MaxEvents.m = [];
% 
% % Call qInvestInitial (could opt to change name) (could remove conversion
% % to BlockMax strings since we unconvert in qInvestInitial...)
% 
% % Stuff below needs work...
% 
% % Can put into OverMax format (see OverMaxT below for titles)
% OverMaxT = array2table(OverMax,'VariableNames',{'InfCase','MaxLE','DLF','BrStInd','MaxDamage'});
% 
% % Manual Save for the big stuff, but OutInfo autosave to selected folder
% % Save structure variable with essential simulation information
% 
% % Below code from VBSim... must be modified
% 
% % Get simulation time
% Time = GetSimTime();
% 
% % In the future could add InfCaseName
% OverMaxT = array2table(OverMax,'VariableNames',{'InfCase','SimNum','BatchNum','MaxLE','DLF','BrStInd','MaxDamage'});
% 
% % Reshape OverMax results for output
% OverMax = sortrows(OverMax); OverMax = OverMax(:,4);
% OverMax = reshape(OverMax,BaseData.NumSims(g),Num.InfCases);
% 
% % Duplicate if only 1 Sim to avoid statistical errors
% if BaseData.NumSims(g) == 1
%     OverMax = [OverMax; OverMax];
% end
% 
% % Get ESIM and Ratio
% ESIM = 1.1*prctile(OverMax,99); Ratio = ESIM./ESIA.Total;
% 
% % Print Summary Stats to Command Window
% VBPrintSummary(BaseData(g,:),BatchSize,TrData,Num,VirtualWIM,Time,Lane.TrDistr)
% 
% % Create folders where there are none
% CreateFolders(BaseData.Folder{g},BaseData.VWIM(g),BaseData.Apercu(g),BaseData.Save(g))
% 
% TName = datestr(now,'mmmdd-yy HHMMSS');
% 
% % Convert VirtualWIM to Table and save if necessary
% if BaseData.VWIM(g) == 1
%     PD = array2table(VirtualWIM,'VariableNames',VWIMCols);
%     save(['VirtualWIM' BaseData.Folder{g} '/WIM_' TName], 'PD')
% end
% 
% % Covert Apercu to Table and save if necessary
% if BaseData.Apercu(g) == 1
%     PD = array2table(ApercuOverMax,'VariableNames',[VWIMCols 'InfCase']);
%     save(['Apercu' BaseData.Folder{g} '/AWIM_' TName], 'PD')
% end
% 
% 
%     
% OutInfo.Name = TName; OutInfo.BaseData = BaseData(g,:);
% OutInfo.ESIA = ESIA; OutInfo.ESIM = ESIM;
% OutInfo.OverMax = OverMax; OutInfo.OverMaxT = OverMaxT;
% OutInfo.ILData = ILData;
% 
% % Plot BlockMax, find Design Values, Ed, using Beta, rather than 99th percentile
% for r = 1:Num.InfCases
%     [pd,OutInfo.x_values(:,r),OutInfo.y_valuespdf(:,r),y_valuescdf] = GetBlockMaxFit(OutInfo.OverMax(:,r),'Lognormal',BaseData.Plots(g));
%     [ECDF,ECDFRank,PPx,PPy,Fity,OutInfo.LNFitR2] = GetLogNormPPP(OutInfo.OverMax(:,r),false);
%     [OutInfo.EdLN(r), OutInfo.AQ(r), ~] = GetBlockMaxEd(OverMax(:,r),BaseData.Period(g),'Lognormal',ESIA.Total(r),ESIA.EQ(:,r),ESIA.Eq(:,r),0.6,0.5);
% end
% % Apply Model Factor to EdLN and AQ
% OutInfo.AQ = 1.1*OutInfo.AQ;  OutInfo.EdLN = 1.1*OutInfo.EdLN;
% 
% if BaseData.Save(g) == 1
%     save(['Output' BaseData.Folder{g} '/' OutInfo.Name], 'OutInfo')
% end
