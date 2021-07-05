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
% Start just with 2 lane locations
BaseData = VBReadInputFile('Input/VBWIMqInput.xlsx');

% Initialize parpool if necessary and initialize progress bar
if BaseData.Parallel(1) > 0, gcp; clc; end

% Initialize variables and start row counter
MaxEvents = [];

% Each row of BaseData represents one analysis
for g = 1:height(BaseData)

    % Update analysis data for current row of BaseData
    [Num,Lane,ILData,~,~,ESIA] = VBUpdateData(BaseData(g,:));
    
    % Load File
    load(['WIM/',num2str(BaseData.SITE(g)),'.mat']);
    
    %Stage2Prune
    PDs(PDs.GW_TOT<6000,:) = [];
    % Only do the disqualification if we actually have SW10 Classification
    if sum(PDs.CS == 2 | PDs.CS == 3 | PDs.CS == 4 | PDs.CS == 6) < 0.7*height(PDs)
        PDs(PDs.CS == 2 | PDs.CS == 3 | PDs.CS == 4 | PDs.CS == 6,:) = [];
    end

    % Separate for each year...
    UYears = unique(year(PDs.DTS));
    
    % Start Progress Bar
    u = StartProgBar(length(UYears), 1, g, height(BaseData)); tic; st = now;
    
    %parfor r = 1:length(UYears)
    for r = 1:length(UYears)
        
        %MaxEvents1 = nan(500000,14);
        %j = 1;
        MaxEvents1 = [];
        
        PDsy = PDs(year(PDs.DTS) == UYears(r),:);
    
        % Convert PDC to AllTrAx - Spacesave at MaxLength
        MaxLength = (max(arrayfun(@(x) size(x.v,1),ILData))-1)*BaseData.ILRes(g);
        [PDsy, AllTrAx, TrLineUp] = VBWIMtoAllTrAx(PDsy,MaxLength,Lane,BaseData.ILRes(g));
        
        % Make groups out of each unique day
        % In the end make this weekly
        PDsy.Group = findgroups(dateshift(PDsy.DTS,'start','day'));
        
        % Round TrLineUp first row, move unrounded to fifth row
        TrLineUp(:,5) = TrLineUp(:,1); TrLineUp(:,1) = round(TrLineUp(:,1)/BaseData.ILRes(g));
        % Expand TrLineUp to include groups
        TrLineUp(:,6) = PDsy.Group(TrLineUp(:,3));
        
        % TrLineUp [     1            2         3         4          5     ]
        %           AllTrAxIndex  AxleValue   Truck#    LaneID   Station(m)
    
        % Perform search for maximums for each day
        for z = 1:max(PDsy.Group)
            
            % Store starting and end indices
            Starti = max(1,min(TrLineUp(TrLineUp(:,6) == z,1)));
            Endi = min(max(TrLineUp(TrLineUp(:,6) == z,1)),length(AllTrAx));
            
            % Subdivide AllTrAx
            AllTrAxSub = AllTrAx(Starti:Endi,:);
            TrLineUpGr = TrLineUp(TrLineUp(:,6) == z,:);
            
            % Don't bother running if the segment is too small
            if length(AllTrAxSub) < 20000
                continue
            end
            
            % For each InfCase
            %for t = 1:Num.InfCases
            parfor t = 1:Num.InfCases
                
                % Get length of bridge in number of indices
                BrLengthInd = size(ILData(t).v,1);                
                
                % Reset for each t
                AllTrAxSub = AllTrAx(Starti:Endi,:);
                k = 0;
                
                % Eliminate the need for padding or BrStInd index issues
                AllTrAxSub(1:BrLengthInd,:) = 0;
                AllTrAxSub(end-BrLengthInd:end,:) = 0;
            
                % For each analysis
                while k < BaseData.NumAnalyses(g) && sum(AllTrAxSub,'all') > 0
                    
                    % Subject Influence Line to Truck Axle Stream
                    %if BrLengthInd/BaseData.ILRes(g) < 60
                        [MaxLE,DLF,BrStInd,R] = VBGetMaxLE(AllTrAxSub,ILData(t).v,BaseData.RunDyn(g));
                    %else
                        %[MaxLE,DLF,BrStInd,R] = VBGetMaxLE(AllTrAxSub,ILData(t).v,0);
                    %end
                    
                    % Now add to k
                    k = k+1;
                    
                    % Adjust BrStInd for Starti
                    BrStIndx = BrStInd + Starti -1;
                    
                    % Get BrEndIndx
                    BrEndIndx = BrStIndx + BrLengthInd - 1;
                    
                    % Get Bridge Indices
                    BrInds = [BrStIndx:BrEndIndx]';
                    %AxOnBr = sum(AllTrAxt(StripInds,:),2);
                    
                    % Get Key Info to Save
                    TrNums = TrLineUpGr(TrLineUpGr(:,1) >= min(BrInds) & TrLineUpGr(:,1) <= max(BrInds),3);
                    TrNumsU = unique(TrNums);
                    
                    % We use TrNums because they don't depend on Starti shift
                    MaxLETime = PDsy.DTS(TrNums(1));
                    Vehs = PDsy.CLASS(TrNumsU);
                    Spds = PDsy.SPEED(TrNumsU);
                    Lnes = PDsy.LANE(TrNumsU);
                    L1Veh = numel(Vehs(Lnes == 1));
                    L2Veh = numel(Vehs(Lnes == 2));
                    L1Spd = mean(Spds(Lnes == 1));
                    L2Spd = mean(Spds(Lnes == 2));
                    % 99 is coded as empty (0 is taken by unclassified)
                    if isempty(L1Veh); L1Veh = 99; end
                    if isempty(L2Veh); L2Veh = 99; end
                    if isnan(L1Spd); L1Spd = -1; end
                    if isnan(L2Spd); L2Spd = -1; end
                    
                    % When retreiving from AllTrAxt we must shift back
                    L1Load = sum(AllTrAxSub(BrInds-(Starti-1),1));
                    L2Load = sum(AllTrAxSub(BrInds-(Starti-1),2));
                    L1Ax = sum(AllTrAxSub(BrInds-(Starti-1),1)>0);
                    L2Ax = sum(AllTrAxSub(BrInds-(Starti-1),2)>0);
                    
                    % Get ClassT (in m form for now)
                    if min(Vehs) == 0
                        m = 1;
                    elseif sum(Vehs > 39 & Vehs < 50) > 0
                        m = 2;
                    else
                        m = 3;
                    end
                    
                    % Optional Apercu
                    if BaseData.Apercu(g) == 1
                        ApercuTitle = Lane.Sites.SName + " " + num2str(BaseData.SITE(g)) + " Max " + num2str(k);
                        T = VBApercu(PDsy,ApercuTitle,ILData(t),BrStIndx,TrLineUp,MaxLE/ESIA.Total(t),DLF,Lane,BaseData.ILRes(g));
                        %exportgraphics(gcf,"Apercu" + BaseData.Folder(g) + "/" + ApercuTitle + ".jpg",'Resolution',600)
                    end
                    
                    % Save MaxEvents...
                    % Save Times and Datenums and then convert
                    %MaxEvents1(j,:) = [datenum(MaxLETime), BaseData.SITE(g), MaxLE, t, m, k, L1Veh, L2Veh, L1Load, L2Load, L1Ax, L2Ax, L1Spd, L2Spd];
                    %j = j+1;
                    MaxEvents1 = [MaxEvents1; datenum(MaxLETime), BaseData.SITE(g), MaxLE, t, m, k, L1Veh, L2Veh, L1Load, L2Load, L1Ax, L2Ax, L1Spd, L2Spd];
                    
                    % Prepare for next run - Set Axles to zero in AllTrAx (can't delete because indices are locations)
                    AllTrAxSub(BrInds-(Starti-1),:) = 0;
                    
                end % k, analyses
            end % t, InfCases
        end % z, groups
    
        MaxEvents= [MaxEvents; MaxEvents1];
    
        % Update progress bar
        UpProgBar(u, st, g, 1, length(UYears), 1)
        
    end % r, years
    
    
end % g, BaseData

% Trim back
%MaxEvents = reshape(permute(MaxEvents,[1 3 2]),[],14);
%MaxEvents(isnan(MaxEvents(:,1)),:) = [];

% Convert back to datetime !!
% Delete empty rows and convert to table
MaxEvents = array2table(MaxEvents,'VariableNames',{'Datenum', 'SITE', 'MaxLE', 'InfCase', 'm', 'DayRank', 'L1Veh', 'L2Veh', 'L1Load', 'L2Load', 'L1Ax', 'L2Ax', 'L1Sp', 'L2Sp'});
MaxEvents.DTS = datetime(MaxEvents.Datenum,'ConvertFrom','datenum');
MaxEvents.Datenum = [];

% Add Column for All, Class, ClassOW and delete former m
MaxEvents.ClassT(MaxEvents.m == 1) = "All";
MaxEvents.ClassT(MaxEvents.m == 2) = "ClassOW";
MaxEvents.ClassT(MaxEvents.m == 3) = "Class";
MaxEvents.m = [];

% Manual Save