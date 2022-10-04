function [Num,Lane,ILData,TrData,FolDist] = VBUpdateData(BaseData)
%Updates the data for the main variables, based on a signle row in BaseData

% --- Required in BaseData --- Depends on AnalysisType
% LaneTrDistr
% Flow
% Traffic

if strcmp(BaseData.AnalysisType,"WIM")
    
    % We need Lane.Dir and Num.Lanes for other things to work!
    load('Sites.mat'); load('SiteLanes.mat'); 
    
    % Get LaneDir automatically using Sites and SiteLanes
    Lane.Sites = Sites(Sites.SITE == BaseData.SITE,:);
    Lane.Details = SiteLanes(SiteLanes.SITE == BaseData.SITE,:);
    Num.Lanes = Lane.Sites.NumLanes;

    % Assign Lane.Details.ALANE here
    % Start by doing all with NSEW == 2 or 4
    IDs = Lane.Details.NSEW == 2 | Lane.Details.NSEW == 4;
    if sum(IDs) > 0
        Check = sortrows(Lane.Details(IDs,:),{'LANE'});
        % If the first one is a Lane 1 we good
        if Check.LANETYPE(1) == 1
            Check.ALANE = [1:height(Check)]';
        else
            Check = sortrows(Lane.Details(IDs,:),{'LANE'},'descend');
            Check.ALANE = [1:height(Check)]';
        end
        % Now check for other IDs
        IDs2 = Lane.Details.NSEW == 1 | Lane.Details.NSEW == 3;
        if sum(IDs2) > 0
            Check2 = sortrows(Lane.Details(IDs2,:),{'LANE'});
            % If the first one is a Lane 1 we good
            if Check2.LANETYPE(1) == 1
                Check2.ALANE = [height(Check2)+height(Check):-1:height(Check)+1]';
            else
                Check2 = sortrows(Lane.Details(IDs2,:),{'LANE'},'descend');
                Check2.ALANE = [height(Check2)+height(Check):-1:height(Check)+1]';
            end
        else
            Check2 = [];
        end
    else
        Check = [];
        Check2 = sortrows(Lane.Details,{'LANE'});
        % If the first one is a Lane 1 we good
        if Check2.LANETYPE(1) == 1
            Check2.ALANE = [height(Check2)+height(Check):-1:height(Check)+1]';
        else
            Check2 = sortrows(Lane.Details(IDs2,:),{'LANE'},'descend');
            Check2.ALANE = [height(Check2)+height(Check):-1:height(Check)+1]';
        end
    end
    Lane.Details = [Check; Check2];
    
    %[~, ~, Lane.Dir] = unique(Lane.Details.NSEW);
    for i = 1:height(Lane.Details)
        if Lane.Details.NSEW(i) == 2 || Lane.Details.NSEW(i) == 4
            Lane.Details.Dir(i) = 1;
        else
            Lane.Details.Dir(i) = 2;
        end
    end
      
else % Lets try to form Lane.Sites and Lane.Details from Lane.Dir for the VBApercu
    % Get Lane Truck Distribution, Lane.TrDistr, and Lane Directions, Lane.Dir
    % If optional, do try
    
    % Splitting by ',' is no problem, even where there is no ','
    Lane.Details = table();
    Lane.Details.Dir =  cellfun(@str2num,split(BaseData.LaneDir{:},','));
    Lane.Sites = table();
    Lane.Sites.NumLanes = length(Lane.Details.Dir);
    Lane.Sites.STATE = "";
    Lane.Sites.HWY = "";
    %Lane.Details = table();
    Lane.Details.LANE(1:Lane.Sites.NumLanes) = 1:Lane.Sites.NumLanes;
    Lane.Details.ALANE = Lane.Details.LANE;
    Lane.Details.NSEW(Lane.Details.Dir == 1) = 3;
    Lane.Details.NSEW(Lane.Details.Dir == 2) = 4;
    Lane.Details.FROM = repmat("",Lane.Sites.NumLanes,1);
    Lane.Details.DIR = repmat("",Lane.Sites.NumLanes,1);
    
end

try Lane.TrDistr =  cellfun(@str2num,split(BaseData.LaneTrDistr{:},',')); catch end
% Get Num.Lanes from the length of Lane.Dir
Num.Lanes = length(Lane.Details.Dir);

% FolDist can use qualitative measures:
% "Jammed" or "Stopped" : 0 kph   -- Koshini Stopped
% "At-rest" or "Crawling" : 2 kph -- Koshini w/ velocity as input (light and heavy treated the same)
% "Congested" : 30 kph
% "Free-flowing" : 1000 veh/hr

% Alternatively, one can simply use a single number representing speed in kph (for
% 100 and under), OR volume in veh/hr (any value over 100 is assume as volume)

% Update FolDist
if strcmp(BaseData.AnalysisType,"Sim")
    FolDist = array2table(zeros(4,4));
    % Note that we include truck and car transitions, even if not jammed (simpler coding)
    FolDist.Properties.VariableNames = {'TaT', 'TaC', 'CaT', 'CaC'}; % "a" means after TaC is Truck after Car <<<car<<<<TRUCK
    if iscell(BaseData.Flow)
        if strcmp(BaseData.Flow{:},'Jammed') || strcmp(BaseData.Flow{:},'Stopped')
            FolDist.TaT = [0.1 15 2.93 10.8]';
            FolDist.TaC = [0.1 15 2.15 10.9]';
            FolDist.CaT = [0.1 15 2.41 9.18]';
            FolDist.CaC = [0.1 15 2.15 15.5]';
            VehSpd = 0; % kph
        elseif strcmp(BaseData.Flow{:},'At-rest') || strcmp(BaseData.Flow{:},'Crawling')
            VehSpd = 2; % kph
        elseif strcmp(BaseData.Flow{:},'Congested')
            VehSpd = 30; % kph
        elseif strcmp(BaseData.Flow{:},'Free-flowing')
            VehSpd = 1000; % > 100 therefore, veh/hr
        else
            fprintf('\nWarning: Not a recognized FolDist input\n')
        end
        %  end
    else
        VehSpd = BaseData.Flow;
    end
    if VehSpd > 0 && VehSpd < 101
        FolDist.TaT = [VehSpd/15 15+1.1*VehSpd 2.15 9]';
        FolDist.TaC = FolDist.TaT; FolDist.CaT = FolDist.TaT; FolDist.CaC = FolDist.TaT;
    elseif VehSpd > 100
        % Difficult task... what to do with flowing traffic. Has a large
        % effect inside PerLaneRates AND GetFloDist... we opt to represent
        % exponential distribution AS beta! See Free-movingFollowing.xlsx
        beta = VehSpd*0.0126-0.2490; alpha = 1;
        mind = 5.5; maxd = 1000; % m
        FolDist.TaT = [mind maxd alpha beta]';
        FolDist.TaC = FolDist.TaT; FolDist.CaT = FolDist.TaT; FolDist.CaC = FolDist.TaT;
    elseif VehSpd == 0
        FolDist.TaT = [0.1 15 2.93 10.8]';
        FolDist.TaC = [0.1 15 2.15 10.9]';
        FolDist.CaT = [0.1 15 2.41 9.18]';
        FolDist.CaC = [0.1 15 2.15 15.5]';
    end
else
    FolDist = [];
end

% Update TrData
% Load TrLib if necessary
if strcmp(BaseData.AnalysisType,"Sim")
    if ~exist('TrLib','var')
        load('TrLib.mat')
    end
    % Ensure that the chosen traffic exists in TrLib
    if isfield(TrLib,BaseData.Traffic{:})
        % Overwrite
        TrData = TrLib.(BaseData.Traffic{:});
    else
        fprintf('\nWarning: Traffic input not recognized\n\n')
    end
else
    TrData = [];
end

% Update LaneData
% Load ILLib if necessary
if ~exist('ILLib','var')
    load('ILLib.mat')
end
% GetInfLines function is done next - we prepare LaneData for that function

% This is key. We want seamless input here
% Box
% Twin
% Multi
% Slab

% - Must work with IL which have different x steps
% - Num lanes is from traffic input. Use multidimentional arrays. When only 1
% dim is given, this is equivalent to a "0", applying to all lanes
% - In the future we can make "Area average" for computation of AGB area
% loads from the code
% - You can specify groups or individual lines. Dot notation for groups. Take
% everything downstream of the dot!

% BaseData.ILs = ILFamily.SpecificILName,OtherILFamily,OtherILFamily. ...

% A generic AGBBox or AGBTwin.Standard.Mn means all in the family
% Example: AGBBox.Mp.S80 means just Span of 80 for M+ Box Girder
% See ILGuide.xls

% Split input by the commas to get individual IL families
ILs = split(BaseData.ILs{:},',');
[Num.InfCases, ILData] = findIL(ILs,BaseData.ILRes,Num.Lanes);

% NOTE - sometimes the "track average" (average of two wheel positions) is
% not equal to the "area average", and so ESIA calculations which involve
% the placement of area loads will be wrong. We can add extra ILs in these
% locations which correspond to the area average for the purpose of ESIA
% calculation... we can also add custom lines for the purpose of fixing the
% issue mentioned at the bottom (Twin Girder error based on truck
% placement) even though I fundamentally disagree with TM there.

% Flip signs, if necessary...
for i = 1:Num.InfCases 
    % NOTE: We only consider the first lane when decided if we should flip
    % Keep in mind... may not always be true
    
    % Switch signs of all ILs associated with InfCase together if warranted
    if abs(max(ILData(i).v(:,1))) < abs(min(ILData(i).v(:,1)))
        ILData(i).v = -ILData(i).v;
    end
    
    % Find max Influence line value index, k
    %[~, k] = max(ILData(i).v);
    
    clear b, clear c
     
 end

end

