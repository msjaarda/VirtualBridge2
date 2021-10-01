function [Num,Lane,ILData,TrData,FolDist,ESIA] = VBUpdateData(BaseData)
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
    [A, B, Lane.Dir] = unique(Lane.Details.NSEW);
    for i = 1:length(Lane.Dir)
    if Lane.Details.NSEW(i) == 2 | Lane.Details.NSEW(i) == 4
        Lane.Dir(i) = 1;
    else
        Lane.Dir(i) = 2;
    end
    end
else % Lets try to form Lane.Sites and Lane.Details from Lane.Dir for the VBApercu
    % Get Lane Truck Distribution, Lane.TrDistr, and Lane Directions, Lane.Dir
    % If optional, do try
    
    % Splitting by ',' is no problem, even where there is no ','
    Lane.Dir =  cellfun(@str2num,split(BaseData.LaneDir{:},','));
    Lane.Sites = table();
    Lane.Sites.NumLanes = length(Lane.Dir);
    Lane.Sites.CANTON = "";
    Lane.Sites.HWY = "";
    Lane.Details = table();
    Lane.Details.LANE(1:Lane.Sites.NumLanes) = 1:Lane.Sites.NumLanes;
    Lane.Details.ALANE = Lane.Details.LANE;
    Lane.Details.NSEW(Lane.Dir == 1) = 3;
    Lane.Details.NSEW(Lane.Dir == 2) = 4;
    Lane.Details.FROM = repmat("",Lane.Sites.NumLanes,1);
    Lane.Details.DIR = repmat("",Lane.Sites.NumLanes,1);
end

try Lane.TrDistr =  cellfun(@str2num,split(BaseData.LaneTrDistr{:},',')); catch end
% Get Num.Lanes from the length of Lane.Dir
Num.Lanes = length(Lane.Dir);

% FolDist can use qualitative measures:
% "Jammed" or "Stopped" : 0 kph
% "At-rest" or "Crawling" : 2 kph
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
    [~, k] = max(ILData(i).v);
    
    clear b, clear c

    % Interpolate around influence lines to figure out next biggest max
    for j = 1:size(ILData(i).v,2)
        x = 0:BaseData.ILRes:(length(ILData(i).v)-1)*BaseData.ILRes;
        yILv = ILData(i).v(:,j);
        
        if j == 1
            PL = 300;
        elseif j == 2
            PL = 200;
        else
            PL = 0;
        end
        
        % Should be 1.3... changed on 30/9/21
        Conc = zeros(round(1.3/BaseData.ILRes),1);
        if max(size(Conc)) == 1
            Conc(1) = PL*2;
        else
            Conc(1) = PL; Conc(end) = PL;
        end
        
        MaxInfvCONV(j,i) = max(conv(Conc,yILv));
        
        % New convolution 25/7/21
        % Lets create a concentrated load matrix... then add padding
        % depending on length of IL..
        
        % Lucas 11/06/21 : We have to find the worst position of the axles,
        % it could be before or after k (Shear cases with a peak
        % value), or in between of k (+ or - Moment cases).
        
        % No need for previous if statement because interpolating outside
        % of x simply returns NaN
                    
            top(1) = interp1(x,ILData(i).v(:,j),x(k(j)));
            bot(1) = interp1(x,ILData(i).v(:,j),x(k(j))-1.2);
                        
            top(2) = interp1(x,ILData(i).v(:,j),x(k(j))+0.6);
            bot(2) = interp1(x,ILData(i).v(:,j),x(k(j))-0.6);
            
            top(3) = interp1(x,ILData(i).v(:,j),x(k(j))+1.2);
            bot(3) = interp1(x,ILData(i).v(:,j),x(k(j)));
            
            [~,Posmax] = max(top+bot);
            
            b(j) = top(Posmax); 
            c(j) = bot(Posmax);
    end
    
    aprime = (b+c)/2;
    MaxInfv(:,i) = aprime';

end

% Assign integral values into IntInfv (each InfCase)
for i = 1:Num.InfCases
    x = 0:BaseData.ILRes:(length(ILData(i).v)-1)*BaseData.ILRes;
    A = ILData(i).v;
    A(A<0) = 0;
    IntInfv(:,i) = trapz(x,A);
%  IntInfv(:,i) = trapz(ILData(i).v);
%     trapz(ILData(i).v);
%     BPlan = trapz(Infx(Infv(:,i)>=0),Infv(Infv(:,i)>=0,i));
end

% Define ESIA details
LaneWidth = 3; % meters, hard coded
% Initialize concentrated loads, Qk
Qk = zeros(Num.Lanes,1);
% Distributed loads
qk = 2.5*ones(Num.Lanes,1);
Qk(1) = 300; qk(1) = 9; % kN, kN/m2
% If there is more than 1 lane, the second lanes has 200 kN loads
if Num.Lanes > 1
    Qk(2) = 200;
end
% Alpha is 1 to make ratios easier (note that it is 0.9 in the code)
Alpha = 1;

% On 25.03.2021 Matt and Lucas used LucasInfluenceLine to show that this
% method underpredicts ESIA for twin girder bridges because in Lucas' code he
% shifts the point loads Q1 and Q2 to the edge, and I do not. TM did the
% same as Lucas.

% Calculate ESIA for each InfCase
for i = 1:Num.InfCases
    MaxvCONV = MaxInfvCONV(:,i);
    Maxv = MaxInfv(:,i);
    Intv = IntInfv(:,i);
%     ESIA.Total(i) = 1.5*Alpha*(Maxv'*Qk*2+Intv'*qk*LaneWidth);
%     ESIA.EQ(:,i) = Maxv.*Qk*2;
    ESIA.Total(i) = 1.5*Alpha*(sum(MaxvCONV)+Intv'*qk*LaneWidth);
    ESIA.EQ(:,i) = MaxvCONV;
    ESIA.Eq(:,i) = Intv.*qk*LaneWidth;
end

% if AnalysisType == WIM & Chan.
%    % Delete interior  
% end

% This method is not good enough when doing influence lines of only 2m and
% such... we need a better method. I think Convolution will solve it.
% Nothing wrong with the integration side of things, point loads are the
% problem... do convolution with padding... no multiplying by 2...

end

