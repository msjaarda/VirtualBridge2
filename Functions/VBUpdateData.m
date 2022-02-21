function [Num,Lane,ILData,TrData,FolDist,E] = VBUpdateData(BaseData)
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
    %[~, k] = max(ILData(i).v);
    
    clear b, clear c

    % Interpolate around influence lines to figure out next biggest max
    for j = 1:size(ILData(i).v,2)
        %x = 0:BaseData.ILRes:(length(ILData(i).v)-1)*BaseData.ILRes;
        yILv = ILData(i).v(:,j);
        
        if j == 1
            PL = 300;
            PLBR(1) = 210; % Bruhwiler values from 15 october 2019 document
            PLBR(2) = 170;
            PLBR(3) = 170;
            PL12 = 120; % 12 tonnes per axle (load model for 41, 48 and 72)
        elseif j == 2
            PL = 200;
            PLBR(1) = 0;
            PLBR(2) = 130;
            PLBR(3) = 125;
        else
            PL = 0;
            PLBR(1) = 0;
            PLBR(2) = 0;
            PLBR(3) = 0;
        end
        
        % Should be 1.3... changed on 30/9/21
        % Lucas 23/12/21 : FALSE, should be 1.2 + add 1 to the size
        
        % Lucas 23/12/21 : Not correct, for resolution of 1 we must have 2
        % values for conc and not only multiplied by 2. For higher
        % resolution than 2.4 we can admit only one value multiplied by 2.
        Conc = zeros(round(1.2/BaseData.ILRes)+1,1);
        ConcBR1 = PLBR(1);
        ConcBR2 = zeros(round(1.2/BaseData.ILRes)+1,1);
        Conc41 = zeros(round(7.2/BaseData.ILRes)+1,1);
        Conc48 = zeros(round(9.1/BaseData.ILRes)+1,1);
        Conc72 = zeros(round(14/BaseData.ILRes)+1,1);
        
        if max(size(Conc)) == 1
            Conc(1) = PL*2;
            ConcBR2(1) = PLBR(2)*2;
        else
            Conc(1) = PL; Conc(end) = PL;
            ConcBR2(1) = PLBR(2); ConcBR2(end) = PLBR(2);
        end
        
        ConcBR3 = zeros(round(2.4/BaseData.ILRes)+1,1);
        if max(size(ConcBR3)) == 1
            ConcBR3(1) = PLBR(3)*3;
        elseif max(size(ConcBR3)) == 2
            ConcBR3(1) = PLBR(3)*1.5; ConcBR3(end) = PLBR(3)*1.5;
        else
            ConcBR3(1) = PLBR(3); ConcBR3(end) = PLBR(3);
            try
            ConcBR3((end+1)/2) = PLBR(3);
            catch
            ConcBR3((end)/2) = PLBR(3)*.5; ConcBR3((end)/2+1) = PLBR(3)*.5;
            end
        end
        
        if j == 1 %Load model for the 41 crane, du to resolution must add previous load if applied in the same spot
           Conc41(1) = PL12*0;
           Conc41(round(2.4/BaseData.ILRes)+1) = Conc41(round(2.4/BaseData.ILRes)+1) + PL12*0;
           Conc41(round(4/BaseData.ILRes)+1) = Conc41(round(4/BaseData.ILRes)+1) + PL12*0;
           Conc41(round(5.6/BaseData.ILRes)+1) = Conc41(round(5.6/BaseData.ILRes)+1) + PL12*0;
           Conc41(end) = Conc41(end) + PL12*0;
        else
           Conc41 = Conc;
        end
        
        if j == 1 %Load model for the 48 crane, du to resolution must add previous load if applied in the same spot
           Conc48(1) = PL12;
           Conc48(round(2.6/BaseData.ILRes)+1) = Conc48(round(2.6/BaseData.ILRes)+1) + PL12;
           Conc48(round(4.2/BaseData.ILRes)+1) = Conc48(round(4.2/BaseData.ILRes)+1) + PL12;
           Conc48(round(5.9/BaseData.ILRes)+1) = Conc48(round(5.9/BaseData.ILRes)+1) + PL12;
           Conc48(round(7.5/BaseData.ILRes)+1) = Conc48(round(7.5/BaseData.ILRes)+1) + PL12;
           Conc48(end) = Conc48(end) + PL12;
        else
           Conc48 = Conc;
        end
        
        if j == 1 %Load model for the 48 crane, du to resolution must add previous load if applied in the same spot
           Conc72(1) = PL12;
           Conc72(round(3.2/BaseData.ILRes)+1) = Conc72(round(3.2/BaseData.ILRes)+1) + PL12;
           Conc72(round(4.6/BaseData.ILRes)+1) = Conc72(round(4.6/BaseData.ILRes)+1) + PL12;
           Conc72(round(9.8/BaseData.ILRes)+1) = Conc72(round(9.8/BaseData.ILRes)+1) + PL12;
           Conc72(round(11.2/BaseData.ILRes)+1) = Conc72(round(11.2/BaseData.ILRes)+1) + PL12;
           Conc72(round(12.6/BaseData.ILRes)+1) = Conc72(round(12.6/BaseData.ILRes)+1) + PL12;
           Conc72(end) = Conc72(end) + PL12;
        else
           Conc72 = Conc;
        end
        
        MaxInfvCONV(j,i) = max(conv(Conc,yILv));
        [MaxInfvCONVBR1(j,i),MaxInfvCONVBR1Posi(j,i)] = max(conv(ConcBR1,yILv));
        [MaxInfvCONVBR2(j,i),MaxInfvCONVBR2Posi(j,i)] = max(conv(ConcBR2,yILv));
        [MaxInfvCONVBR3(j,i),MaxInfvCONVBR3Posi(j,i)] = max(conv(ConcBR3,yILv));
        [MaxInfvCONV41(j,i),MaxInfvCONV41Posi(j,i)] = max(conv(Conc41,yILv));
        [MaxInfvCONV48(j,i),MaxInfvCONV48Posi(j,i)] = max(conv(Conc48,yILv));
        [MaxInfvCONV72(j,i),MaxInfvCONV72Posi(j,i)] = max(conv(Conc72,yILv));
        
        % New convolution 25/7/21
        % Lets create a concentrated load matrix... then add padding
        % depending on length of IL..
        
        % Lucas 11/06/21 : We have to find the worst position of the axles,
        % it could be before or after k (Shear cases with a peak
        % value), or in between of k (+ or - Moment cases).
        
        % No need for previous if statement because interpolating outside
        % of x simply returns NaN
            
            % Old Lucas
            %top(1) = interp1(x,ILData(i).v(:,j),x(k(j)));
            %bot(1) = interp1(x,ILData(i).v(:,j),x(k(j))-1.2);
                        
            %top(2) = interp1(x,ILData(i).v(:,j),x(k(j))+0.6);
            %bot(2) = interp1(x,ILData(i).v(:,j),x(k(j))-0.6);
            
            %top(3) = interp1(x,ILData(i).v(:,j),x(k(j))+1.2);
            %bot(3) = interp1(x,ILData(i).v(:,j),x(k(j)));
            
            %[~,Posmax] = max(top+bot);
            
            %b(j) = top(Posmax); 
            %c(j) = bot(Posmax);
    end
    
    % Old Lucas
    % aprime = (b+c)/2;
    % MaxInfv(:,i) = aprime';

end

% Assign integral values into IntInfv (each InfCase)
for i = 1:Num.InfCases
    x = 0:BaseData.ILRes:(length(ILData(i).v)-1)*BaseData.ILRes;
    A = ILData(i).v;
    A(A<0) = 0;
    IntInfv(:,i) = trapz(x,A);
    for j = 1:width(A)
    IntInfvBR1(j,i) = IntInfv(j,i)-(j==1)*trapz(x(max(MaxInfvCONVBR1Posi(j,i)-round(2/BaseData.ILRes),1):min(MaxInfvCONVBR1Posi(j,i)+round(2/BaseData.ILRes),end)),A(max(MaxInfvCONVBR1Posi(j,i)-round(2/BaseData.ILRes),1):min(MaxInfvCONVBR1Posi(j,i)+round(2/BaseData.ILRes),end),j));
    IntInfvBR2(j,i) = IntInfv(j,i)-(j==1||j==2)*trapz(x(max(MaxInfvCONVBR2Posi(j,i)-round(3.2/BaseData.ILRes),1):min(MaxInfvCONVBR2Posi(j,i)+round(2/BaseData.ILRes),end)),A(max(MaxInfvCONVBR2Posi(j,i)-round(3.2/BaseData.ILRes),1):min(MaxInfvCONVBR2Posi(j,i)+round(2/BaseData.ILRes),end),j));
    IntInfvBR3(j,i) = IntInfv(j,i)-(j==1||j==2)*trapz(x(max(MaxInfvCONVBR3Posi(j,i)-round(4.4/BaseData.ILRes),1):min(MaxInfvCONVBR3Posi(j,i)+round(2/BaseData.ILRes),end)),A(max(MaxInfvCONVBR3Posi(j,i)-round(4.4/BaseData.ILRes),1):min(MaxInfvCONVBR3Posi(j,i)+round(2/BaseData.ILRes),end),j));
    IntInfv41(j,i) = IntInfv(j,i)-(j==1)*trapz(x(max(MaxInfvCONV41Posi(j,i)-round(8.9/BaseData.ILRes),1):min(MaxInfvCONV41Posi(j,i)+round(4.3/BaseData.ILRes),end)),A(max(MaxInfvCONV41Posi(j,i)-round(8.9/BaseData.ILRes),1):min(MaxInfvCONV41Posi(j,i)+round(4.3/BaseData.ILRes),end),j));
    IntInfv48(j,i) = IntInfv(j,i)-(j==1)*trapz(x(max(MaxInfvCONV48Posi(j,i)-round(11/BaseData.ILRes),1):min(MaxInfvCONV48Posi(j,i)+round(4.9/BaseData.ILRes),end)),A(max(MaxInfvCONV48Posi(j,i)-round(11/BaseData.ILRes),1):min(MaxInfvCONV48Posi(j,i)+round(4.9/BaseData.ILRes),end),j));
    IntInfv72(j,i) = IntInfv(j,i)-(j==1)*trapz(x(max(MaxInfvCONV72Posi(j,i)-round(15/BaseData.ILRes),1):min(MaxInfvCONV72Posi(j,i)+round(1.5/BaseData.ILRes),end)),A(max(MaxInfvCONV72Posi(j,i)-round(15/BaseData.ILRes),1):min(MaxInfvCONV72Posi(j,i)+round(1.5/BaseData.ILRes),end),j));
    end   
%  IntInfv(:,i) = trapz(ILData(i).v);
%     trapz(ILData(i).v);
%     BPlan = trapz(Infx(Infv(:,i)>=0),Infv(Infv(:,i)>=0,i));
end

% Define ESIA details
LaneWidth = 3; % meters, hard coded
% Initialize concentrated loads, Qk
% Qk = zeros(Num.Lanes,1);
% Distributed loads
qk = 2.5*ones(Num.Lanes,1);
qkBR = qk;
qk05 = qk; % kN/m2 with alpha = 1 (needed alpha = 0.5)
%Qk(1) = 300; 
qk(1) = 9; % kN/m2
qkBR(1) = 3.6; % kN/m2
qk05(1) = 9; % kN/m2 with alpha = 1 (needed alpha = 0.5)
% If there is more than 1 lane, the second lanes has 200 kN loads
%if Num.Lanes > 1
%    Qk(2) = 200;
%end
% Alpha is 1 to make ratios easier (note that it is 0.9 in the code)
Alpha = 1;

% On 25.03.2021 Matt and Lucas used LucasInfluenceLine to show that this
% method underpredicts ESIA for twin girder bridges because in Lucas' code he
% shifts the point loads Q1 and Q2 to the edge, and I do not. TM did the
% same as Lucas.

% Calculate ESIA for each InfCase
for i = 1:Num.InfCases
    MaxvCONV = MaxInfvCONV(:,i);
    MaxvCONVBR1 = MaxInfvCONVBR1(:,i);
    MaxvCONVBR2 = MaxInfvCONVBR2(:,i);
    MaxvCONVBR3 = MaxInfvCONVBR3(:,i);
    MaxvCONV41 = MaxInfvCONV41(:,i);
    MaxvCONV48 = MaxInfvCONV48(:,i);
    MaxvCONV72 = MaxInfvCONV72(:,i);
    % Maxv = MaxInfv(:,i); Old Lucas
    Intv = IntInfv(:,i);
    IntvBR1 = IntInfvBR1(:,i);
    IntvBR2 = IntInfvBR2(:,i);
    IntvBR3 = IntInfvBR3(:,i);
    Intv41 = IntInfv41(:,i);
    Intv48 = IntInfv48(:,i);
    Intv72 = IntInfv72(:,i);
%   E.Total(i) = 1.5*Alpha*(Maxv'*Qk*2+Intv'*qk*LaneWidth);
%   E.EQ(:,i) = Maxv.*Qk*2;
    E.Total(i) = 1.5*Alpha*(sum(MaxvCONV)+Intv'*qk*LaneWidth);
    E.EQ(:,i) = MaxvCONV;
    E.Eq(:,i) = Intv.*qk*LaneWidth;
    [E.EBRU.Total(i),posi] = max([1.5*(sum(MaxvCONVBR1)+IntvBR1'*qkBR*LaneWidth);1.5*(sum(MaxvCONVBR2)+IntvBR2'*qkBR*LaneWidth);1.5*(sum(MaxvCONVBR3)+IntvBR3'*qkBR*LaneWidth)]);
    E.EBRU.LoadMod{i} = append('Model n°',int2str(posi));
    if posi == 1
    E.EBRU.EQ(:,i) = MaxvCONVBR1;
    E.EBRU.Eq(:,i) = IntvBR1.*qkBR*LaneWidth;
    elseif posi == 2
    E.EBRU.EQ(:,i) = MaxvCONVBR2;
    E.EBRU.Eq(:,i) = IntvBR2.*qkBR*LaneWidth;
    elseif posi == 3
    E.EBRU.EQ(:,i) = MaxvCONVBR3;
    E.EBRU.Eq(:,i) = IntvBR3.*qkBR*LaneWidth;
    else
        error('error with EBRU');
    end
    E.E41.Total(i)= 1.5*Alpha*(sum(MaxvCONV41)+Intv41'*qk05*LaneWidth);
    E.E41.EQ(:,i) = MaxvCONV41;
    E.E41.Eq(:,i) = Intv41.*qk05*LaneWidth;
    E.E48.Total(i)= 1.5*Alpha*(sum(MaxvCONV48)+Intv48'*qk05*LaneWidth);
    E.E48.EQ(:,i) = MaxvCONV48;
    E.E48.Eq(:,i) = Intv48.*qk05*LaneWidth;
    E.E72.Total(i)= 1.5*Alpha*(sum(MaxvCONV72)+Intv72'*qk05*LaneWidth);
    E.E72.EQ(:,i) = MaxvCONV72;
    E.E72.Eq(:,i) = Intv72.*qk05*LaneWidth;
    
end

% if AnalysisType == WIM & Chan.
%    % Delete interior  
% end

% This method is not good enough when doing influence lines of only 2m and
% such... we need a better method. I think Convolution will solve it.
% Nothing wrong with the integration side of things, point loads are the
% problem... do convolution with padding... no multiplying by 2...

end

