function [PDCx, AllTrAx, TrLineUp] = VBWIMtoAllTrAx(PDCx,SpaceSaver,Lane,ILRes)
% WIMTOALLTRAX Translates WIM or VWIM data into AllTrAx and TrLineUp
% Also returns PDC (in the form of PDCx) with some mods

[PDCx, TrLineUp] = VBWIM2TrLineUp(PDCx,SpaceSaver,Lane);

% Get Lanes
Lanes = unique(PDCx.LANE);

% The way that the indexing and accumarray is working, we have wasted stuff
% at the start of the AllTrAx... and it is much too long (when using VWIM)
TrLineUp(:,1) = round(TrLineUp(:,1)/ILRes);

% Make a separate axle stream vector for each lane, and last one for all
% Put max() function in incase one lane has no representation in TrLineUp
AllTrAx = zeros(max(TrLineUp(:,1)),max(length(Lane.Details.Dir),length(Lanes)));

for i = 1:length(Lanes)
    A = accumarray(TrLineUp(TrLineUp(:,4)==Lanes(i),1),TrLineUp(TrLineUp(:,4)==Lanes(i),2));
    AllTrAx(1:length(A),i) = A;
end

% Return TrLineUp first row unrounded
TrLineUp(:,1) = WBv;

end




