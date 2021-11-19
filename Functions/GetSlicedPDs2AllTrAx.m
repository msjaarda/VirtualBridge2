function [TrLineUpGr,PDsy] = GetSlicedPDs2AllTrAx(PDsy,MaxLength,Lane,ILRes)
%GETSLICEDPDS2ALLTRAX Effective way of only have Sliced vars in mem

% Convert PDsy to AllTrAx
[PDsy, TrLineUp] = VBWIMtoAllTrAx_Exp(PDsy,MaxLength,Lane,ILRes);

% Make groups out of each unique day
if Lane.Sites.SITE == 460
    PDsy.Group = findgroups(dateshift(PDsy.DTS,'start','week'));
else
    PDsy.Group = findgroups(dateshift(PDsy.DTS,'start','day'));
end

% Round TrLineUp first row, move unrounded to fifth row
TrLineUp(:,5) = TrLineUp(:,1); TrLineUp(:,1) = round(TrLineUp(:,1)/ILRes);
% Expand TrLineUp to include groups
TrLineUp(:,6) = PDsy.Group(TrLineUp(:,3));
% TrLineUp [ 1: AllTrAxIndex  2: AxleValue  3: Truck#  4: LaneID  5: Station(m)  6: Group  ]

% In order to prevent broadcast variables, and instead have sliced
% variables, particularly for TrLineUp, and AllTrAx
for z = 1:max(PDsy.Group)
    TrLineUpGr{z} = TrLineUp(TrLineUp(:,6) == z,:);
end
                   
            
end

