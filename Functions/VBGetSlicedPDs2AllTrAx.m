function [TrLineUpGr,PDsy] = VBGetSlicedPDs2AllTrAx(PDsy,MaxLength,Lane,ILRes,SliceUnit)
%VBGETSLICEDPDS2ALLTRAX Effective way of only have Sliced vars in mem
% SliceUnit is normally 'week' (fatigue) or 'day' ULS

% Convert PDsy to AllTrAx
[PDsy, TrLineUp] = VBWIM2TrLineUp(PDsy,MaxLength,Lane);

% Make groups out of each unique SliceUnit
PDsy.Group = findgroups(dateshift(PDsy.DTS,'start',SliceUnit));

% Round TrLineUp first row, move unrounded to fifth row
TrLineUp(:,5) = TrLineUp(:,1); TrLineUp(:,1) = round(TrLineUp(:,1)/ILRes);
% Expand TrLineUp to include groups
TrLineUp(:,6) = PDsy.Group(TrLineUp(:,3));
% TrLineUp [ 1: AllTrAxIndex  2: AxleValue  3: Truck#  4: LaneID  5: Station(m)  6: Group  ]
TrLineUp = array2table(TrLineUp,'VariableNames',{'ATAIndex','AxleValue','TrNum','LaneID','mStation','Group'});

% In order to prevent broadcast variables, and instead have sliced
% variables, particularly for TrLineUp, and AllTrAx
for z = 1:max(PDsy.Group)
    TrLineUpGr{z} = TrLineUp(TrLineUp.Group == z,:);
end         
            
end

