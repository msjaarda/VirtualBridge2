%TrLineUpt is a function that reorganize the trucks axles position in order to
%match the Max value found during simulation

function TrLineUpt = CorrectPosApercu(MaxLE,BrStInd,BrLengthInds,AllTrAxt,TrLineUpt,Max,ILData)

BestMaxLE = [MaxLE,0];
AllTrAxtempo = AllTrAxt(BrStInd:(BrStInd + BrLengthInds-1),:);
PosiAllTr = find(AllTrAxtempo~=0);
LoadAllTr = AllTrAxtempo(PosiAllTr);
PosiAllTr(:,2) = PosiAllTr(:,1);PosiAllTr(:,1)=PosiAllTr(:,2)-1;PosiAllTr(:,3)=PosiAllTr(:,2)+1; % Warning : by adding or removing 1, an axle can change lane, repair done bellow
PosiAllTr(mod(PosiAllTr(:,1),BrLengthInds)==0,1) = PosiAllTr(mod(PosiAllTr(:,1),BrLengthInds)==0,2); % axle gone on lower lane or outside
PosiAllTr(mod(PosiAllTr(:,3),BrLengthInds)==1,3) = PosiAllTr(mod(PosiAllTr(:,3),BrLengthInds)==1,2); % axle gone on upper lane or outside
if height(PosiAllTr)<=16 % For the moment, impossible to play with too many axles
NameFunction = 'PosiAllCases = combvec(';
for j=1:height(PosiAllTr)
    NameFunction = append(NameFunction,'PosiAllTr(',int2str(j),',:),');
end
NameFunction(end) = [];
NameFunction = append(NameFunction,')'';');
eval(NameFunction);
for j=1:height(PosiAllCases)
    if width(unique(PosiAllCases(j,:)))<width(PosiAllCases(j,:))
        % If two or more axles are on the same position, skip
    else
        AllTrAxtempo = zeros(BrLengthInds,width(AllTrAxtempo));
        AllTrAxtempo(PosiAllCases(j,:)) = LoadAllTr;
        MaxLEtemp = 0;
        for k=1:width(AllTrAxtempo)
            %MaxLEtemp = MaxLEtemp + max(AllTrAxtempo(:,k)'*ILData(:,k),AllTrAxtempo(:,k)'*flip(ILData(:,k)));
            MaxLEtemp = MaxLEtemp + AllTrAxtempo(:,k)'*ILData(:,k);
        end
        if MaxLEtemp-Max ==0
            BestMaxLE = [MaxLEtemp,j];
            break
        end
        if abs(MaxLEtemp-Max)<=abs(BestMaxLE(1)-Max)
            BestMaxLE = [MaxLEtemp,j];
        end
        MaxLEtemp = 0; %obliged to run it 2 times due to possible fliped infl
        for k=1:width(AllTrAxtempo)
            MaxLEtemp = MaxLEtemp + AllTrAxtempo(:,k)'*flip(ILData(:,k));
        end
        if MaxLEtemp-Max ==0
            BestMaxLE = [MaxLEtemp,j];
            break
        end
        if abs(MaxLEtemp-Max)<=abs(BestMaxLE(1)-Max)
            BestMaxLE = [MaxLEtemp,j];
        end
    end
end
AllTrAxtempo = zeros(BrLengthInds,width(AllTrAxtempo));
AllTrAxtempo(PosiAllCases(BestMaxLE(2),:)) = LoadAllTr;
OriginalPosi = PosiAllTr(:,2);
%OriginalPosi = PosiAllTr(:,2) + BrStInd-1;
NewPosi = PosiAllCases(BestMaxLE(2),:)';
%NewPosi = PosiAllCases(BestMaxLE(2),:)' + BrStInd-1;
for j=1:width(AllTrAxtempo)
    PosiLineUpt = ismember(TrLineUpt(:,1),OriginalPosi(OriginalPosi>(j-1)*(BrLengthInds) & OriginalPosi<=(j)*(BrLengthInds))-(j-1)*BrLengthInds+ BrStInd-1).*TrLineUpt(:,4)==j;
    try CheckDir = TrLineUpt(PosiLineUpt); CheckDir = CheckDir(end)-CheckDir(1); catch CheckDir=0; end
    if CheckDir >=0 %if CheckDir>=0 good direction, otherwhise need to flip to go in the other direction
    TrLineUpt(PosiLineUpt) = NewPosi(NewPosi>(j-1)*(BrLengthInds) & NewPosi<=(j)*(BrLengthInds))-(j-1)*BrLengthInds+ BrStInd-1;
    else
    TrLineUpt(PosiLineUpt) = flip(NewPosi(NewPosi>(j-1)*(BrLengthInds) & NewPosi<=(j)*(BrLengthInds))-(j-1)*BrLengthInds+ BrStInd-1);
    end
end
end
end