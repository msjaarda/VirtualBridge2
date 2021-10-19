function LaneTrDistr = CheckLaneTrDistr(LaneTrDistr)
%CheckLaneTrDistr This function chackes if LaneTrDistr is unique. If yes it
%returns only one set of vector, if not returns a 0 vector.

comp = strcmpi(LaneTrDistr,LaneTrDistr(1,1)) == 0;
if sum(sum(comp)) == 0
    LaneTrDistr = LaneTrDistr(1,1);
else
    LaneTrDistr = strsplit(char(LaneTrDistr(1,1)),',');
    [~,TrSize] = size(LaneTrDistr);
    LaneTrDistr = append(repmat('0,',1,TrSize-1),'0');
    LaneTrDistr = {LaneTrDistr};
end

end

