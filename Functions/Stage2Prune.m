function [PDs] = Stage2Prune(PDs)
%STAGE2PRUNE

PDs(PDs.GW_TOT<6000,:) = [];
% Only do the disqualification if we actually have SW10 Classification
try
    if sum(PDs.CS == 2 | PDs.CS == 3 | PDs.CS == 4 | PDs.CS == 6) < 0.7*height(PDs)
        PDs(PDs.CS == 2 | PDs.CS == 3 | PDs.CS == 4 | PDs.CS == 6,:) = [];
    end
catch
end

PDs(PDs.SPEED > 130,:) = [];
PDs(PDs.SPEED < 10,:) = [];

end

