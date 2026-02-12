function [PDs] = Stage2Prune(PDs)
%STAGE2PRUNE

PDs(PDs.GW_TOT<6000,:) = [];
%this disqualification of 6t (instead of 3.5 t) is made to have better consistency with SW10
%classification.
try
    if sum(PDs.CS == 2 | PDs.CS == 3 | PDs.CS == 4 | PDs.CS == 6) < 0.7*height(PDs)
        PDs(PDs.CS == 2 | PDs.CS == 3 | PDs.CS == 4 | PDs.CS == 6,:) = [];
    end
catch
end

PDs(PDs.SPEED > 120,:) = [];
PDs(PDs.SPEED < 15,:) = [];

end

