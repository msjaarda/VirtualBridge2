function [PD] = Stage2Prune(PD)
%STAGE2PRUNE

PD(PD.GW_TOT<6000,:) = [];
% Only do the disqualification if we actually have SW10 Classification
if sum(PD.CS == 2 | PD.CS == 3 | PD.CS == 4 | PD.CS == 6) < 0.7*height(PD)
    PD(PD.CS == 2 | PD.CS == 3 | PD.CS == 4 | PD.CS == 6,:) = [];
end

PD(PD.SPEED > 120,:) = [];


end

