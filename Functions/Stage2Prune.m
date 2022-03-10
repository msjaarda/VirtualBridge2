function [PD] = Stage2Prune(PD)
%STAGE2PRUNE

PD(PD.GW_TOT<6000,:) = [];
% Only do the disqualification if we actually have SW10 Classification
if sum(PD.CS == 2 | PD.CS == 3 | PD.CS == 4 | PD.CS == 6) < 0.7*height(PD)
    PD(PD.CS == 2 | PD.CS == 3 | PD.CS == 4 | PD.CS == 6,:) = [];
end
% 
          PDs(PDs.SPEED > 120,:) = [];
%         if BaseData.SITE(g) == 405 || BaseData.SITE(g) == 406
%             PDs(PDs.LANE == 1 & PDs.GAPT > 99.8,:) = [];
%         end
%         if sum(BaseData.SITE(g)== SiteGroups.Bi4L) == 0 && sum(BaseData.SITE(g)== SiteGroups.Uni3L) == 0
%             % Get Duplicates
%             PDs = FindDup2(PDs,0,0);
%             % Delete Duplicates - from L1
%             PDs(PDs.Dup & PDs.LANE == 1,:) = [];
%         end

end

