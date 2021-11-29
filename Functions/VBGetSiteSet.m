function [Sites] = VBGetSiteSet(SITE,StopSim)
%VBGetSiteSet Summary of this function goes here
%   Detailed explanation goes here


    load('SiteGroups.mat')
    if SITE == 11
        Sites = SiteGroups.('Uni2L');
    elseif SITE == 111
        Sites = SiteGroups.('Uni3L');
    elseif SITE == 12
        Sites = SiteGroups.('Bi2L');
        % Don't use Simplong for StopSim
        if StopSim
            Sites(Sites == 441) = [];
        end
    elseif SITE == 1122
        Sites = SiteGroups.('Bi4L');
    elseif SITE == 110
        Sites = SiteGroups.('LSVAUni2L');
    else
        Sites = SITE;
    end
    
    
end

