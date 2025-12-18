function [Sitesx] = VBGetSiteSet(SITE,LightVehs,Core,Country)
%VBGetSiteSet

% We have eliminated SiteGroups.mat (placed it in Legacy)
% Therefore, now we used simply Sites.mat

% Country shall be alpha2 code (CH, DE, AT, US)

% FOR LightVehs
% 0 is NOT INCLUDED, 1 IS INCLUDED, 2 IS EXCLUSIVELY INCLUDED

% Can give Country or SITE == 0... if so, it does not filter

% Light gives you the option to ONLY (2), include (1), or disclude (not given or 0) light

load('Sites.mat')
SitesM = Sites;
if all(ismember(SITE,Sites.SITE))
    Sitesx = SITE;
else
    if SITE ~= 0
        % Take only layouts
        Sites = Sites(Sites.Layout == SITE,:);
    end
    if LightVehs == 0
        Sites = Sites(Sites.Light ~= 1,:);
    elseif LightVehs == 2
        Sites = Sites(Sites.Light == 1,:);
    end
    if Core == 1
        Sites = Sites(Sites.Core == 1,:);
    end
    if height(SitesM(strcmp(SitesM.COUNTRY,Country),:)) > 1
        Sites = Sites(strcmp(Sites.COUNTRY,Country),:);
    end
    Sitesx = Sites.SITE;

end   
    
end

