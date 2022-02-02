function [Ed, AQ, Aq] = GetBlockMaxEd_60t(Data,BlockM,Dist,ESIAT,ESIAEQ,ESIAEq,AQ1,AQ2)
%GETBLOCKMAXEd Fits, and optionally plots, BlockMaximumData
%   Data    - simply the block maximum data (max moment, shear, etc. during period)
%   BlockM  - string, 'Daily', 'Weekly', 'Monthly', 'Yearly', or 'Lifetime'
%   Dist    - string, 'Nomral, 'Lognormal'

% --- Calculate Real Design Value Ed ---

if strcmp(BlockM,'Yearly')
    n = 1;
elseif strcmp(BlockM,'Weekly')
    n = 1/50;
elseif strcmp(BlockM,'Daily')
    n = 1/(5*50);
elseif strcmp(BlockM,'Monthly')
    n = 1/12;
elseif strcmp(BlockM,'Lifetime')
    n = 50;
else
    n = BlockM;
end
    
% Beta.Yearly = 4.700; % Here we use the Beta annual - so we should use annual max effects
% Beta.Weekly = 5.444; % See Tail Fitting > Beta Conversion
% Beta.Daily = 5.724;
% Beta.Lifetime = 3.830;

Beta = norminv(1-n*0.0000013);
Alpha = 0.7;

Em = mean(Data);
Stdev = std(Data);
COV = Stdev/Em;
Delta2 = log(COV^2+1);

if strcmp(Dist,'Normal')
    Ed = Em*(1+Alpha*Beta*COV);
    % FYI ESIAT has the 1.5 in it already...
    AQ = Ed/(ESIAT);
    Aq = ((Ed/1.5)-AQ1*ESIAEQ(1)-AQ2*ESIAEQ(2))/(sum(ESIAEq));
elseif strcmp(Dist,'Lognormal')
    Ed = Em*exp(Alpha*Beta*sqrt(Delta2)-0.5*Delta2);
    AQ = Ed/(ESIAT);
    Aq = ((Ed/1.5)-AQ1*ESIAEQ(1)-AQ2*ESIAEQ(2))/(sum(ESIAEq));
elseif strcmp(Dist,'Extreme Value')
    Ed = Em*(1 + COV*(0.45 + 0.78*log(-log(normpdf(Alpha*Beta)))));    %            exp(Alpha*Beta*sqrt(Delta2)-0.5*Delta2);
    AQ = Ed/(ESIAT);
    Aq = 1;   
end

end