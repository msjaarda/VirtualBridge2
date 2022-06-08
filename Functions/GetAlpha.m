function Alpha = GetAlpha(BlockM)

% As of 8/6/22, Alpha depends on Reference Period...
% This is because otherwise the return period is not correct...
% See Beta Conversion spreadsheet... yearly design value corresponds to RP of 1,996 years
Alphax.Daily = 0.806;
Alphax.Weekly = 0.783;
Alphax.Monthly = 0.759;
Alphax.Yearly = 0.70; % All calibrated around the yearly according to SIA 269/1
Alphax.Lifetime = 0.512;

Alpha = Alphax.(BlockM);

end

