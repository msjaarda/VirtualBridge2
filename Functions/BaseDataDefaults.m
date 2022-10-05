function BaseData = BaseDataDefaults(BaseData)
%BaseDataDefaults

if ~ismember('RunFat', BaseData.Properties.VariableNames)
    BaseData.RunFat(:) = 0;
end
if ~ismember('Apercu', BaseData.Properties.VariableNames)
    BaseData.Apercu(:) = 0;
end
if ~ismember('Analysis', BaseData.Properties.VariableNames)
    BaseData.Analysis(:) = 1;
end
if ~ismember('BunchFactor', BaseData.Properties.VariableNames)
    BaseData.BunchFactor(:) = 1;
end
if ~ismember('Parallel', BaseData.Properties.VariableNames)
    BaseData.Parallel(:) = 0;
end
if ~ismember('RunDyn', BaseData.Properties.VariableNames)
    BaseData.RunDyn(:) = 1;
end
if ~ismember('Save', BaseData.Properties.VariableNames)
    BaseData.Save(:) = 1;
end
if ~ismember('NumAnalyses', BaseData.Properties.VariableNames)
    BaseData.NumAnalyses(:) = 1;
end
if ~ismember('ClassType', BaseData.Properties.VariableNames)
    BaseData.ClassType(:) = "All";
end
if ~ismember('VWIM', BaseData.Properties.VariableNames)
    BaseData.VWIM(:) = 0;
end
if ~ismember('AnalysisType', BaseData.Properties.VariableNames)
    BaseData.AnalysisType(:) = "Sim";
end
if ~ismember('Plots', BaseData.Properties.VariableNames)
    BaseData.Plots(:) = false;
end
if ~ismember('Period', BaseData.Properties.VariableNames)
    BaseData.Period(:) = "Yearly";
end
% Can remove?
if ~ismember('LSVA', BaseData.Properties.VariableNames)
    BaseData.LSVA(:) = false;
end
if ~ismember('Stage2P', BaseData.Properties.VariableNames)
    BaseData.Stage2P(:) = false;
end
if ~ismember('StopSim', BaseData.Properties.VariableNames)
    BaseData.StopSim(:) = false;
end
% Can remove?
if ~ismember('Detailed', BaseData.Properties.VariableNames)
    BaseData.Detailed(:) = false;
end
if ~ismember('Plat', BaseData.Properties.VariableNames)
    BaseData.Plat(:) = false;
end
if ~ismember('Country', BaseData.Properties.VariableNames) % Anything that doesn't produce a match does nothing... can do "ALL COUNTRIES" for the sake of plots!
    BaseData.Country(:) = "CH";
end
if ~ismember('LightVehs', BaseData.Properties.VariableNames) % 0 is not included, 1 is included, 2 is ONLY
    BaseData.LightVehs(:) = false;
end
if ~ismember('PUN', BaseData.Properties.VariableNames)
    BaseData.PUN(:) = 0;
end

end

