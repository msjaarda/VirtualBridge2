% findILNamesStr gives the names of the fields inside a structure

function [NumInfCases, ILData] = findILNamesStr(CombInfo)

% Initialize
NumInfCases = 0; ILData = struct();
ILs = fields(CombInfo);

for i = 1:length(ILs)
     
    TName = ['CombInfo.' ILs{i}];
    [NumInfCases, ILData] = recurser(NumInfCases,ILData,TName,CombInfo);
    
end
end

function [NumInfCases, ILData] = recurser(NumInfCases,ILData,TName,CombInfo)
    
    EndName = strsplit(TName,'.');
    EndName = EndName{end};
    % Detect if we are at the end
    if (EndName(1) == 'S')&&(~isempty(str2num(EndName(2:end))))
        NumInfCases = NumInfCases + 1;
        % Add to IL.Name
        ILData(NumInfCases).Name = TName;
    else
        FNames = fieldnames(eval(TName));
        for j = 1:length(FNames)
            T1Name = [TName '.' FNames{j}];
            [NumInfCases, ILData] = recurser(NumInfCases,ILData,T1Name,CombInfo);
        end
    end
end
