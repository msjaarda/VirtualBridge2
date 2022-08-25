% findIL gets influence lines at the end of all branches

function [NumInfCases, ILData] = findIL(ILs,ILRes,NumLanes)

% Initialize
NumInfCases = 0; ILData = struct();

for i = 1:length(ILs)
     
    TName = ['ILLib.' ILs{i}];
    [NumInfCases, ILData] = recurser(NumInfCases,ILData,TName,ILRes,NumLanes);
    
end
end

function [NumInfCases, ILData] = recurser(NumInfCases,ILData,TName,ILRes,NumLanes)
    load Misc\ILLib.mat
    % Detect if we are at the end
    if isnumeric(eval([TName '(:,1)']));
        NumInfCases = NumInfCases + 1;
        ILx = eval([TName '(:,1)']);
        ILv = eval([TName '(:,2:end)']);
        % Round to refinement level (ILRes)
        RoundedILx = ILx(1):ILRes:ILx(end); RoundedILx = RoundedILx';
        % Now we interpolate the influence lines and populate IL.v
        ILData(NumInfCases).v = interp1(ILx,ILv,RoundedILx);
        % Add to IL.Name
        ILData(NumInfCases).Name = TName;
        if size(ILData(NumInfCases).v,2) < NumLanes
            % Can choose to duplicate, or add zeros. The rule will be that
            % if size(ILData.v{Num.InfCases},2) == 1, we duplicate, but if
            % it is greater, we add zeros with a notification
            if size(ILData(NumInfCases).v,2) == 1
                ILData(NumInfCases).v = repmat(ILData(NumInfCases).v,1,NumLanes);
            else
                for t = size(ILData(NumInfCases).v,2) + 1:NumLanes
                    ILData(NumInfCases).v(:,t) = 0;
                    % Turn back on later... was annoying!
                    %fprintf('\nWARNING Lane mismatch for IL: %s, ILs with zeros added',TName)
                end
            end
        end
    else
        FNames = fieldnames(eval(TName));
        for j = 1:length(FNames)
            T1Name = [TName '.' FNames{j}];
            [NumInfCases, ILData] = recurser(NumInfCases,ILData,T1Name,ILRes,NumLanes);
        end
    end
end
