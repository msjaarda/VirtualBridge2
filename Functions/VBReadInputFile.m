function [BaseData] = VBReadInputFile(InputFile)
% Reads InputFile (VBInput) and returns BaseData for simulations

% Get all sheetNames in InputFile (should only be 1)
[~, sheetNames] = xlsfinfo(InputFile);

% Convert sheetNames from Cell Array (sheetNames) to String Array (names)
% Initialize names                 Populate names
names = strings(size(sheetNames)); [names{:}] = sheetNames{:};

% Assign BaseData as first sheet
BaseData = readtable(InputFile,'Sheet',names(1));

% We add defaults to BaseData, so that we don't have to include them
% everytime in the InputFile
BaseData = BaseDataDefaults(BaseData);

% Start with Folder
% Folder is optional, but when not included (or when given as 0), must be '/'
if ismember('Folder', BaseData.Properties.VariableNames)
    if iscell(BaseData.Folder)
        for i = 1:height(BaseData)
            if isempty(findstr('/',BaseData.Folder{i}))
                if length(BaseData.Folder{i}) > 1
                    BaseData.Folder{i} = ['/' BaseData.Folder{i}];
                else
                    BaseData.Folder{i} = ['/'];
                end
            end
        end
    else
        % This '/' is given so that the filepath still reads correctly with no folder
        clear BaseData.Folder
        BaseData.Folder = cell(height(BaseData),1);
        BaseData.Folder(:) = {'/'};
    end
else
    % This '/' is given so that the filepath still reads correctly with no folder
    clear BaseData.Folder
    BaseData.Folder = cell(height(BaseData),1);
    BaseData.Folder(:) = {'/'};
end

warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary')

end

