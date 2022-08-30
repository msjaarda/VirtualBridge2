function Name = ConvertLayoutName(Layout)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if Layout == 11
    Name = 'Uni2L';
elseif Layout == 111
    Name = 'Uni3L';
elseif Layout == 12
    Name = 'Bi2L';
elseif Layout == 1122
    Name = 'Bi4L';
elseif Layout == 110
    Name = 'LSVAUni2L';
else
    load('Sites.mat')
    Name = strcat(Sites.SName(Sites.SITE == Layout),num2str(Layout));
end

                                    
end

