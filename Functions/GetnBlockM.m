function n = GetnBlockM(BlockM)

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

end

