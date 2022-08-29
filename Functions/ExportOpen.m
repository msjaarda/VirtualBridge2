function ExportOpen
%EXPORTOPEN 
%   Exports all open figures and saves them according to their titles

FolderName = 'C:\Users\mjsja\Desktop';   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');

j = length(FigList);
for i = 1:length(FigList)
    FigHandle = FigList(i);
    FigName = get(FigHandle, 'Name');
    TName = datestr(now+i/86400,'mmmdd-yy HHMMSS');
    %savefig(FigHandle, fullfile(FolderName, FigName, '.fig'));
%     if isempty(FigName) || contains(FigName,'Figure')
        exportgraphics(FigHandle,[FolderName '\' TName '.jpg'],'Resolution',600);
%     else
%         try
%             exportgraphics(FigHandle,[FolderName '\' FigName '.jpg'],'Resolution',600);
%         catch
%             exportgraphics(FigHandle,[FolderName '\Figure' num2str(j) '.jpg'],'Resolution',600);
%         end
%     end
    j = j-1;
end

end

