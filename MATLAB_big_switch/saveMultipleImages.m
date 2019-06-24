%FolderName = tempdir;   % Your destination folder

FolderName = 'C:\Users\INSERT\Desktop\Measurements\20180518\figure\';
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');
  saveas(FigHandle, [FolderName, [FigName,'_',num2str(iFig)], '.png'], 'png')
end

