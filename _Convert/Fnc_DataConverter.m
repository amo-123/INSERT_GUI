function [dirOut,fileOut]=Fnc_DataConverter(fileIn,WorkDir,Head,Type)
tic;

extension=fileIn((end-3):end);
if ~strcmp(extension,'data')
    disp('Improper extension. File of ''data'' extension is expected.');
    return;
end

if ~isdir([WorkDir,'\Corr\digitalraw'])
    mkdir([WorkDir,'\Corr\digitalraw']);
end

dirOut=[WorkDir,'\Corr\digitalraw\'];
fileOut=[dirOut,'digitalraw_H',num2str(Head,'%02i'),'_',Type,'.dat'];

Count = MexDataConverter(regexprep(fileIn,'\','/'),...
                         regexprep(fileOut,'\','/'));


disp([ 'Number of data packet converted: ', num2str(Count) ]);


EllapsedTime = toc;
 
disp([ 'Ellapsed time: ', num2str(EllapsedTime) ]);

end