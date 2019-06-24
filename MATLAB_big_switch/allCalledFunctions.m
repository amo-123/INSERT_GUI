% copy all the functions called in a .m script to a folder

% main .m file
main = 'big_switch.m';

% destination folder
destination = 'C:\Users\Lab FIORINI\Desktop\per mediso\MedisoFnc';

% all the functions called in main
fList = matlab.codetools.requiredFilesAndProducts(main);

for i=1:length(fList)
    copyfile(fList{i} , destination);
end
