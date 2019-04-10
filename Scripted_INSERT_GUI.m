%folder =uigetdir();

%files = dir(fullfile(folder, '*.data'));

a = load('batchinconv.mat');
b = load('batchtoconvhandle.mat');
%%
f1=struct2cell(a.handles);
f2=struct2cell(b.handles);
%%
f=isequal(f1,f2);
%%
[ia,ib] = setdiff(f1,f2);
%%
a=rmfield(a,f);
b=rmfield(b,f);
