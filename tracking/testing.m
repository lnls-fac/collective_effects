path = '/opt/tracySirius';
arqs = [];
paths = [];
while true
   [file path] = uigetfile('*','Escolha um arquivo:',path,'MultiSelect','on');
   if isequal(file,0), break; end
   if ~iscell(file), file2 = {file}; else file2 = file; end
   path2 = {path};
   paths = [paths repmat(path2(:),1,length(file2))];
   arqs = [arqs file2];
end