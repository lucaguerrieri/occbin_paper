% uncomment a location below -- adapt to local installation

 %location = 'home_matteo';
 location = 'home_luca';
 %location = 'work_luca';
 %location = 'work_luca2';
 
%  location = 'work_matteo';

restoredefaultpath

if strmatch(location,'home_matteo','exact')    
    dir1='C:\U\TQS\luca\dynare4\4.2.1\matlab';
    dir2='C:\E\occbin\occbin_20121219\toolkit_files';
    dir3='C:\E\occbin\occbin_20121219\utils_files';
elseif strmatch(location,'home_luca','exact')
    dir1='/users/jason/documents/matlab/dynare/4.1.0/matlab';
    dir2='../toolkit_files';
    dir3='../utils_files';
elseif strmatch(location,'work_luca','exact')
    dir1='G:\ofs\research2\Guerrieri\matlab\Dynare4\4.2.4\matlab';
    dir2='..\toolkit_files';
    dir3='..\utils_files';
elseif strmatch(location,'work_luca2','exact')
    dir1='G:/ofs/research2/Guerrieri/matlab/Dynare4/4.2.4/matlab';
    dir2='../toolkit_files';
    dir3='../utils_files';
elseif strmatch(location,'work_matteo','exact')  
    dir1='U:\TQS\luca\dynare4\4.2.4\matlab';
    dir2='C:\E\Consumption\borrcon_dyn\toolkit_files';
else 
    error('Specify path to Dynare installation')
end

path(path,dir1);
path(path,dir2);
path(path,dir3);

dynare_config

