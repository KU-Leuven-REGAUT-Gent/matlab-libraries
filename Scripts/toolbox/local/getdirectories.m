% %% PREFIX PROJECTS DIRECTORY
 prefix{1}='%LIB:Projects%';
 prefix{2}='D:\Users\Dimitri\Documenten\OneDrive - KU Leuven\2. Collegas\Projects\MATLAB';
 prefix{3}='D:\Users\Dimitri\Documenten\Libraries';
% prefix{4}='E:\Users\Frederic\Projects';
% 
% %% MATLAB PROJECTS LIBRARY DIRECTORY
 librarydirectory = find_dir(prefix, '\MATLAB\');
% 
% %% MATLAB PROJECTS ROOT DIRECTORY
 global rootdirectory
 rootdirectory = find_dir(prefix, '\MATLAB');
% 
% %% CLEANING UP
 clear prefix;

