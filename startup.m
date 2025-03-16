clc

if ispc
    sslash = '\';
elseif isunix
    sslash = '/';
end

% Add directories
addpath(genpath(strcat(pwd,sslash,'src')));
addpath(genpath(strcat(pwd,sslash,'utils')));
addpath(genpath(strcat(pwd,sslash,'examples')));
addpath(genpath(strcat(pwd,sslash,'YetAnotherFEcode')));

disp('              _____ _____     ')
disp('  _   _  __ _|  ___| ____|___ ')
disp(' | | | |/ _` | |_  |  _| / __|')
disp(' | |_| | (_| |  _| | |__| (__ ')
disp('  \__, |\__,_|_|   |_____\___|')
disp('  |___/       YetAnotherFEcode')
fprintf('\n\n') 
