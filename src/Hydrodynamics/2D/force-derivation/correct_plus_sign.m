txtFile = 'C:\Users\mathi\OneDrive\Bureau\YetAnotherFEcode\external\Hydrodynamics\2D\force-derivation\00_tensor_initial.txt'; 
% Read file and store each line in a cell
txt = fileread(txtFile);
txtCell = strsplit(txt, newline);
% Find rows that start with '+'
hasPattern = regexp(txtCell, '^\s*\+.*'); 
rowIdx = find(~cellfun(@isempty, hasPattern))  %row numbers that begin with '+'
% Replace the rows with new text
content = string(txtCell(rowIdx));
modifiedInitial = erase(string(txtCell(rowIdx-1))," ..."); % erase the space and the ... at the end of the line
modifiedInitial = strtrim(modifiedInitial); % delete the backspace
combinedContent = strcat(modifiedInitial,content); % combine the two parts to form the updated line
for i=1:size(rowIdx,2)
    txtCell(rowIdx(:,i)-1) = {sprintf('%s',combinedContent(i))};
end
txtCell(rowIdx)  = []; 

%txtCell = cellfun(@(x)sprintf('%s\n',x), txtCell, 'UniformOutput', false); %add 'newrow' char
% Write text to new file (or you could overwrite the old one).
newFile = 'C:\Users\mathi\OneDrive\Bureau\YetAnotherFEcode\external\Hydrodynamics\2D\force-derivation\01_tensor_plus_sign_corrected.txt'; 
fid = fopen(newFile, 'wt'); 
fprintf(fid,[txtCell{:}]);
fclose(fid);



