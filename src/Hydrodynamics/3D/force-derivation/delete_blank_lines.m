txtFile = 'C:\Users\mathi\OneDrive\Bureau\make fishes swim\YetAnotherFEcode\external\Hydrodynamics\force-derivation\original_tensor.txt'; 
% Read file and store each line in a cell
txt = fileread(txtFile);
txtCell = strsplit(txt, newline)'; 
% Find rows that start with '# ####'
hasPattern = regexp(txtCell, '\\\s');          %this searches for "#' followed by 1 space and at least 1 number.
rowIdx = find(~cellfun(@isempty, hasPattern));  %row numbers that will be replaced
% Replace the rows with new text
newText = strsplit(sprintf(''));   %Here's where you create the new text; ignore last element. 
txtCell(rowIdx)  = newText(1); 

%txtCell = cellfun(@(x)sprintf('%s\n',x), txtCell, 'UniformOutput', false); %add 'newrow' char
% Write text to new file (or you could overwrite the old one).
newFile = 'C:\Users\mathi\OneDrive\Bureau\make fishes swim\YetAnotherFEcode\external\Hydrodynamics\force-derivation\wo_blank_lines_tensor.txt'; 
fid = fopen(newFile, 'wt'); 
fprintf(fid, [txtCell{setdiff(1:end,[rowIdx,rowIdx+1])}]);
fclose(fid);