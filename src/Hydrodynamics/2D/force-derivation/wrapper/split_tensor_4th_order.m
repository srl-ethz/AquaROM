txtFile = 'C:\Users\mathi\OneDrive\Bureau\make fishes swim\YetAnotherFEcode\external\Hydrodynamics\force-derivation\04_tensor_all_parts.txt'; 
% Read file and store each line in a cell
txt = fileread(txtFile);
parts = {'','','','','',''}; % will store the different parts of the tensor
txtCell = strsplit(sprintf(txt), newline);
partNumber = 1;
for i=1:size(txtCell,2)
    line = strsplit(string(txtCell(i)),']],[['); % separate the line if '],[' appears
    if size(line,2)~=1
        disp(i)
        parts{partNumber} = {sprintf('%s',strcat(string(parts{partNumber}),line(1)))};
        partNumber = partNumber + 1;
        parts{partNumber} = {sprintf('%s',strcat(string(parts{partNumber}),line(2)))};
    else
        parts{partNumber} = {sprintf('%s',strcat(string(parts{partNumber}),string(txtCell(i))))};
    end
end 


% Store results in 6 different files
newFile = 'C:\Users\mathi\OneDrive\Bureau\make fishes swim\YetAnotherFEcode\external\Hydrodynamics\force-derivation\05_4th_tensor_part1.txt'; 
fid = fopen(newFile, 'wt'); 
fprintf(fid,string(parts{1}));
fclose(fid);

newFile = 'C:\Users\mathi\OneDrive\Bureau\make fishes swim\YetAnotherFEcode\external\Hydrodynamics\force-derivation\06_4th_tensor_part2.txt'; 
fid = fopen(newFile, 'wt'); 
fprintf(fid,string(parts{2}));
fclose(fid);


newFile = 'C:\Users\mathi\OneDrive\Bureau\make fishes swim\YetAnotherFEcode\external\Hydrodynamics\force-derivation\07_4th_tensor_part3.txt'; 
fid = fopen(newFile, 'wt'); 
fprintf(fid,string(parts{3}));
fclose(fid);


newFile = 'C:\Users\mathi\OneDrive\Bureau\make fishes swim\YetAnotherFEcode\external\Hydrodynamics\force-derivation\08_4th_tensor_part4.txt'; 
fid = fopen(newFile, 'wt'); 
fprintf(fid,string(parts{4}));
fclose(fid);


newFile = 'C:\Users\mathi\OneDrive\Bureau\make fishes swim\YetAnotherFEcode\external\Hydrodynamics\force-derivation\09_4th_tensor_part5.txt'; 
fid = fopen(newFile, 'wt'); 
fprintf(fid,string(parts{5}));
fclose(fid);

newFile = 'C:\Users\mathi\OneDrive\Bureau\make fishes swim\YetAnotherFEcode\external\Hydrodynamics\force-derivation\10_4th_tensor_part6.txt'; 
fid = fopen(newFile, 'wt'); 
fprintf(fid,string(parts{6}));
fclose(fid);


