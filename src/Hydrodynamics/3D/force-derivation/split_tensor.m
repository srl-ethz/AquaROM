txtFile = 'C:\Users\mathi\OneDrive\Bureau\YetAnotherFEcode\external\Hydrodynamics\3D\force-derivation\04_tensor_all_parts.txt'; 
% Read file and store each line in a cell
txt = fileread(txtFile);
parts = {'','','','','','','','','','','',''}; % will store the different parts of the tensor
txtCell = strsplit(sprintf(txt), newline);
partNumber = 1;
for i=1:size(txtCell,2)
    line = strsplit(string(txtCell(i)),'],['); % separate the line if '],[' appears
    if size(line,2)~=1
        disp(i)
        parts{partNumber} = {sprintf('%s',strcat(string(parts{partNumber}),line(1)))};
        partNumber = partNumber + 1;
        parts{partNumber} = {sprintf('%s',strcat(string(parts{partNumber}),line(2)))};
    else
        parts{partNumber} = {sprintf('%s',strcat(string(parts{partNumber}),string(txtCell(i))))};
    end
end 


% Store results in 12 different files
newFile = 'C:\Users\mathi\OneDrive\Bureau\YetAnotherFEcode\external\Hydrodynamics\3D\force-derivation\05_tensor_part1.txt'; 
fid = fopen(newFile, 'wt'); 
fprintf(fid,string(parts{1}));
fclose(fid);

newFile = 'C:\Users\mathi\OneDrive\Bureau\YetAnotherFEcode\external\Hydrodynamics\3D\force-derivation\06_tensor_part2.txt'; 
fid = fopen(newFile, 'wt'); 
fprintf(fid,string(parts{2}));
fclose(fid);


newFile = 'C:\Users\mathi\OneDrive\Bureau\YetAnotherFEcode\external\Hydrodynamics\3D\force-derivation\07_tensor_part3.txt'; 
fid = fopen(newFile, 'wt'); 
fprintf(fid,string(parts{3}));
fclose(fid);


newFile = 'C:\Users\mathi\OneDrive\Bureau\YetAnotherFEcode\external\Hydrodynamics\3D\force-derivation\08_tensor_part4.txt'; 
fid = fopen(newFile, 'wt'); 
fprintf(fid,string(parts{4}));
fclose(fid);


newFile = 'C:\Users\mathi\OneDrive\Bureau\YetAnotherFEcode\external\Hydrodynamics\3D\force-derivation\09_tensor_part5.txt'; 
fid = fopen(newFile, 'wt'); 
fprintf(fid,string(parts{5}));
fclose(fid);

newFile = 'C:\Users\mathi\OneDrive\Bureau\YetAnotherFEcode\external\Hydrodynamics\3D\force-derivation\10_tensor_part6.txt'; 
fid = fopen(newFile, 'wt'); 
fprintf(fid,string(parts{6}));
fclose(fid);

newFile = 'C:\Users\mathi\OneDrive\Bureau\YetAnotherFEcode\external\Hydrodynamics\3D\force-derivation\11_tensor_part7.txt'; 
fid = fopen(newFile, 'wt'); 
fprintf(fid,string(parts{7}));
fclose(fid);

newFile = 'C:\Users\mathi\OneDrive\Bureau\YetAnotherFEcode\external\Hydrodynamics\3D\force-derivation\12_tensor_part8.txt'; 
fid = fopen(newFile, 'wt'); 
fprintf(fid,string(parts{8}));
fclose(fid);


newFile = 'C:\Users\mathi\OneDrive\Bureau\YetAnotherFEcode\external\Hydrodynamics\3D\force-derivation\13_tensor_part9.txt'; 
fid = fopen(newFile, 'wt'); 
fprintf(fid,string(parts{9}));
fclose(fid);


newFile = 'C:\Users\mathi\OneDrive\Bureau\YetAnotherFEcode\external\Hydrodynamics\3D\force-derivation\14_tensor_part10.txt'; 
fid = fopen(newFile, 'wt'); 
fprintf(fid,string(parts{10}));
fclose(fid);


newFile = 'C:\Users\mathi\OneDrive\Bureau\YetAnotherFEcode\external\Hydrodynamics\3D\force-derivation\15_tensor_part11.txt'; 
fid = fopen(newFile, 'wt'); 
fprintf(fid,string(parts{11}));
fclose(fid);

newFile = 'C:\Users\mathi\OneDrive\Bureau\YetAnotherFEcode\external\Hydrodynamics\3D\force-derivation\16_tensor_part12.txt'; 
fid = fopen(newFile, 'wt'); 
fprintf(fid,string(parts{12}));
fclose(fid);


