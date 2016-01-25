% Function for reading the binary file from the biax, the Marone lab calls
% this the lk file
function [data,outname] = ReadBinBiax(runname)
% Input: runname, a string of the form 'pxxxx'. The code expects the data
% file to be in foo/pxxx/data/pxxx_data.bin. Edit the variable mechfname
% and the filename argument for fidbiax to suit your needs.
% Output: 1) data: data matrix with all relevant columns appended with the biax
% column as the first column (refer line 65), 
%         2) outname: header info of the lk file. Edit the outname variable
%         to change to your choice of file structure.
% Code copies John Leeman's python code
mechfname = ['./' runname '/data/' runname]; % Assumes a structure foo/pxxx/data/pxxx_data.bin for the look file. Change it to whatever you like to suit your file structure
outname   = ['./' runname '/data/' runname '.hdr']; % Outputs header in a structure foo/pxxx/data/pxxx.hdr. Again, change it to whatever you like to suit your file structure.
fidbiax   = fopen([mechfname  '_data.bin'],'rb');% Note the foo_data.bin extension, either change this or change the file name for the lk file.
fidheader = fopen(outname,'w');
% Unpack information at the top of the file about the experiment
name = fread(fidbiax,[1 20],'uint8'); % First Line is the name, 20 bytes name
ind0        = find(name==0,1,'first'); % Where is the first blank space encountered
name(ind0:end) = []; % Remove everything after the first blank space is encountered 
fprintf(fidheader,'Filename: %20s\n',name);
% The rest of the header information is written in big endian format
% Number of records (32 bit integer)
num_recs = fread(fidbiax,1,'uint','b');
fprintf(fidheader,'Number of records: %d\n',num_recs);
% Number of columns (32 bit integer)
num_cols = fread(fidbiax,1,'uint','b');
fprintf(fidheader,'Number of columns: %d\n',num_cols);
% Sweep (int) - No longer used
swp =  fread(fidbiax,1,'uint','b');
fprintf(fidheader,'Swp : %d\n',swp);
% Date/time(int) - No longer used
dtime = fread(fidbiax,1,'uint','b');
fprintf(fidheader,'dtime : %d\n',dtime);
% For each possible column (32 maximum columns) unpack its header
% information and store it.  Only store column headers of columns
% that contain data.  Use termination at first NUL.
chname = cell(32,1);
chunits = cell(32,1);
nelem   = zeros(32,1);
for i=1:32        
        % Channel name (13 characters)
        chname_temp = fread(fidbiax,[1 13],'uint8'); % Temporary names to correct weird problems with space encoding
        ind0        = find(chname_temp==0,1,'first'); % Where is the first blank space encountered
        chname{i}   = char(chname_temp(1:ind0-1)); % Remove everything after the first blank space is encountered      
        % Channel units (13 characters)
        chunits_temp = fread(fidbiax,[1 13],'uint8'); % Temporary names to correct weird problems with space encoding
        ind0        = find(chunits_temp==0,1,'first'); % Where is the first blank space encountered
        chunits{i}   = char(chunits_temp(1:ind0-1)); % Remove everything after the first blank space is encountered 
        % This field is now unused, so we just read past it (int)
        gain = fread(fidbiax,1,'uint','b');
        % This field is now unused, so we just read past it (50 characters)
        comment = fread(fidbiax,[1 50],'uint8=>char');
        % Number of elements (int)
        nelem(i) = fread(fidbiax,1,'uint','b');
end
% Remove blank columns
chname(nelem==0,:)=[]; % Strange deletion behavior for matlab cells 
chunits(nelem==0,:)=[]; % Some more strange behavior
nelem(nelem==0)=[]; % Back to less strange strange
% Print column units and headings in the header file
fprintf(fidheader,'\n\n----------------------------------------------------------------------------\n');
fprintf(fidheader,'|%6s|%17s|%17s|%17s|\n','Column','Name','Unit','Records');
        fprintf(fidheader,'----------------------------------------------------------------------------\n');
for i = 1:length(nelem)
    fprintf(fidheader,'|%6s|%17s|%17s|%17s|\n',num2str(i),chname{i},chunits{i},num2str(nelem(i)));
end
fprintf(fidheader,'----------------------------------------------------------------------------\n');
% Read the actual data from file
data = fread(fidbiax,[nelem(1) length(nelem)],'double');
data = [(1:nelem(1))' data]; % Append the biax column count as the first column

fclose(fidbiax);
fclose(fidheader);