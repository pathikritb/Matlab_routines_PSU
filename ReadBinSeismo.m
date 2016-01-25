function[seismo,size] = ReadBinSeismo(fname)
%% fname is a character array
fid = fopen(fname,'rb');
size = fread(fid,2,'int');
seismo = fread(fid,size(2),'double');
fclose(fid);