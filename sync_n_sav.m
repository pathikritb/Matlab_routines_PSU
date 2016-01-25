% Create a giant binary file with 8 bytes time, save it for memmapfile
% access. Borrows from Jacques Riviere's sync codes
% Produces synced times, biax columns and corresponding (stacked if
% appropriate) waveforms.
%% Initialize variables.
function outname = sync_n_sav(runname,recname,idxfile_path,out)
% Inputs: 
% (1) runname = 'p4458', the runname of the file, will be used to gather
%       biax data and acoustic data. Put all the acoustic data in a folder called /your_path/runname_acoustics
%       or change code accordingly for your own directory structure
% (2) recname = 'CF2', the control file name, expected path /your_path/runname_acoustics/recname
% (3) idxfile_path = $path for file with idx for CF syncing: expected a file
%       with contents like:
%       %CF name %1st ind (first index from sync record) %2nd ind (last index from sync record) %min ind (index for minimum sync record)
%       CF1 232215  536856  396248
%       CF2 606891  1527209 673653
%       If idxfile_path is absent or []: default path = /your_path/runname_acoustics/index_CF.txt
% (4) out = $path for output file, default : /your_path/runname_acoustics/recname/
% Output:
% outname: $path for output file
%%---------------------------------------------------------------------------------------------------%%
% Filename and argcheck
acname = ['./' runname '_acoustics']; % setting path to acoustics file
if (nargin<3)||isempty(idxfile_path)
    idxfile_path = ['./' runname '_acoustics/index_CF.txt'];
    out          = [acname '/' recname '/' recname '_seis_git.comb'];
elseif (nargin<4)||isempty(out)
    out          = [acname '/' recname '/' recname '_seis_git.comb'];
end
if isempty(idxfile_path)
        idxfile_path = ['./' runname '_acoustics/index_CF.txt'];
    elseif isempty(out)
        out      = [acname '/' recname '/' recname '_seis_git.comb'];
end

% Processing begins
recnumber = str2double(recname(3)); % number used to save data (e.g. 2 for control file 2)
fID = fopen(idxfile_path);
C = textscan(fID,'%s %u32 %u32 %u32','HeaderLines',1);
fclose(fID);
idx1    = C{1,2}(recnumber); % Set first  biax_col index
idx2    = C{1,3}(recnumber);  % Set last  biax_col index
idxmin_peak = C{1,4}(recnumber); % Set biax_col index for minimum sync record
% Read in Biax Data
%[LP_Disp,Shr_stress,nor_disp,Nor_stress,Time,CenBlk_slip,sync,Samp_Freq,mu,VarName10] = importfile(filename);
[data,~] = ReadBinBiax(runname); % Path: Read directly from binary file
Time     = data(:,6);
if strcmp(runname,'p4457')
    WFlength = 8192/2;% p4457: 8192, but only half saved
    numberofWFperfile = 640; % p4457: 640, p4458: 2560
elseif strcmp(runname,'p4458') || strcmp(runname,'p4459')
    WFlength = 2048;% p4458: 2048
    numberofWFperfile = 2560; % p4457: 640, p4458: 2560
else
    WFlength = 2048;% p4458: 2048
    numberofWFperfile = 2560; % p4457: 640, p4458: 2560
end

% Syncing acoustic data
N = idx2 - idx1 + 1;
acRate = 1000; % in Hz, number of WF per second
t_btw_trig = numberofWFperfile/acRate*1000; % time in ms

% build a trigger time vector using the large trigger chosen above as a reference
trigger_time1 = fliplr(Time(idxmin_peak):-t_btw_trig/1000:Time(idx1));
trigger_time2 = Time(idxmin_peak):t_btw_trig/1000:Time(idx2);

% the sample corresponding to idxmin_peak is both in 1 and 2 so we remove it from
% 1. We also remove the last one from 2 because the last trigger is sent
% but no data is actually taken (verasonics was freezed meanwhile)
trigger_time = [trigger_time1(1:end-1) trigger_time2(1:end-1)]; 
% this number should equal the number of acoustic files
Ntrigger = length(trigger_time);
totalnumberoffiles = Ntrigger;display(totalnumberoffiles)
biax_num_ac = cell(Ntrigger,1); % Each .ac file corresponds to some number of biax columns
num_elem_biax_ac = ones(Ntrigger,1); % How many biax_columns in the particular file
parfor i = 1:Ntrigger-1
    i
    biax_num_ac{i} = find(Time>=trigger_time(i) & Time<trigger_time(i+1));
    num_elem_biax_ac(i) = length(biax_num_ac{i});
end
time_last    = trigger_time(Ntrigger)+(numberofWFperfile)/acRate;
%biax_num_vec = unique(biax_num_vec);
biax_num_ac{Ntrigger} = find(Time>=trigger_time(Ntrigger) & Time<=time_last);
num_elem_biax_ac(Ntrigger) = length(biax_num_ac{Ntrigger});
%TimeTrig = acPeriod*stackN*ones(1,numberofWFperfile/stackN);
fileID = fopen(out,'w');
fclose(fileID);
fileID = fopen(out,'a');
ncolACmat = max(num_elem_biax_ac);% How many columns maximum needed for writing? Maximum number of biax_cols in a WF file
for ii = 1:1:totalnumberoffiles % normally it should be "1:1:totalnumberoffiles" here (I had only 5000 files on my computer...)
    ii
    ACmat  = zeros(WFlength+2,ncolACmat);
    TimeTrig   = round(Time(biax_num_ac{ii})*1000)/1000;
    ACfilename = [acname '/CF' num2str(recnumber) '/WF_CF' num2str(recnumber) '_' num2str(ii) '.ac'];
    fid = fopen(ACfilename,'r');
    if strcmp(runname,'p4457')
        ACdat = fread(fid,[WFlength*2,numberofWFperfile],'int16');
    else
        ACdat = fread(fid,[WFlength,numberofWFperfile],'int16');
    end
    indstart  = 1;
    for jj = 1:num_elem_biax_ac(ii)
        if jj < num_elem_biax_ac(ii)
            diff_time = TimeTrig(jj+1) - TimeTrig(jj);
            nstack    = round(diff_time*acRate);
            ACmat(3:WFlength+2,jj) = sum(ACdat(1:WFlength,indstart:indstart+nstack-1)/nstack,2);
            indstart  = indstart+nstack;
        else
            nstack = length(indstart:numberofWFperfile);
            ACmat(3:WFlength+2,jj) = sum(ACdat(1:WFlength,indstart:indstart+nstack-1)/nstack,2); 
        end      
    end
    fclose(fid);   
    ACmat(:,num_elem_biax_ac(ii)+1:ncolACmat) = [];
    ACmat(1,:)   = TimeTrig;
    ACmat(2,:)   = biax_num_ac{ii};
    fwrite(fileID,ACmat,'single');% Each ACmat: 1st row: Times, 2nd row: Biax Column, 3rd to (WFlength+2)th row : Waveform Vector
    % ACmat(:,i) = [Time_of_WF(i);Biax_Column(i);WF vector(i)]
end
outname = out;