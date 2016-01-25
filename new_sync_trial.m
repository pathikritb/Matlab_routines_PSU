% Create a giant binary file with 8 bytes time
%% Initialize variables.
clear
clc
runname = 'p4458';
acname = ['./' runname '_acoustics'];
recname = 'CF2';
recnumber = str2double(recname(3)); % number used to save data (e.g. 2 for control file 2)
if strcmp(runname,'p4457')
    CF1 = 232215:536856;%1500000; % CF1 for p4457, adjusted for usable data
    CF2 = 606891:1527209; % CF2 for p4457
    idxmin_peak_vec = [396248,673653]; % central peaks in sync for each control file
elseif strcmp(runname,'p4458')
    CF1 = 208309:1233694;%1500000; % CF1 for p4458, adjusted for usable data
    CF2 = 1311457:2764167; % CF2 for p4458
    CF3 = 2822024:3804488; % CF3 for p4458
    idxmin_peak_vec = [1065560,1674979,3260187]; % central peaks in sync for each control file
end
syncind = eval(recname);
%[LP_Disp,Shr_stress,nor_disp,Nor_stress,Time,CenBlk_slip,sync,Samp_Freq,mu,VarName10] = importfile(filename);
[data,outname] = ReadBinBiax(runname); % Path: Read directly from binary file
indlim   = [1 length(data)];% After 1526272, the data is useless!!
ind      = indlim(1):indlim(2);
if strcmp(runname,'p4457')
    siz      = size(data);
    biaxcol  = data(ind,1);
    slip_lp  = data(ind,2);
    Time     = data(ind,6);
    mu       = data(ind,10);
    sigma    = data(ind,5);
    closure  = data(ind,4);
    cent_blk = data(ind,7);
    sync     = data(ind,8);
    Samp_Freq = data(ind,9);
    shr_stress = data(ind,3);
elseif strcmp(runname,'p4458') || strcmp(runname,'p4459')
    siz      = size(data);
    biaxcol  = data(ind,1);
    slip_lp  = data(ind,2);
    time     = data(ind,6);
    mu       = data(ind,10);
    sigma    = data(ind,5);
    closure  = data(ind,4);
    cent_blk = data(ind,7);
    surf_slip = data(ind,8);
    shr_stress = data(ind,3);
    Time     = data(:,6);
    sync     = data(ind,9);
    Samp_Freq = [time(2)-time(1);diff(time)];
end
% CF2 (indexes for Control File 2)
idx1 = syncind(1)-indlim(1)+1; % 208309 for CF1, 1311457 for CF2, 2822024 for CF3; index corresponding to first trigger of the sequence on figure 1
idx2 = syncind(end)-indlim(1)+1; % 1233694 for CF1, 2764167 for CF2, 3804488 for CF3; index corresponding to last trigger of the sequence on figure 1
idxmin_peak = idxmin_peak_vec(recnumber)-indlim(1)+1; % 1065560 for CF1, 1674979 for CF2, 3260187 for CF3; index corresponding to a large trigger in the middle of the sequence
filenamedata = ['rec' num2str(recnumber) '.mat']; % filename of the mat file

N = idx2 - idx1 + 1;
if strcmp(runname,'p4457')
    WFlength = 8192/2;% p4457: 8192, but only half saved
    numberofWFperfile = 640; % p4457: 640, p4458: 2560
else
    WFlength = 2048;% p4458: 2048
    numberofWFperfile = 2560; % p4457: 640, p4458: 2560
end
acRate = 1000; % in Hz, number of WF per second
acPeriod = 1/acRate;
% mechanical data
mechtime = Time(idx1:idx2);
shearstress = shr_stress(idx1:idx2);
sync1 = sync(idx1:idx2);
 
% normalize sync between 0 and 1
MAX = max(sync1);
MIN = min(sync1);
AMP = MAX-MIN;
sync1 = (sync1-mean(sync1))/AMP;
sync1 = sync1 - max(sync1);
idealsync = zeros(N,1);
t_btw_trig = numberofWFperfile/acRate*1000; % time in ms

% build a trigger time vector using the large trigger chosen above as a reference
trigger_time1 = fliplr(Time(idxmin_peak):-t_btw_trig/1000:Time(idx1));
trigger_time2 = Time(idxmin_peak):t_btw_trig/1000:Time(idx2);

% the sample corresponding to idxmin_peak is both in 1 and 2 so we remove it from
% 1. We also remove the last one from 2 because the last trigger is sent
% but no data is actually taken (verasonics was freezed meanwhile)
trigger_time = [trigger_time1(1:end-1) trigger_time2(1:end-1)]; 

check = diff(trigger_time); % check time between triggers (should be constant...)

% this number should equal the number of acoustic files
Ntrigger = length(trigger_time);
totalnumberoffiles = Ntrigger;display(totalnumberoffiles)
acN = totalnumberoffiles*numberofWFperfile; % total number of waveforms (WF)
trigs = -0.5*ones(Ntrigger,1); % -0.15 to be adjusted depending on the experiment (adjust it see both signals clearly on figure 2)

biax_num_ac = cell(Ntrigger-1,1);
num_elem_biax_ac = ones(Ntrigger,1);
parfor i = 1:Ntrigger-1
    i
    biax_num_ac{i} = find(Time>=trigger_time(i) & Time<trigger_time(i+1));
    num_elem_biax_ac(i) = length(biax_num_ac{i});
end
time_last    = trigger_time(Ntrigger)+(numberofWFperfile)/acRate;
%biax_num_vec = unique(biax_num_vec);
biax_num_ac{Ntrigger} = find(Time>=trigger_time(Ntrigger) & Time<=time_last);
num_elem_biax_ac(Ntrigger) = length(biax_num_ac{Ntrigger});
biax_num_vec = vertcat(biax_num_ac{:});
hh = 1;
WF = 0; % init waveform
%TimeTrig = acPeriod*stackN*ones(1,numberofWFperfile/stackN);
fileID = fopen([acname '/' recname '/' recname '_seis_v2.comb'],'w');
fclose(fileID);
fileID = fopen([acname '/' recname '/' recname '_seis_v2.comb'],'a');
ncolACmat = max(num_elem_biax_ac);
ACmat  = zeros(WFlength+2,ncolACmat);
for ii = 1:1:totalnumberoffiles % normally it should be "1:1:totalnumberoffiles" here (I had only 5000 files on my computer...)
    ii
    ACmat  = zeros(WFlength+2,ncolACmat);
    TimeTrig   = Time(biax_num_ac{ii});
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
%     if stackN > 1
%         for jj = 1:numberofWFperfile/stackN
%             ACmat(3:WFlength+1,jj) = sum(ACdat(1:WFlength,stackN*(jj-1)+1:stackN*jj)/stackN,2);
%         end
%     else
%         ACmat(3:WFlength+1,:) = ACdat;
%     end
    fclose(fid);   
    ACmat(:,num_elem_biax_ac(ii)+1:ncolACmat) = [];
    ACmat(1,:)   = TimeTrig;
    ACmat(2,:)   = biax_num_ac{ii};
    fwrite(fileID,ACmat,'single');
end
% build acoustic time vector (after staking 10, i.e. 1 sample every 10ms)
%acTime = 0:stackN*acPeriod:stackN*acPeriod*(acN/stackN-1); % /10 because of stacking 10 waveforms
%acTime = acTime + Time(idx1);
acTime = Time(biax_num_vec);