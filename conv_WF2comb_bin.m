% Create a giant binary file with 8 bytes time
%% Initialize variables.
runname = 'p4457';
acname = ['./' runname '_acoustics'];
recname = 'CF1';

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
    Shr_stress = data(ind,3);
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
idx1 = 232215-indlim(1)+1; % 208309 for CF1, 1311457 for CF2, 2822024 for CF3; index corresponding to first trigger of the sequence on figure 1
idx2 = 536856-indlim(1)+1; % 1233694 for CF1, 2764167 for CF2, 3804488 for CF3; index corresponding to last trigger of the sequence on figure 1
idxmin_peak = 396248-indlim(1)+1; % 1065560 for CF1, 1674979 for CF2, 3260187 for CF3; index corresponding to a large trigger in the middle of the sequence
recnumber = 1; % number used to save data (e.g. 2 for control file 2)
filenamedata = ['rec' num2str(recnumber) '.mat']; % filename of the mat file

N = idx2 - idx1 + 1;

% read acoustic data
numberofWFperfile = 640; % p4457: 640, p4458: 2560, p4459: 2560
WFlength = 8192/2;% p4457: 8192; p4458: 2048; p4459: 2048
acRate = 1000; % in Hz, number of WF per second
acPeriod = 1/acRate;
stackN  = 10; % If you need to stack, put the numbers you want to go into the stack, p4457: 10, p4458: 1

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

energy = zeros(1,acN/stackN); %/10 because of stacking (I stack 10 waveforms to get one sample every 10ms)
hh = 1;
WF = 0; % init waveform
TimeTrig = acPeriod*stackN*ones(1,numberofWFperfile/stackN);
fileID = fopen([acname '/' recname '/' recname '_seis.comb'],'w');
fclose(fileID);
fileID = fopen([acname '/' recname '/' recname '_seis.comb'],'a');
ACmat  = zeros(WFlength+1,numberofWFperfile/stackN);
for ii = 1:1:totalnumberoffiles % normally it should be "1:1:totalnumberoffiles" here (I had only 5000 files on my computer...)
    ii
    ACfilename = [acname '/CF' num2str(recnumber) '/WF_CF' num2str(recnumber) '_' num2str(ii) '.ac'];
    fid = fopen(ACfilename,'r');
    ACdat = fread(fid,[WFlength,numberofWFperfile],'int16');
    if stackN > 1
        for jj = 1:numberofWFperfile/stackN
            ACmat(2:WFlength+1,jj) = sum(ACdat(1:WFlength,stackN*(jj-1)+1:stackN*jj)/stackN,2);
        end
    else
        ACmat(2:WFlength+1,:) = ACdat;
    end
    if ii == 1
        TimeTrig = Time(idx1) - acPeriod*stackN + cumsum(acPeriod*stackN*ones(1,numberofWFperfile/stackN));
        Tlast    = TimeTrig(end);
    else
        TimeTrig = Tlast + cumsum(acPeriod*stackN*ones(1,numberofWFperfile/stackN));
        Tlast    = TimeTrig(end);
    end
    fclose(fid);   
    ACmat(1,:)   = TimeTrig;
    fwrite(fileID,ACmat,'single');
end
% build acoustic time vector (after staking 10, i.e. 1 sample every 10ms)
acTime = 0:stackN*acPeriod:stackN*acPeriod*(acN/stackN-1); % /10 because of stacking 10 waveforms
acTime = acTime + Time(idx1);