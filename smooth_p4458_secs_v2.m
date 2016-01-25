% smooth sections of p4457 data
%% first load cpts_data; IndCpts 9:22 (refer the cpts file) run1, IndCpts 24:38 run2
clear
clc
runname = 'p4458';
sel  = 'CF3';
Fs      = 25e6;
dt      = 1/Fs;
[data,outname] = ReadBinBiax(runname); % Path: Read directly from binary file
acoustic = 1;
siz      = size(data);
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
ind = eval(sel);
cpts = load(['./' runname '/data/' runname '_cpts.txt']);
nworkers = 1000; % Number of workers on system
select = cpts(cpts(:,2)>=ind(1)&cpts(:,2)<=ind(end),1);
vel    = cpts(:,4);
cpts = cpts(select,:);
cpts(:,2) = cpts(:,2) - ind(1) + 1;
acname = ['./' runname '_acoustics'];
%WFlength = 8192;
%% Central dcdt has an issue, reverses direction at 595238 for p4457
% Dealt with later
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
%% Sync of acoustic data
ac_mat_name = ['./' runname '_acoustics/' sel '_withstack_v2.mat'];
if acoustic == 1 && exist(ac_mat_name, 'file') ~= 2
% CF2 (indexes for Control File 2)
    idx1 = ind(1); % 208309 for CF1, 1311457 for CF2, 2822024 for CF3; index corresponding to first trigger of the sequence on figure 1
    idx2 = ind(end); % 1233694 for CF1, 2764167 for CF2, 3804488 for CF3; index corresponding to last trigger of the sequence on figure 1
    idxmin_peak = idxmin_peak_vec(str2double(sel(3)));% 1065560 for CF1, 1674979 for CF2, 3260187 for CF3; index corresponding to a large trigger in the middle of the sequence
%     idx1 = 606891; % Jacques chose 606891 for CF2, index corresponding to first trigger of the sequence on figure 1
%     idx2 = 1527209; % Jacques chose 1527209 for CF2, index corresponding to last trigger of the sequence on figure 1
%     idxmin_peak = 673653; % index corresponding to a large trigger in the middle of the sequence
%     recnumber = 2; % number used to save data (e.g. 2 for control file 2)
    N = idx2 - idx1 + 1;    
    if strcmp(runname,'p4457')
        WFlength = 8192/2;% p4457: 8192, but only half saved
        numberofWFperfile = 640; % p4457: 640, p4458: 2560
    else
        WFlength = 2048;% p4458: 2048
        numberofWFperfile = 2560; % p4457: 640, p4458: 2560
    end
    seism_st = 800; % What element of the recorded waveform to start on
    seiswid  = 1200; % How many wform elements to use starting at 870
    acRate = 1000; % in Hz, number of WF per second
    acPeriod = 1/acRate;
    %% Sync begins
    memfile = [acname '/' sel '/' sel '_seis_v2.comb'];
    bytes = getfield(dir(memfile), 'bytes');
    totNWF = bytes/(4*(WFlength+2));
%% Stack sections
    if isempty(gcp('nocreate'))
        parpool(10);
    end
    nstack = round(totNWF/nworkers);
    if (nstack*nworkers==totNWF)
        indstack = 1:nstack:totNWF+nstack;
    else
        indstack = 1:nstack-1:totNWF+nstack-1;
    end
    indstack(end) = totNWF;
    lenstack = zeros(nworkers,1);
    for i = 1:nworkers
        lenstack(i) = indstack(i+1) - indstack(i) + 1;
    end
    lenstack(nworkers) = totNWF - sum(lenstack(1:nworkers-1));
    time_ac = cell(nworkers,1);
    biax_ac = cell(nworkers,1);
    max_WFP  = cell(nworkers,1);
    max_WFS  = cell(nworkers,1);
    P_delay  = cell(nworkers,1);
    S_delay  = cell(nworkers,1);
    StackP   = cell(nworkers,1);
    StackS   = cell(nworkers,1);
    lagcell  = cell(nworkers,1);
    corrcell = cell(nworkers,1);
    pk_corrP = cell(nworkers,1);
    pk_indP  = cell(nworkers,1);
    pk_corrS = cell(nworkers,1);
    pk_indS  = cell(nworkers,1);
    locdat  = [];
    %kk      = zeros(nworkers,1);
    %% Make a global template for cross-correlation
    fulstack = zeros(seiswid,nworkers);
    parfor i = 1:nworkers
        %kk(i) = 1;
        disp(['Stacking ' num2str(i) 'th worker'])
        if i == 1
            Offset = 0;
        else
            Offset = (WFlength+2)*indstack(i)*4; % Number of bytes to skip in SP for time
        end
        fseism_comb    = memmapfile(memfile,...
            'Format',{'single',[WFlength+2 lenstack(i)],'seism'},...
            'Repeat',1,...
            'Offset',Offset);% CF2_seis.comb only has half of the waveforms, also every 10 stacked
        locdat     = double(fseism_comb.Data.seism(seism_st:seism_st+seiswid-1,:));
        fulstack(:,i)   = fulstack(:,i) + sum(locdat(:,3:end),2)/(lenstack(i)-2);
    end
    fulstack = sum(fulstack,2)/nworkers;
    %% Parallel access of segments of data map
    parfor i = 1:nworkers
        %kk(i) = 1;
        i
        if i == 1
            Offset = 0;
        else
            Offset = (WFlength+2)*indstack(i)*4; % Number of bytes to skip in SP for time
        end
        fseism_comb    = memmapfile(memfile,...
            'Format',{'single',[WFlength+2 lenstack(i)],'seism'},...
            'Repeat',1,...
            'Offset',Offset);% CF2_seis.comb only has half of the waveforms, also every 10 stacked
        locdat     = double(fseism_comb.Data.seism(seism_st:seism_st+seiswid-1,:));
        [peak_posP,corr_matP,lag_matP,stackP,peak_corrP,peak_indP] = analyze_WFS(locdat,runname,2,seism_st,'p',3,dt,fulstack); % 3 P-cycles
        [peak_posS,corr_matS,lag_matS,stackS,peak_corrS,peak_indS] = analyze_WFS(locdat,runname,str2double(sel(3)),seism_st,'s',3,dt,fulstack); % Only 1 S-cycle due to clipping of amplitudes
        %lagcell{i}  = lag_matS;
        %corrcell{i} = corr_matS;
        time_ac{i}  = round(double(fseism_comb.Data.seism(1,:))*1000)/1000;
        biax_ac{i}  = double(fseism_comb.Data.seism(2,:));
        max_WFP{i}  = peak_posP(:,3);
        max_WFS{i}  = peak_posS(:,3);
        P_delay{i}  = peak_posP(:,4);
        S_delay{i}  = peak_posS(:,4);
        StackP{i}   = stackP;
        StackS{i}   = stackS;
        pk_corrP{i} = peak_corrP;
        pk_indP{i}  = peak_indP;
        pk_corrS{i} = peak_corrS;
        pk_indS{i}  = peak_indS;
        
        %disp([num2str(sum(kk)/nworkers*100) '% done----------------'])
    end
    %%
    time_ac_vec = horzcat(time_ac{:})';
    biax_ac_vec = horzcat(biax_ac{:})';
    max_WFP_vec = vertcat(max_WFP{:})';
    max_WFS_vec = vertcat(max_WFS{:})';
    P_delay_vec = vertcat(P_delay{:})';
    S_delay_vec = vertcat(S_delay{:})';
    pk_corrP_vec = vertcat(pk_corrP{:})';
    pk_indP_vec = vertcat(pk_indP{:})';
    pk_corrS_vec = vertcat(pk_corrS{:})';
    pk_indS_vec = vertcat(pk_indS{:})';
    t_sub_ac = Time(idx1:idx2);
    biax_ac_vec = biax_ac_vec - ind(1) + 1;
    save(ac_mat_name,...
    'time_ac_vec','biax_ac_vec','max_WFS_vec','max_WFP_vec','P_delay_vec','S_delay_vec' ...
        ,'StackP','StackS','pk_corrP_vec','pk_indP_vec','pk_corrS_vec','pk_indS_vec');
elseif acoustic == 1 && exist(ac_mat_name,'file') == 2
        load(ac_mat_name);
end

%% Smoothing etc.
%time_ac_vec = horzcat(time_ac{:})';
%max_WF_vec  = horzcat(max_WF{:})';
%t_sub_ac = Time(idx1:idx2);
if ind(end) ~= length(data)
    Indcpts = [1;cpts(:,2); length(mu)];
    if select(1)~=1
        vel     = [vel(select(1)-1);vel(select)];
    else
        vel     = [NaN;vel(select)];
    end
else
    Indcpts = [1;cpts(:,2)];
    if select(1)~=1
        vel     = [vel(select(1)-1);vel(select(1:end-1))];
    else
        vel     = [NaN;vel(select(1:end-1))];
    end
end
nsmooth = round(diff(Indcpts)/50);
%nsmooth = [1;nsmooth];
peakseq = zeros(length(vel)+1,1);
for i = 1:length(vel)-1
    if vel(i+1)>vel(i)
        peakseq(i+1) = 1;
    elseif vel(i+1)<vel(i) && vel(i+1)~=0
        peakseq(i+1) = -1;
    end
end         
% Use loess for non-peaked data, Use moving for peaked data
p_closure = polyfit(slip_lp,closure,1);% Trend in closure
closure_tr = p_closure(1)*slip_lp + p_closure(2);% Compute trend in closure
closure_det = closure-closure_tr;% Detrend closure
slip_lp_s = make_smooth(slip_lp,Indcpts,nsmooth,0,peakseq,'loess');
surf_slip_s = make_smooth(surf_slip,Indcpts,nsmooth,0,peakseq,'loess');
mu_s      = make_smooth(mu,Indcpts,nsmooth,1,peakseq,'moving');
cents     = make_smooth(cent_blk,Indcpts,nsmooth,0,peakseq,'loess');
sigmas    = make_smooth(sigma,Indcpts,nsmooth,0,peakseq,'loess');
closures  = make_smooth(-closure_det,Indcpts,nsmooth,0,peakseq,'loess');
closures_corr = -(closures*1e-6 - (0.163/15e9)*((sigmas-4.2) - 4*0.15*sigmas.*(mu_s-0.76))*1e6)*1e6;