function [peak_pos,corr_mat,lag_mat,template,peak_corr,peak_ind] = analyze_WFS(seism_mat,runname,recname,n_start,phase,ncyc,dt,fulstack)
%% This is a code to perform cross-correlation analysis on waveforms
% Inputs: 
%        seism_mat: Matrix of waveform, each waveform is one column, as
%                   many columns as there are waveforms
%        runname  : Exp. name, character array: e.g. 'p4457' etc.
%        recname  : Control File number, double: 1, 2, 3 etc.
%        n_start  : Which sample of the original WF is stored in the first
%        row of seism_mat
%        phase    : string, which phase: 'p' or 's'.
%        ncyc     : number of cycles
%        dt       : Sampling interval for the WF
%        fulstack : User defined stack, optional
% Outputs:
%        peak_pos : col1: index for peaks, col2: index for minima, 
%                   col3: peak2peak amplitude, col4: Arrival time for phase
%        corr_mat : matrix of correlation vectors
%        lag_mat  : matrix of lag vectors for each correlation vector

pick  = load(['./' runname '_acoustics/' runname '_' phase '-pick.txt']);
nchannel = size(seism_mat,2); % Number of seismograms
peak_pos  = zeros(nchannel,4); % Matrix for peak_positions, col1: peak-P, col2: min-P, col3: peak-S, col4: min-S
rown_pick = find(pick(:,1)==recname); % Which row of picks file
npick = pick(rown_pick,2)-n_start+1; %% First index of phase
nwind = pick(rown_pick,3)-pick(rown_pick,2)-1; %% Total width
nwindow = floor(nwind/(4*ncyc)); % Window length for smoothing: 1 cycle
g = gausswin(nwindow,4); % Ensure that smoothing window goes to zero within nwindow samples by choosing a large alpha (here 4)
g = g/sum(g);
%seism_mat = seism_mat*0.1;
seism_mat_filt  = conv2(seism_mat,g,'same');
%% Cross-correlation
%Now compute cross-correlations
    maxlags     = nwind;
    corr_mat   = zeros(nchannel,2*maxlags+1);
    lag_mat    = zeros(nchannel,2*maxlags+1);
    corrmaxlag = zeros(nchannel,1);
    peak_corr  = zeros(nchannel,1);
    peak_ind   = zeros(nchannel,1);
    npickin   = npick;
    if nargin<8
        [template, ~, npick] = make_template(seism_mat_filt(:,3:end),npick,nwind,[],1,0,2,ncyc);
    else
         template  = fulstack(npick:npick+nwind);
    end
    if npickin~=npick
        disp(['NPICK adjusted by ' num2str(abs(npick-npickin)) ' indices']);
    end
    section   = seism_mat_filt(npick-100:npick+nwind,:); % Choose the section to be filtered
    for i = 1 : nchannel
        xs                       = section(:,i);
        ys                       = [template(:);zeros(100,1)]; % zero pad the end
        [corrs,lags]             = xcorr(xs,ys,maxlags,'coeff');
        [peak,pk_ind]            = max(corrs);
        corr_mat(i,:)            = corrs;
        lag_mat(i,:)             = lags;
        corrmaxlag(i)            = lags(pk_ind);
        peak_corr(i)             = peak;
        peak_ind(i)              = pk_ind;    
    end           
corrmaxind              = corrmaxlag + npick - 1;
p_s_delay               = corrmaxind;
%[pks,locs]             = findpeaks(abs(seism_mat_filt(i,corrmaxindp(i):corrmaxindp(i)+ceil(nwind_p))),'sortstr','descend');
[pks,locs]              = max(seism_mat_filt(corrmaxind:corrmaxind+nwind,:),[],1);
[mins,lmins]            = min(seism_mat_filt(corrmaxind:corrmaxind+nwind,:),[],1);
peak_pos(:,1:2)         = [locs' lmins'];
Ain                     = max(template)-min(template); % Change on template
p_peak_mat              = (pks'-mins')/Ain; % Normalised change in transmitted amplitude, for our frequencies, proportional to change in interface stiffness
%p_peak_mat(i)          = seism_mat(i,npick_p+locs-1)-seism_mat(i,npick_p+lmins-1); %srec4
travel_time_p           = parabolic_extrapolation(npick,n_start,corr_mat,lag_mat,dt);
peak_pos(:,3:4)         = [p_peak_mat travel_time_p];
%waveform                = hann(nwind,'periodic').*(seism_mat_filt(corrmaxind:corrmaxind+nwind,:) - mins)/(pks-mins);% nwind point Hanning window applied to waveform
