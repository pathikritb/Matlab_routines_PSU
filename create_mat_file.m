% Access memmapfile and then analyze waveform data.
% Expects Wave Form memmapfiles created by sync_n_save or similarly
% formatted binary files. Offers to save as mat file if file does not exist
% Uses Parallel Computing Toolbox to run parfor loops.
function outname = create_mat_file(runname,recname,idxfile_path,memfile,out)
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
% (4) memfile = $path for binary file with acoustic data, expects format:- 
%       N blocks each of: 1st row: Times, 2nd row: Biax Column, 3rd to (WFlength+2)th row : Waveform Vector
%       ith column of Block(j) = [Time_of_WF(i);Biax_Column(i);WF vector(i)]
%       This format is naturally created by sync_n_sav
%       N is the number of trigger files
% (5) out = $path for output file, default : /your_path/runname_acoustics/recname/
% Output:
% outname: $path for output file
%%---------------------------------------------------------------------------------------------------%%
% Filename and argcheck
acname = ['./' runname '_acoustics']; % setting path to acoustics file
if (nargin<3)
    idxfile_path = ['./' runname '_acoustics/index_CF.txt'];
    memfile      = [acname '/' recname '/' recname '_seis_git.comb'];
    out          = ['./' runname '_acoustics/' recname '_withstack_git.mat'];
elseif (nargin<4) 
    memfile      = [acname '/' recname '/' recname '_seis_git.comb'];
    out          = ['./' runname '_acoustics/' recname '_withstack_git.mat'];
elseif (nargin<5)
    out          = ['./' runname '_acoustics/' recname '_withstack_git.mat'];
end
if isempty(idxfile_path)
        idxfile_path = ['./' runname '_acoustics/index_CF.txt'];
    elseif isempty(memfile)
        memfile      = [acname '/' recname '/' recname '_seis_git.comb'];
    elseif isempty(out)
        out      = [acname '/' recname '/' recname '_seis_git.comb'];
end
% Processing begins
Fs      = 25e6; % Sampling of waveforms, not pinging rate
dt      = 1/Fs; % Sampling period
recnumber = str2double(recname(3)); % number used to save data (e.g. 2 for control file 2)
fID = fopen(idxfile_path);
C = textscan(fID,'%s %u32 %u32 %u32','HeaderLines',1);
fclose(fID);
idx1    = C{1,2}(recnumber); % Set first  biax_col index
idx2    = C{1,3}(recnumber);  % Set last  biax_col index
ind = idx1:idx2;
nworkers = 1000; % Number of workers on system for use in parallel, remember nworkers >> nthreads for best results
if strcmp(runname,'p4457')
    WFlength = 8192/2;% p4457: 8192, but only half saved
elseif strcmp(runname,'p4458') || strcmp(runname,'p4459')
    WFlength = 2048;% p4458: 20480
else
   WFlength = 2048;% p4458: 2048
end
%% Sync of acoustic data
% Does the requisite mat fle exist
if exist(out, 'file') ~= 2
    replace = 1;
elseif exist(out,'file') == 2
        replace = input([out 'already exists. Rplace file? 1-Yes,0-No']);
        if (replace~=1 || replace~=0)
            error('Wrong option, use either 0 or 1')
        end
end

if (replace == 1)
    seism_st = 800; % What element of the recorded waveform to start on
    seiswid  = 1200; % How many wform elements to use starting at seism_st
    %% Sync begins
    bytes = getfield(dir(memfile), 'bytes'); % Total size of memfile
    totNWF = bytes/(4*(WFlength+2));% Then number of waveforms assuming format as specfied in help
%% Stack sections
    if isempty(gcp('nocreate'))
        parpool; % create parallel pool if not already open
    end
    nstack = round(totNWF/nworkers); % size of each parallel block
    if (nstack*nworkers==totNWF)
        indstack = 1:nstack:totNWF+nstack; % Find indices for each parallel segment
    else
        indstack = 1:nstack-1:totNWF+nstack-1; % Corrrect if totNWF is not exactly divisible by nworkers
    end
    indstack(end) = totNWF;
    lenstack = zeros(nworkers,1);
    for i = 1:nworkers
        lenstack(i) = indstack(i+1) - indstack(i) + 1; % Exact number of waveforms in each stack
    end
    lenstack(nworkers) = totNWF - sum(lenstack(1:nworkers-1));
    % Create cell files for each variable
    time_ac = cell(nworkers,1); % Times corr. to each amplitude
    biax_ac = cell(nworkers,1); % Biax Col corr. to each amplitude
    max_WFP  = cell(nworkers,1); % Max P amplitude
    max_WFS  = cell(nworkers,1); % Max S amplitude
    P_delay  = cell(nworkers,1); % P-delay
    S_delay  = cell(nworkers,1); % S-delay
    StackP   = cell(nworkers,1); % Stack used for P
    StackS   = cell(nworkers,1); % Stack used for S
    lagcell  = cell(nworkers,1); % All cross-correlation lags
    corrcell = cell(nworkers,1); % All cross-correlation values
    pk_corrP = cell(nworkers,1); % Maximum correlations for P
    pk_indP  = cell(nworkers,1); % Index for peak P correlation
    pk_corrS = cell(nworkers,1); % Maximum correlations for P
    pk_indS  = cell(nworkers,1); % Index for peak S correlation
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
        [peak_posP,~,~,stackP,peak_corrP,peak_indP] = analyze_WFS(locdat,runname,2,seism_st,'p',3,dt,fulstack); % 3 P-cycles
        [peak_posS,~,~,stackS,peak_corrS,peak_indS] = analyze_WFS(locdat,runname,str2double(recname(3)),seism_st,'s',3,dt,fulstack); % Only 1 S-cycle due to clipping of amplitudes
        %lagcell{i}  = lag_matS;
        %corrcell{i} = corr_matS;
        time_ac{i}  = double(fseism_comb.Data.seism(1,:));
        biax_ac{i}  = uint32(fseism_comb.Data.seism(2,:));
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
    %% Convert cells to vectors
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
    biax_ac_vec = biax_ac_vec - ind(1) + 1;
    save(out,...
    'time_ac_vec','biax_ac_vec','max_WFS_vec','max_WFP_vec','P_delay_vec','S_delay_vec' ...
        ,'StackP','StackS','pk_corrP_vec','pk_indP_vec','pk_corrS_vec','pk_indS_vec'); % Save in mat file
end
outname = out;