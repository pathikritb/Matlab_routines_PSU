function [template, fullstack, n_pick] = make_template(seism_mat,n_pick,n_wind,vectemp,stack_opt,plot_opt,stackdim,adj_opt)

%% Inputs: seism_mat: Matrix of unfiltered seismograms
%          n_pick   : Pick for the phase, expecting an eyeball phase pick
%                       for this
%          n_wind   : How many cycles in the template, in terms of number
%                       of samples from the pick.
%          vectemp  : vector of indices for seismograms to be used for the stack, only useful for a non-global stack
%          stack_opt: Options for stacking for the template: 1: sum all seismograms, then pick a phase from that, 
%                       2: interspersed samples of nstacp seismograms each from different portions of the data.
%          plot_opt : Plot the stacks? 1-Yes, 0-No, plotting not
%                       available for stack_opt = 1
%          stackdim : Which dimension to stack along? Optional argument.
%                       1-along rows, 2-along columns
%          adj_opt  : Optional input, Do we want to adjust n_pick based on
%                       the closest zero crossing before the user defined n-pick, user can check
%                       if this is appropriate by plotting template on
%                       fullstack. 0-No adjustment necessary, also default.
%                       n-Yes, adjust. n is the number of cycles in n_wind
%% Outputs: 
%           template: The particular phase template
%           fullstack : The stacked, full seismogram
%           n_pick  : optional output, only important if adj_opt > 0
    if nargin>6 && (stackdim==2)
        seism_mat = seism_mat';
    end
    len          = size(seism_mat);
    n_seism      = len(2);
    ntemp        = length(vectemp);
    seism_master = zeros(ntemp-1,n_seism); 
    t_offset     = zeros(ntemp-1,1);
    if stack_opt == 1
        fullstack = sum(seism_mat(1:len(1),:),1)/len(1); % The stacked full seismogram
        %template: Between about n_pick-th and (n_pick+n_wind)-th points on each master waveform
    elseif stack_opt == 2
        for i = 1 : ntemp-1
            nstack             = vectemp(i+1) - vectemp(i);
            seism_master(i,:) = sum(seism_mat(vectemp(i):vectemp(i+1),:))/nstack;
            if plot_opt == 1
                if i == 1
                    t_offset(i)     = max(seism_master(i,:))*0;
                elseif i > 1
                    t_offset(i)     = max(seism_master(i,:))*1.5 + t_offset(i-1);
                end
            end
        end
        fullstack = sum(seism_master(:,:),1)/ntemp; % The stacked full seismogram
        %template: Between about n_pick-th and (n_pick+n_wind)-th points on each master waveform
    end
    
    if nargin<=7
            % adj_opt does not exist
            template  = fullstack(n_pick:n_pick+n_wind);
        elseif nargin>7 && adj_opt==0
            % adj_opt does exist, but no adjustment
            template  = fullstack(n_pick:n_pick+n_wind);
        elseif nargin>7 && adj_opt~=0
            % adj_opt does exist, adjust
            ncyc = adj_opt;
            widsearch = round(round(n_wind/ncyc)/4); % 1/4th of period is our search window;
            for i = n_pick-widsearch:n_pick+widsearch
                if (fullstack(i)*fullstack(i+1)<0) && (n_pick~=i+1)
                    n_pick = i+1;
                end
            end
            template  = fullstack(n_pick:n_pick+n_wind);
     end
        
   
    if plot_opt == 1
        figure
        plotx = repmat(1:n_seism,ntemp-1,1);
        ploty = seism_master+repmat(t_offset,1,n_seism);
        plotyn = seism_mat(vectemp(1:end-1),:)+repmat(t_offset,1,n_seism);
        plot(plotx(:,1:end)',ploty(:,1:end)','k-','linewidth',2)
        hold on
        plot(plotx(:,1:end)',plotyn(:,1:end)','b--')
        plot(plotx(:,n_pick:n_pick+n_wind)',ploty(:,n_pick:n_pick+n_wind)','r-','linewidth',2)
    end