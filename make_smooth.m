function [argsmooth,peaks] = make_smooth(arg,Indcpts,nsmooth,peak,peakseq,method)
% A generic function to smoothen any component of velocity step/SHS data smooth by each
% section
% Inputs: 1) arg: Argument variable 2) Indcpts: Points at which departure
% from steady state is imposed 3) nsmooth: Smoothing windows in each
% section 4) peak: Is the data peaked in each section? 0-No, 1-Yes. If data is peaked the
% code breaks the data in each section into two segments--from start to
% peak, from peak to end and smooths each separately.
% 5) peakseq: A vector of the same size as vec_jump. 1 -- find
% maximum, -1 -- find minimum
% 6) method - 'gaussian' or any of the methods under smooth in Matlab

argsmooth_struct = cell(length(Indcpts)-1,1);
% g = gausswin(nsmooth);
% g = g/sum(g);
argsmooth = zeros(size(arg));
peakind   = zeros(length(Indcpts)-1,1);
if length(Indcpts)==1
    Indcpts = [1 Indcpts length(arg)];
end
for i = 2:length(Indcpts)    
    arg_sec = arg(Indcpts(i-1):Indcpts(i)); % extract relevant section
    switch peak
        case 0 % If there are no peaks
            if strcmp(method,'gaussian')
                g = gausswin(nsmooth(i-1));
                g = g/sum(g);
                argsmooth_struct{i-1} = conv2(arg_sec,g,'same');
            elseif strcmp(method,'sgolay')
                argsmooth_struct{i-1} = smooth(arg_sec,nsmooth(i-1),method,min(4,nsmooth(i-1)-1)); % lowess smoothing 
                argsmooth(Indcpts(i-1):Indcpts(i)) = argsmooth_struct{i-1}; % smoothed data being assigned 
            else
                argsmooth_struct{i-1} = smooth(arg_sec,nsmooth(i-1),method); % lowess smoothing 
                argsmooth(Indcpts(i-1):Indcpts(i)) = argsmooth_struct{i-1}; % smoothed data being assigned 
            end
        case 1 % If there are peaks
            if peakseq(i-1) == 1
                [~,peakind(i-1)] = max(arg_sec-min(arg_sec)); % Choose max for high peaks
            elseif peakseq(i-1) == -1
                [~,peakind(i-1)] = min(arg_sec-min(arg_sec)); % Choose min for troughs
            elseif peakseq(i-1) == 0
                peakind(i-1)     = 1; % No peaks, just choose peakind to be 1
            end
            peakind(i-1) = peakind(i-1) + Indcpts(i-1) - 1; % Map peakind to global indices
            widback = floor((peakind(i-1)-Indcpts(i-1))/10); % Width behind peak
            widfrnt  = floor((Indcpts(i)-peakind(i-1))/10); % Width in front of peak
            if widback < 30
                nsmoothback = 1; % For less than 30 points, don't smooth
            else
                nsmoothback = min(nsmooth(i-1),widback); % appropriate smoothing window for initial part
            end
            if widfrnt < 30
                nsmoothfrnt = 1; % For less than 30 points, don't smooth
            else
                nsmoothfrnt = min(nsmooth(i-1),widfrnt); % appropriate smoothing window for ending part
            end       
            argback = arg(Indcpts(i-1):peakind(i-1)-1);
            argfrnt = arg(peakind(i-1):Indcpts(i)); 
            if strcmp(method,'gaussian')
                gback = gausswin(nsmoothback);
                gfrnt  = gausswin(nsmoothfrnt);
                gback = gback/sum(gback);
                gfrnt = gfrnt/sum(gfrnt);                            
                argsmooth_struct{i-1} = [conv2(argback,gback,'same');...
                    conv2(argfrnt,gfrnt,'same')]; % Applying Gaussian smoothing
            elseif strcmp(method,'sgolay')
                argsmooth_struct{i-1} = [smooth(argback,nsmoothback,method,min(4,nsmoothback-1));...
                    smooth(argfrnt,nsmoothfrnt,method,min(4,nsmoothfrnt-1))]; % lowess smoothing 
            else
                argsmooth_struct{i-1} = [smooth(argback,nsmoothback,method);...
                    smooth(argfrnt,nsmoothfrnt,method)]; % lowess smoothing 
            end
            if nsmooth(i-1) == 1
                argsmooth(Indcpts(i-1):Indcpts(i)) = arg_sec; % If nsmooth(i-1) = 1, no smoothing at all
            else
                argsmooth(Indcpts(i-1):Indcpts(i)) = argsmooth_struct{i-1}; % smoothed data being assigned
            end
        otherwise
    end
            
end
peaks = peakind;
