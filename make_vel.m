function vel = make_vel(disp,time,Indcpts,nsmooth)
vel_struct = cell(length(Indcpts),1);
disp_struct = cell(length(Indcpts),1);
% g = gausswin(nsmooth);
% g = g/sum(g);
vel = zeros(size(disp));

for i = 2:length(Indcpts)
    disp_struct{i} = disp(Indcpts(i-1):Indcpts(i));
    wind_len       = length(disp(Indcpts(i-1):Indcpts(i)));
    disp_struct{i} = smooth(disp_struct{i},nsmooth(i-1),'lowess');
    if (wind_len<=10000)
        center_diff    = diag(zeros(wind_len,1))+diag(ones(wind_len-1,1),1)+diag(-ones(wind_len-1,1),-1);
        center_diff(1,1) = -1;
        center_diff(wind_len,wind_len) = 1;
        vel_struct{i}    = 0.5*center_diff*disp_struct{i}/mean(diff(time(Indcpts(i-1):Indcpts(i))));
        vel_struct{i} = smooth(vel_struct{i},nsmooth(i-1),'lowess');
        vel(Indcpts(i-1):Indcpts(i)) = vel_struct{i};
    else
        parfor j = Indcpts(i-1) + 3 : Indcpts(i) - 3
            first_ord  = 39*(disp(j+1) - disp(j-1));
            second_ord = 12*(disp(j+2) - disp(j-2));
            third_ord  = -5*(disp(j+3) - disp(j-3));
            denom      = 96*(time(j+1) - time(j-1));
            vel(j) = (first_ord+second_ord+third_ord)/denom;
        end
    end
    vel(Indcpts(i-1):Indcpts(i)) = smooth(vel(Indcpts(i-1):Indcpts(i)),nsmooth(i-1),'lowess');        
end

    