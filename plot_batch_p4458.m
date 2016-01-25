% Plot useful stuff from the data from p4457/p4458 type datsets
stiff_eval = 1;
ind_reload = find(vel==0)+1;
if strcmp(runname,'p4457') && strcmp(sel,'CF2')
    ind_reload(1) = [];
end
Vrs = vel(ind_reload(end-1));
ind_reload(ind_reload>length(vel))=[];
ind_reload(isnan(vel(ind_reload)))=[];
ind_holds = ind_reload - 1;
ind_reload = Indcpts(ind_reload);
ind_holds  = Indcpts(ind_holds);
ind_1 = Indcpts(find(Indcpts==ind_holds(1))-1);
ind_end = Indcpts(find(Indcpts==ind_reload(end))+1);
ind_pl = [ind_holds';ind_reload'];
ind_pl = [ind_1;ind_pl(:);ind_end];
plind  = ind_1:ind_end-1;
if strcmp(runname,'p4458') && strcmp(sel,'CF3')
    [time_ac_loc,itu] = unique(time_ac_vec+0);
elseif strcmp(runname,'p4457') && strcmp(sel,'CF1')
    [time_ac_loc,itu] = unique(time_ac_vec-0.10);
else
    [time_ac_loc,itu] = unique(time_ac_vec);
end
biax_ac_loc = biax_ac_vec(itu);    
max_WFP_loc = max_WFP_vec(itu);% Multiply by Vp/Vs ratio to make sure scaling is perfect, % ge change in contact stiffness
max_WFS_loc = max_WFS_vec(itu); % ge change in contact stiffness
max_WF_loc = max_WFS_loc;
P_delay_loc = P_delay_vec(itu);
S_delay_loc = S_delay_vec(itu);
lp_ac   = slip_lp_s(biax_ac_loc);
stiffness = 0.008;
slipf     = slip_lp_s - mu_s/stiffness;
if stiff_eval
    nwidth    = 15;
    npwidth   = 500;
    c = {'c','r','g','c','r','g','b','k'};
    stiff = zeros(length(ind_reload),1);
    off_ind = 200;
    musec = zeros(off_ind+npwidth,length(ind_reload));
    lpsec = zeros(off_ind+npwidth,length(ind_reload));
    lpnssec = zeros(off_ind+npwidth,length(ind_reload));
    tsec  = zeros(off_ind+npwidth,length(ind_reload));
    figure
    for i = 3:length(ind_reload)
        indstiff = ind_reload(i)-off_ind + 1:ind_reload(i)+npwidth;
        lp_sec = slip_lp_s(indstiff)/1e3;        
        lp_ns  = cent_blk(indstiff)/1e3;
        mu_sec = mu_s(indstiff);
        slip_sec = slipf(indstiff)/1e3;
        t_sec  = time(indstiff); 
        clsec = closures(indstiff);
        Vf_sec   = make_vel(slip_sec*1e3,t_sec,[1 off_ind length(indstiff)],[10 10]);       
        indac  = time_ac_loc>=t_sec(1) & time_ac_loc<=t_sec(end);
        tac    = time_ac_loc(indac); 
        lp_ac_loc = lp_ac(indac)/1e3;
        slip_ac_loc = (lp_ac_loc*1e3-mu_s(indac)/stiffness)/1e3;
        WF_loc    = max_WF_loc(indac);
        musec(:,i) = mu_sec;
        lpsec(:,i) = lp_sec;
        lpnssec(:,i) = lp_ns;
        tsec(:,i) = t_sec;
        subplot(1,length(ind_reload)-2,i-2)
        %[AX,H1,H2] = plotyy(t_sec,mu_sec,t_sec,slip_sec-slip_sec(off_ind+1));
        %[AX,H1,H2] = plotyy(t_sec,mu_sec,t_sec(Vf_sec>0),log10(Vf_sec(Vf_sec>0)));
        [AX,H1,H2] = plotyy(t_sec,mu_sec,t_sec,log10(Vf_sec));
        hold(AX(1),'on');
        hold(AX(2),'on');
        plot(AX(1),tac,WF_loc,'b--','linewidth',2);
        if strcmp(runname,'p4458')
            surf_slip_sec = surf_slip(indstiff)/1e3;
            %plot(AX(2),t_sec,surf_slip_sec-surf_slip_sec(off_ind+1),'r--','linewidth',2);
        end
        set(H1,'color','b','linewidth',2);
        set(H2,'linewidth',2)
        ylim = linspace(min(WF_loc)-0.05*min(WF_loc),max(WF_loc)+0.05*max(WF_loc),5);
        ylim = linspace(min(min(mu_sec),round(ylim(1)*10)/10),round(ylim(end)*10)/10,5);
        xlim = linspace(min(t_sec),max(t_sec),5);
        set(AX(1),'ycolor','b','ylim',[ylim(1) ylim(end)],'ytick',ylim,'xlim',[xlim(1) xlim(end)]);
        set(AX(2),'xlim',[xlim(1) xlim(end)]);
        p = polyfit(lp_sec(1:nwidth)-lp_sec(off_ind+1),mu_sec(1:nwidth)-mu_sec(off_ind+1),1);
        lfit = polyval(p,lp_sec(1:nwidth)-lp_sec(off_ind+1));
        %plot(lp_sec(1:nwidth)-lp_sec(1),lfit-lfit(1),'-','color',c{i},'linewidth',3)
        stiff(i) = p(1);       
        xlabel(AX(1),'Time [s]')
        ylabel(AX(2),'log_{10}(V_f)');
        legend(AX(1),'\Delta\mu','S-Amp.')
        set(findall(gcf,'-property','FontSize'),'FontSize',16) 
    end
end
%%
figure
t_prox = time;
wts = 1./diff(time(ind_pl));
%aa = [0.681559 0.781988 0.746114;0.675393 0.793715 0.7462;0.668662 max(mu_s) 0.7455];
mu_av = zeros(length(ind_holds),1);
for i = 1:length(ind_holds)
    mu_av(i) = sum(mu_s(ind_holds(i)-50:ind_holds(i)))/50;
end
%mu_av = aa(:,3);
mu_peaks = zeros(length(ind_holds),3);
%mu_peaks(:,1) = aa(:,1) - mu_av;
mu_peaks(:,1) = mu_s(ind_reload)-mu_av;
k = 0;
for i=1:length(ind_pl)-1
    t_prox(ind_pl(i):ind_pl(i+1)-1) = (t_prox(ind_pl(i):ind_pl(i+1)-1)-t_prox(ind_pl(i)))*wts(i)+(i-1);
    if ismember(ind_pl(i),ind_reload)
        k = k + 1;
        mu_peaks(k,2) = max(mu_s(ind_pl(i):ind_pl(i+1)-1));%;aa(k,2);
        mu_peaks(k,3) = round(time(ind_reload(k))-time(ind_holds(k)));
    end
end
tplot = t_prox(plind);
plot(tplot,mu(plind));
hold on;
plot(tplot,mu_s(plind),'linewidth',2);
title(['V_{s/r} = ' num2str(Vrs) ' \mums^{-1}'])
for i = 1:length(ind_holds)
    text(t_prox(ind_holds(i))+0.15,0.065,[num2str(mu_peaks(i,3)) 's hold']);
end
xlabel('Scaled time');ylabel('Friction Coeff.')    
set(gca,'ylim',[0.6 0.9],'xlim',[min(tplot) max(tplot)])
figure
semilogx(mu_peaks(:,3)*Vrs,mu_peaks(:,1),'s');hold on,semilogx(mu_peaks(:,3)*Vrs,mu_peaks(:,2)-mu_av,'o');
title(['V_{s/r} = ' num2str(Vrs) ' \mums^{-1}'])
ylabel('\Delta\mu');xlabel('t_{hold}*V_{s/r}');
timacloc = cell(length(ind_holds),1);
max_WFloc = cell(length(ind_holds),1);
figure
for i=1:length(ind_holds)
    subplot(2,3,i)
    pl_loc = ind_holds(i)-round((ind_reload(i)-ind_holds(i))/100):ind_holds(i)+round((ind_reload(i)-ind_holds(i))+1);
    t_loc = time(pl_loc);
    timacloc{i} = time_ac_loc(time_ac_loc>=t_loc(1) & time_ac_loc<=t_loc(end));
    max_WFloc{i} = max_WF_loc(time_ac_loc>=t_loc(1) & time_ac_loc<=t_loc(end));
    [AX,H1,H2] = plotyy(t_loc,mu_s(pl_loc),timacloc{i},max_WFloc{i});
    xlabel(AX(1),'Time [s]')
    set(AX(1),'xlim',[min(t_loc) max(t_loc)],'xscale','log');
    set(AX(2),'xlim',[min(t_loc) max(t_loc)],'xscale','log');
    title([num2str(round(time(ind_reload(i))-time(ind_holds(i)))) 's hold'])
    if (mod(i,3)==1)
        ylabel(AX(1),'Friction coeff.');
    elseif (mod(i,3)==0||i==length(ind_holds))
        ylabel(AX(2),'Transmissivity');
    end
end
xlabel(AX(1),'Time [s]')
figure
ind_pl_st  = zeros(length(ind_holds),1);
max_WFPloc = cell(length(ind_holds),1);
max_WFSloc = cell(length(ind_holds),1);
Pdelayloc  = cell(length(ind_holds),1);
Sdelayloc  = cell(length(ind_holds),1);
for i=1:length(ind_holds)
    subplot(2,3,i)
    ind_pl_st(i) = 0;%round((ind_reload(i)-ind_holds(i))/100);
    pl_loc = ind_holds(i)-ind_pl_st:ind_holds(i)+round((ind_reload(i)-ind_holds(i))+1);
    t_loc = time(pl_loc);
    timacloc{i} = time_ac_loc(time_ac_loc>=t_loc(1) & time_ac_loc<=t_loc(end));
    max_WFloc{i} = max_WF_loc(time_ac_loc>=t_loc(1) & time_ac_loc<=t_loc(end));
    max_WFPloc{i} = max_WFP_loc(time_ac_loc>=t_loc(1) & time_ac_loc<=t_loc(end));
    max_WFSloc{i} = max_WFS_loc(time_ac_loc>=t_loc(1) & time_ac_loc<=t_loc(end));
    Pdelayloc{i} = P_delay_loc(time_ac_loc>=t_loc(1) & time_ac_loc<=t_loc(end));
    Sdelayloc{i} = S_delay_loc(time_ac_loc>=t_loc(1) & time_ac_loc<=t_loc(end));
    [AX,H1,H2] = plotyy(t_loc,slip_lp(pl_loc),timacloc{i},max_WFloc{i});
    xlabel(AX(1),'Time [s]')
    set(AX(1),'xlim',[min(t_loc) max(t_loc)],'xscale','log','NextPlot','add');
    semilogx(t_loc,slip_lp_s(pl_loc),'g-','linewidth',2)
    set(AX(2),'xlim',[min(t_loc) max(t_loc)],'xscale','log');
    title([num2str(round(time(ind_reload(i))-time(ind_holds(i)))) 's hold'])
    if (mod(i,3)==1)
        ylabel(AX(1),'Load Pt.');
    elseif (mod(i,3)==0||i==length(ind_holds))
        ylabel(AX(2),'Transmissivity');
    end
end
xlabel(AX(1),'Time [s]')
figure
legendstr = {};
cseq      = {'r','g','b','c','m','k','y'};
Lsol      = [];
Thickness_block = 0.163; % Thickness of assembly
P_vel_bulk      = 5500; % Bulk P-velocity in m/s at 4 MPa
S_vel_bulk      = 3200; % Bulk S-velocity in m/s at 4 MPa
BulkPdelay      = Thickness_block/P_vel_bulk*1e6; % P-Delay due to travel in bulk in microsec
BulkSdelay      = Thickness_block/S_vel_bulk*1e6; % S-Delay due to travel in bulk in microsec
for i = length(ind_holds):-1:1
    %scale = max_WFloc{i}(end) - max_WFloc{i}(ind_pl_st(i)+1);
    scaleP = max_WFPloc{i}(ind_pl_st(i)+1);%max_WFPloc{length(ind_holds)}(end) - max_WFPloc{length(ind_holds)}(ind_pl_st(i)+1);
    scaleS = max_WFSloc{i}(ind_pl_st(i)+1);%max_WFSloc{length(ind_holds)}(end) - max_WFSloc{length(ind_holds)}(ind_pl_st(i)+1);
    Ploc   = (max_WFPloc{i}(1:end)-max_WFPloc{i}(ind_pl_st(i)+1))/scaleP*100;
    Sloc   = (max_WFSloc{i}(1:end)-max_WFSloc{i}(ind_pl_st(i)+1))/scaleS*100;
    tloc   = timacloc{i}(1:end)-timacloc{i}(ind_pl_st(i)+1);
    Vp_Vs  = Sdelayloc{i}./Pdelayloc{i};
    Vp_Vs_0 = Vp_Vs(ind_pl_st(i)+1);
    Vp_Vs  = (Vp_Vs - Vp_Vs_0)/Vp_Vs_0*100;
    if i == length(ind_holds)
        AX1 = subplot(1,2,1);
        L = semilogx(AX1,tloc,Ploc,[cseq{i} '-'],tloc,Sloc,[cseq{i} '--']);
        set(L(1),'linewidth',2);
        set(L(2),'linewidth',2);
        set(AX1,'NextPlot','Add','Fontsize',16);
        xlabel(AX1,'Time [s]')
        ylabel(AX1,'% \Delta Transmissivity')
        Lsol   = [Lsol L(1)];
        legendstr{i} = [num2str(round(time(ind_reload(i))-time(ind_holds(i)))) 's hold'];
        AX2 = subplot(1,2,2); 
        semilogx(AX2,tloc,Vp_Vs,[cseq{i} '-'],'linewidth',2)
        xlabel(AX2,'Time [s]')
        ylabel(AX2,'% change in V_p/V_s')
        set(AX2,'NextPlot','Add','Fontsize',16,'yaxislocation','right');
    else
        L = semilogx(AX1,tloc,Ploc,[cseq{i} '-'],tloc,Sloc,[cseq{i} '--']);
        semilogx(AX2,tloc,Vp_Vs,[cseq{i} '-'],'linewidth',2)
        set(L(1),'linewidth',2);
        set(L(2),'linewidth',2);
        Lsol   = [Lsol L(1)];
        legendstr{i} = [num2str(round(time(ind_reload(i))-time(ind_holds(i)))) 's hold'];
    end
end
legend(AX1,fliplr(Lsol),legendstr)
%% n-plot %% Plot parameters
figure
x1 = time;
x2 = time_ac_loc;
xlim = [time(ind_holds(1)) time(ind_reload(end))];
ind1 = x1>=xlim(1) & x1<=xlim(2);
ind2 = x2>xlim(1) & x2<=xlim(2);
y2   = max_WFS_loc(ind2);
x2   = x2(ind2);
y1   = mu_s-0.76;
y1   = y1(ind1);
x1   = x1(ind1);
y6   = closures_corr(ind1);
y3   = max_WFP_loc(ind2);%P_delay_loc(ind2)*10^6;
y4   = S_delay_loc(ind2)*10^6;
y4   = y4-y4(1);
y5   = P_delay_loc(ind2)*10^6;
y5   = y5 - y5(1);
y7   = closures(ind1);
%y3  = (y3 - min(y3))/(max(y3)-min(y3)) - 1;
%y2  = (y2 - min(y2))/(max(y2)-min(y2));
%y4   = vf;
%y4   = y4 - min(y4);
xmat = {x1',x2',x2',x1'};
ymat = {(y6-y6(1))',y3',y2',y1'};
nplot = 4;
color = {'k','r','b','k'};
%ylabmat  = {'V [mms^{-1}]','S-amp','P-amp','\Delta\mu'};
ylabmat = {'Closure [\mum]','P-amp','S-amp','\Delta\mu'};
xlab = 'Time [s]';
yon    = [1,1,1,1];
ticprecision = {'%4.3f','%4.3f','%5.4f','%4.3f'};
axisloc = {'right','left','left','right'};
[Axcomb,Lcomb] = make_n_plot(nplot,[0.2 0.1 0.6 0.8],ymat,xmat,color,axisloc,xlab,ylabmat,yon,ticprecision);
set(Lcomb,'linewidth',2);
set(Axcomb,'Fontsize',20,'Fontname','Times','xlim',xlim)
set(Axcomb(1),'NextPlot','Add')
plot(Axcomb(1),x1,y7-y7(1),'c-','linewidth',2)
legend(Axcomb(1),'Poisson Corrected','Raw')
%% 
Pdelay_loc = cell(length(ind_holds),1);
Sdelay_loc = cell(length(ind_holds),1);
figure
for i=1:length(ind_holds)
    subplot(2,3,i)
    pl_loc = ind_holds(i)+round((ind_reload(i)-ind_holds(i))/100):ind_holds(i)+round((ind_reload(i)-ind_holds(i))+1);
    t_loc = time(pl_loc);
    timacloc{i} = time_ac_loc(time_ac_loc>=t_loc(1) & time_ac_loc<=t_loc(end));
    Pdelay_loc{i} = P_delay_loc(time_ac_loc>=t_loc(1) & time_ac_loc<=t_loc(end))*1e6;
    Sdelay_loc{i} = S_delay_loc(time_ac_loc>=t_loc(1) & time_ac_loc<=t_loc(end))*1e6;
    semilogx(timacloc{i}/1000,Pdelay_loc{i}-Pdelay_loc{i}(1),timacloc{i}/1000,(Sdelay_loc{i}-Sdelay_loc{i}(1)));
    xlabel('Time [10^3s]')
    ylabel('Diff. Travel time [\mus]')
    set(gca,'xlim',[min(t_loc) max(t_loc)]/1000);
    title([num2str(round(time(ind_reload(i))-time(ind_holds(i)))) 's hold'])
    legend(gca,'P-\DeltaTr','S-\DeltaTr');
end