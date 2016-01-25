function [Axcomb,Lcomb] = make_n_plot(nplot,position,ymat,xmat,color,axisloc,xlab,ylabmat,yon,ticprecision)

% Function to generate stacked plots on the top of each other
% Inputs: nplot : Number of Y-axis
%         position : Four element vector for master axis, std 
%                    Matlab convention [left bottom width height]
%         ymat : cell for all ydata
%         xmat : cell for all xdata
%         color : cell with color scheme, each entry specifies color of the
%                 respective axis in order
%         axisloc : cell where each element is either 'right' or 'left',
%                   these specify where each y-axis is.
%         xlab : String for xlabel
%         ylabmat : Cell of strings for each ylabel
%         ticprecision : cell string for precision specs for each Y-axis.
%         yon  : Do you want to leave ylabel: On - 1, Off - 0
%         scal
% Outputs : Axcomb : Vector of all axes handles, includes the outer box,
%                       hence 1 more element than no. of lines
%           Lcomb  : Vector of all lines plotted, 1 per inner axis
% 

%% Section for size specification
left = position(1);
bottom = position(2);
totwid = position(3);
totht  = position(4);
htpart = totht/nplot;
%% Generate figure 
set(gcf,'units','normalized',...
       'DefaultAxesXMinorTick','on','DefaultAxesYminorTick','on');
Axcomb = [];
Lcomb = [];

for i = 1:nplot
    ymin = min(ymat{i});
    ymax = max(ymat{i});
    ywid = ymax-ymin;
    Ytic = ymin:ywid/3:ymax;
    Ylim = [ymin-ywid/10 ymax+ywid/10];
    if yon(i)
        YticStr = sprintf([ticprecision{i} '\n'],Ytic);
    else
        YticStr = [];
    end
    H = axes('position',[left bottom+(i-1)*htpart totwid htpart]);% First subaxis from bottom
        L = plot(xmat{i},ymat{i},'Parent',H,'color',color{i});
        set(H,'box','on','Yaxislocation',axisloc{i},'Ycolor',color{i},...
            'Xcolor','none','Yminortick','on','Xminortick','on',...
            'Xtick',[],'Ytick',Ytic,'ylim',Ylim,'Yticklabel',YticStr);
    if yon(i)
        ylabel(H,ylabmat{i})
    end
    Axcomb = [Axcomb H];
    Lcomb = [Lcomb L];
end
% Final axes overlay on top for bottom xtick 
Hmaster = axes('position',[left bottom totwid totht]);%[left bottom width height
set(Hmaster,'box','on','Color','none','Ytick',[],'Xminortick','on','Ycolor','none');
xlabel(Hmaster,xlab)
Axcomb = [Axcomb Hmaster];
