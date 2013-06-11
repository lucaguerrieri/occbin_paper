function makechart(titlelist,legendlist,figlabel,zdata1,zdata2,zdata3,nobsmat,ylabels)

figure

titlelist = char(strrep(cellstr(titlelist),'_','.'));

ndsets=3;       % default, changed below as applicable
if nargin==4
    zdata2=nan*zdata1;
    zdata3=nan*zdata1;
    ndsets =1;
    nobsmat = size(zdata1,1)*ones(size(zdata1,2));
    ylabelflag = 0;
elseif (nargin == 5)
    zdata3 =nan*zdata1;
    ndsets=2;
    nobsmat = size(zdata1,1)*ones(size(zdata1,2));
    ylabelflag = 0;
elseif (nargin == 6)
    nobsmat = size(zdata1,1)*ones(size(zdata1,2));
    ylabelflag = 0 ;
elseif (nargin == 7)
    ylabelflag = 0;
elseif (nargin == 8)
    ylabelflag = 1;
elseif ((nargin>8) | (nargin <=3))
    error ('makechart takes 4 to 8 arguments')
end

if isempty(zdata2) & isempty(zdata3)
    ndsets = 1;
    zdata2=nan*zdata1;
    zdata3=nan*zdata1;
elseif isempty(zdata3)
    ndsets = 2
    zdata3=nan*zdata1;
end
    
nobs = size(zdata1,1);
xvalues = (1:nobs)';

nvars = size(titlelist,1);
if nvars==1 
    nrows=1;
    ncols = 1;
elseif nvars==2
    nrows =2;
    ncols = 1;
elseif (nvars == 3 | nvars ==4)
    nrows = 2;
    ncols =2;
elseif (nvars==5 |nvars ==6)
    nrows = 3;
    ncols = 2;
elseif (nvars==7 | nvars==8)
    nrows = 4;
    ncols = 2;
else 
    error('too many variables (makechart)')
end

for i = 1:nvars
    subplot(nrows,ncols,i)
    if i==1
        h1=plot(xvalues(1:nobsmat(2)),zdata1(1:nobsmat(2),i),'k-',xvalues(1:nobsmat(2)),zdata2(1:nobsmat(2),i),'r--',xvalues(1:nobsmat(2)),zdata3(1:nobsmat(2),i),'b:');
    [x0 x1 y10 y11] = pickaxes(xvalues(1:nobsmat(i)),zdata1(1:nobsmat(i),i));
    [x0 x1 y20 y21] = pickaxes(xvalues(1:nobsmat(i)),zdata2(1:nobsmat(i),i));
    [x0 x1 y30 y31] = pickaxes(xvalues(1:nobsmat(i)),zdata3(1:nobsmat(i),i));
        set(h1,'linewidth',2);
        hold on
        h2=plot(xvalues(nobsmat(2)+1:nobsmat(1)),zdata1(nobsmat(2)+1:nobsmat(1),2),'k-',xvalues(nobsmat(2)+1:nobsmat(1)),zdata2(nobsmat(2)+1:nobsmat(1),2),'r--',xvalues(nobsmat(2)+1:nobsmat(1)),zdata3(nobsmat(2)+1:nobsmat(1),2),'b:');
        set(h2,'linewidth',2);
        hold off
    else
     h1=plot(xvalues(1:nobsmat(i)),zdata1(1:nobsmat(i),i),'k-',xvalues(1:nobsmat(i)),zdata2(1:nobsmat(i),i),'r--',xvalues(1:nobsmat(i)),zdata3(1:nobsmat(i),i),'b:');
    [x0 x1 y10 y11] = pickaxes(xvalues(1:nobsmat(i)),zdata1(1:nobsmat(i),i));
    [x0 x1 y20 y21] = pickaxes(xvalues(1:nobsmat(i)),zdata2(1:nobsmat(i),i));
    [x0 x1 y30 y31] = pickaxes(xvalues(1:nobsmat(i)),zdata3(1:nobsmat(i),i));
        set(h1,'linewidth',2);
    end
    
    y0 = min([y10,y20,y30]);
    y1 = max([y11,y21,y31]);
    if y0==y1
        y1=y0+1;
    end
    
    axis([x0 x1 y0 y1])
    
    if i==1
        legend(legendlist)
        text('String',figlabel,'Units','normalized','Position',[1.2 1.24],...
       'FontSize',14,'FontWeight','bold','HorizontalAlignment','center');
    end
 
%     set(gca,'XTick',xtick)
%     set(gca,'XTickLabel',xticklabel)
    if ylabelflag
       ylabel(ylabels(i,:));
    end
    if i==nvars | i==nvars-1
        xlabel('Quarters');
    end
    title([num2str(i),'. ',titlelist(i,:)]);
    
end

% sets printing preferences
printpref
