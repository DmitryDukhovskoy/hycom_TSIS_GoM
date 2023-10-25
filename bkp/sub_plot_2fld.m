function sub_plot_2fld(pos,ssh,LON,LAT,HH,c1,c2,stt,IN,xl1,xl2,yl1,yl2,clb,varargin);
% Plot ssh field
% or other 2D fields
% Colormap:
%c1=-0.5;
%c2=0.5;

cmp = [];
nV = length(varargin);
if nV>0
  for k=1:nV
    vfld = varargin{k};
    if strmatch(vfld,'cmp') % 
      cmp = varargin{k+1}; %
    end
  end
end

if isempty(cmp);
		nc=100;
		clB=flipud(colormap_blue(nc));
		clR=colormap_red(nc);
		clB(nc-3,:)=[1 1 1];
		clB(nc-2,:)=[1 1 1];
		clB(nc-1,:)=[1 1 1];
		clB(nc,:)  =[1 1 1];
		clR(4,:)   =[1 1 1];
		clR(3,:)   =[1 1 1];
		clR(2,:)   =[1 1 1];
		clR(1,:)   =[1 1 1];
		cmp=[clB;clR];
		cmp=smooth_colormap(cmp,5,3);
end
nint=length(cmp);
cnt=(c1:(c2-c1)/nint:c2);  % 



axes('Position',pos);
pcolor(LON,LAT,ssh); shading flat;
caxis([c1 c2]);
colormap(cmp);
hold on;
if ~isempty(HH)
  contour(LON,LAT,HH,[0 0],'Color',[0 0 0]);
end
%  title(['SSH, ',regn,' - ',expt,' Year: ',int2str(year),'; D= ',sday]);
%date_str=sprintf('%4.4i/%2.2i/%2.2i : %2.2ih',DV(1:4));
%stt=sprintf('GLBv0.08_%3.3i, SSH %s',expt,date_str);
title(stt,'Interpreter','none');

% ==================
% Plot ssh
% =================
smm=ssh;
smm(IN==0)=nan;
%contour(LON,LAT,ssh,[0:0.1:1],'w','linewidth',1.);
%contour(LON,LAT,smm,[0.17 0.17],'k','linewidth',1.2);
%contour(LON,LAT,ssh,[-1:0.1:-0.01],'w--','linewidth',1.);

if clb==1
  chb = colorbar('SouthOutside');
  set(chb,'Fontsize',12,...
	         'Position',[0.25 0.11 0.5 0.02],...
          'Ticks',[c1:0.5:c2],...
          'TickLength',0.02);
elseif clb==2
  chb = colorbar('SouthOutside');
  set(chb,'Fontsize',12,...
          'Position',[0.25 0.08 0.5 0.02],...
          'Ticks',[c1:0.5:c2],...
          'TickLength',0.02);
elseif clb==3
  chb = colorbar('SouthOutside');
  set(chb,'Fontsize',12,...
          'Position',[0.25 0.06 0.5 0.02],...
          'Ticks',[c1:0.5:c2],...
          'TickLength',0.02);
end
axis('equal');

%set(gca,'xlim',[-98 -56.08],'ylim',[7.05 31.92]);
%set(gca,'xlim',[-97.71 -78],'ylim',[15 31.1]);
set(gca,'xlim',[xl1 xl2],'ylim',[yl1 yl2]);
set(gca,'tickdir','out',...
	'xtick',[-96:2:-50],...
	'ytick',[4:2:40],...
	'Fontsize',12);

return
