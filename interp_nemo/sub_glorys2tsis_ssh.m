% function is identical to glorys2tsis_ssh.m
%
% Interpolate GLORYS into HYCOM-TSIS
% Horizontal grid
% Use bilinear interpolation:
% F= a0+a1*x+a2*y+a3*x*y - basis functions [a0,a1,a2,a3]
% are obtained in calculate_glorys_weights.m
function sshi = sub_glorys2tsis_ssh(dnmb,pthnemo,pthtopo,pthdata,pthglorys,LAT,LON,HH);

[mm,nn] = size(HH);
[XM,YM]=meshgrid([1:nn],[1:mm]);


% Get GLORYS field:
[fglnm,flglrs] = sub_find_GLORYS_file(pthglorys,dnmb);

LONN = nc_varget(flglrs,'longitude');
LATN = nc_varget(flglrs,'latitude');
[LNN,LTN] = meshgrid(LONN,LATN);

ZZN  = nc_varget(flglrs,'depth');

ssh = squeeze(nc_varget(flglrs,'zos'));
ssh0 = ssh;
ssh = sub_fill_land(ssh);

mm = size(LATN,1);
nn = size(LONN,1);
% 
% Interpolation points, polynomial coeff. 
fglrint = sprintf('%sGLRS_TO_HYCOM_TSIS_interp_pnts.mat',pthdata)
fprintf('Loading %s\n',fglrint);
GLR = load(fglrint);

zz0 = 0;
sshi = sub_interp_glorys2hycom(LNN,LTN,LON,LAT,HH,GLR,ssh,zz0);


% ------------------------
% Check & plot:
%
f_chck=0;
if f_chck==1;
% Get NEMO grid:
		f_get_grid=0;
		fgrd = sprintf('%sNEMO_grid.mat',pthnemo);
		if f_get_grid==1
				[LONN,LATN,ZZN] = sub_get_NEMO_grid(dnmb);
				LONN = double(LONN);
				LATN = double(LATN);
				ZZN = double(ZZN);
				fprintf('Saving grid %s\n',fgrd);
				save(fgrd,'LONN','LATN','ZZN');
		else
				fprintf('Loading NEMO grid %s\n',fgrd);
				load(fgrd);
		end

  ln1=min(min(LONN));
  ln2=max(max(LONN));
  lt1=min(min(LATN));
  lt2=max(max(LATN));

  lnmn = min(min(LON));
  lnmx = max(max(LON));
  ltmn = min(min(LAT));
  ltmx = max(max(LAT));


  figure(1); clf;
  axes('Position',[0.07 0.2 0.4 0.7]);
  pcolor(LNN,LTN,ssh0); shading flat;
  hold on;
% NEMO domain:
  plot([ln1 ln2],[lt1 lt1],'k-');
  plot([ln2 ln2],[lt1 lt2],'k-');
  plot([ln1 ln2],[lt2 lt2],'k-');
  plot([ln1 ln1],[lt1 lt2],'k-');

  caxis([-0.4 0.4]);
  colorbar('SouthOutside');
  axis('equal');
%  set(gca,'xlim',[950 1500],...
%          'ylim',[1060 1350]);
  set(gca,'xlim',[lnmn lnmx],...
          'ylim',[ltmn ltmx]);
  set(gca,'xtick',[],...
          'ytick',[]);
  title('GLORYS');

  axes('Position',[0.53 0.2 0.4 0.7]);
  pcolor(LON,LAT,sshi); shading flat;
  hold on;
% NEMO domain:
  plot([ln1 ln2],[lt1 lt1],'k-');
  plot([ln2 ln2],[lt1 lt2],'k-');
  plot([ln1 ln2],[lt2 lt2],'k-');
  plot([ln1 ln1],[lt1 lt2],'k-');

  caxis([-0.4 0.4]);
  colorbar('SouthOutside');
  axis('equal');
%  set(gca,'xlim',[1 1401],...
%          'ylim',[1 891]);
  set(gca,'xlim',[lnmn lnmx],...
          'ylim',[ltmn ltmx]);
  set(gca,'xtick',[],...
          'ytick',[])
%          'xlim',[1 ah2],...
%          'ylim',[1 ah1]);
  title('HYCOM-TSIS');

end

return

