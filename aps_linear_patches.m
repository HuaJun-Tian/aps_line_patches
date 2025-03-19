function []=aps_linear_patches(stamps_processed,non_defo_flag,hgt_matfile,ll_matfile,phuw_matfile,lon_range,lat_range,isodata_parameter)
% [ph_tropo_linear] = aps_linear_patches(save_path)
% Scipt to compute the tropospheric delay map from a linear relation
% between phase and topography in a patch, optional a non-deforming polygon can be
% specified by setting crop_flag to 'y' in the parms_aps list.
% The computed tropospheric delay is in the same units as the
% inputed interferogram phases.
%
% All required inputs will be read from the aps_parm list. This includes:
% non_defo_flag 'y' or 'n' to use a non-deformaing region.
%               Polygon of the non-defor ming area. By default 'n' the whole
%               interferogram is used. Change to 'y' by using setparm_aps.
%               Note that this variable is a matrix with in its columns the
%               longitude and latitude of the non-deforming area.
% hgt_matfile   Path to the interferograms
%               Interferogram phases, given as a column vector or matrix
%               with in its columens the different interferograms.
%               Stamps structure will automatically be recognised. Use
%               setparm_aps so change the data path pointing to the .mat
%               file.
% hgt_matfile   Path to the heights file.
%               Colum vector of the topography.
%               Stamps structure will automatically be recognised. Use
%               setparm_aps so change the data path pointing to the .mat
%               file.
% ll_matfile    Path to the longitude and latitude file
%               Matrix with in its columns the longitude and latitude.
%               Stamps structure will automatically be recognised. Use
%               setparm_aps so change the data path pointing to the .mat
%               file. In case the poly argument is specified, both need to
%               have the same units.
%
% OUTPUTS:
% ph_tropo_linear   The topography correlated delays either estimated from a
%                   linear relationship over the whole interferogram of using
%                   a non-deforming region. By default the output of this
%                   function is stored in the 'tca2.mat' or 'tca_sb2.mat'
%                   for StaMPS SM and SB option. In case no StaMPs
%                   structure is used the data is saved in 'tca2.mat'.
%
%   HuA JUN 2023/5
%


test_fig = 0;   % debug flag which plots the scatter cloud and the
% estimated line
n_fig_line = 7; % number of ifgs per row for the plots
fontsize =10;   % figure fontsize
addpath(genpath('./'));
save_path=['.'];



% loading the data
phuw = load(phuw_matfile);
lonlat =load(ll_matfile);
hgt = load(hgt_matfile);
psver =2;

%%################################################
%%################################################

phuw = phuw.phuw_single;
lonlat = lonlat.lonlat_single;
hgt = hgt.hgt_single;


%% Loading of the data
% file names of the output data
apsname = [save_path filesep 'tca' num2str(psver) '.mat'];
% apssbname = [save_path filesep 'tca_sb' num2str(psver) '.mat'];


% the number of interferograms
n_dates = size(phuw,2);

n_points=size(phuw,1);
%% use a non-deforming area
if non_defo_flag==1
    %     non_defo = load('non_defo.mat');
    %     poly = non_defo.poly;

    % search those points within the coseismic region
    ixnon_points = [1:size(hgt,1)]';
    ix_temp=find(lonlat(:,1)>=lon_range(1) & lonlat(:,1)<=lon_range(2) & ...
        lonlat(:,2)>=lat_range(1) & lonlat(:,2)<=lat_range(2));
    ixnon_points=ixnon_points(ix_temp);
    ix_points=[1:size(hgt,1)]';
    ix_points(ix_temp)=[];
    clear ix_temp
else
    % use all points
    ix_points = [1:size(hgt,1)]';
    ixnon_points=[];
end

% figure()
% scatter(lonlat(:,1),lonlat(:,2),[],phuw(:,3));colorbar;
fprintf('interp method: 1 linear, 2 kring\n');
% imethod=input();


%% Compute the linear relation between phase and topography for each interferogram
% and compute the tropospheric delay for the full interferogram from it.
% initialisation of the tropospheric delay matrix
ph_tropo_linear = zeros([size(hgt,1) n_dates]);
hgt_range = [min(hgt) max(hgt)]';
if hgt_range(2)>10
    % height are in m
    hgt = hgt/1000;       % km
    hgt_range = hgt_range/1000.0;
end

xy = (llh2local(lonlat',mean(lonlat)))';
xy = [[1:size(xy,1)]' xy];

if numel(lon_range) >1
xy_range=(llh2local([lon_range; lat_range],mean(lonlat)))';

else
    xy_range=1;
end
% ps.xy = xy;
% ps.n_ifg = size(phuw,2);
% ps.n_ps = size(phuw,1);


ix_points_or = ix_points;
ixnon_points_or = ixnon_points;
aps_val=load(apsname);
if ~exist('aps_patches','dir')
    mkdir('aps_patches')
end
% if (~isfield(aps_val,'ph_tropo_gacos'))
% [nb, err, nPts, centers, dhgt, polys, xLims, yLims] = quadtree(xy, hgt, 0.5, 1000, 1); % Run Quadtree on DEM
% drawnow
% end
for k=1:n_dates

    ixnon_points = ixnon_points_or;
    ix_points = ix_points_or;

    % removing NaN's from the phase
    ixnon_points = unique([ixnon_points ; find(isnan(phuw(:,k))==1)]);

    ix_points = [1:n_points]';
    ix_points(ixnon_points)=[];
    if (isfield(aps_val,'ph_tropo_gacos'))
        gacos_tropo=aps_val.ph_tropo_gacos(:,k);
        %         [nb, err, nPts, centers, dhgt, polys, xLims, yLims] = quadtree(xy(ix_points,:), gacos_tropo(ix_points), 0.35, 1000, 1); % Run Quadtree on previous atm
        %         drawnow
%         [centers,polys,class,index_class]=ISODATA([xy(ix_points,:), gacos_tropo(ix_points),hgt(ix_points),phuw(ix_points,k)],12,1000,5.0e-4,1.8,2,200,[phuw(ix_points,k),hgt(ix_points)],lonlat(ix_points,1:2));
%         [centers,polys,class,index_class]=ISODATA([xy(ix_points,:),gacos_tropo(ix_points)],30,3000,5.0e-5,2.8,2,200,[phuw(ix_points,k),hgt(ix_points)],lonlat(ix_points,1:2));
%          polys=1;
%         save test;
        if isstruct(isodata_parameter) & isodata_parameter.flag==1

            [centers,polys,class,index_class]=ISODATA([xy(ix_points,:),gacos_tropo(ix_points)],isodata_parameter.K,isodata_parameter.theta_N,isodata_parameter.theta_S,isodata_parameter.theta_c,isodata_parameter.L,200,[phuw(ix_points,k),hgt(ix_points)],lonlat(ix_points,1:2));
        else
            [centers,polys,class,index_class]=ISODATA([xy(ix_points,:),gacos_tropo(ix_points)],30,3000,5.0e-5,2.8,2,isodata_parameter.I,[phuw(ix_points,k),hgt(ix_points)],lonlat(ix_points,1:2));

        end
        [scalePstep,offsetPstep,Xstep,Ystep,atm_tropo_S,atm_tropo_O]=strafication_correction_HJ(xy(ix_points,:),phuw(ix_points,k),hgt(ix_points),polys,centers,k,class,index_class,lonlat(ix_points,1:2),xy_range);
        clear gacos_tropo scalePstep offsetPstep Xstep Ystep class index_class polys centers;
    else
        [scalePstep,offsetPstep,Xstep,Ystep,atm_tropo_S,atm_tropo_O]=strafication_correction_HJ(xy(ix_points,:),phuw(ix_points,k),hgt(ix_points),polys,centers,k);
    end

    imethod=2; % 2 use linear interpolation, else using krigging
    if imethod ==2
        %         xy_temp=xy
        %         scalePinterp0=griddata(Xstep(~isnan(scalePstep)),Ystep(~isnan(scalePstep)),scalePstep(~isnan(scalePstep)),xy(:,2),xy(:,3),'linear');
        %         scalePinterp=griddata(xy(~isnan(scalePinterp0),2),xy(~isnan(scalePinterp0),3),scalePinterp0(~isnan(scalePinterp0)),xy(:,2),xy(:,3),'nearest');
        %         offsetPinterp0=griddata(Xstep(~isnan(offsetPstep)),Ystep(~isnan(offsetPstep)),offsetPstep(~isnan(offsetPstep)),xy(:,2),xy(:,3),'linear');
        %         offsetPinterp=griddata(xy(~isnan(offsetPinterp0),2),xy(~isnan(offsetPinterp0),3),scalePinterp0(~isnan(offsetPinterp0)),xy(:,2),xy(:,3),'nearest');
        scalePinterp=griddata(xy(ix_points,2),xy(ix_points,3),atm_tropo_S,xy(:,2),xy(:,3),'linear');
        % %       scalePinterp=griddata(xy(~isnan(scalePinterp0),2),xy(~isnan(scalePinterp0),3),atm_tropo_S(~isnan(scalePinterp0)),xy(:,2),xy(:,3),'nearest');
        offsetPinterp=griddata(xy(ix_points,2),xy(ix_points,3),atm_tropo_O,xy(:,2),xy(:,3),'linear');
        % %       offsetPinterp=griddata(xy(~isnan(offsetPinterp0),2),xy(~isnan(offsetPinterp0),3),atm_tropo_O(~isnan(offsetPinterp0)),xy(:,2),xy(:,3),'nearest');

        %         clear

        %
    else
        if imethod==1
            keyboard

            S=fitVariogram_HJ([xy(ix_points,2),xy(ix_points,3),phuw(ix_points,k)]);
            scalePinterp=kriging(S,xy(ix_points,2),xy(ix_points,3),atm_tropo_S,xy(:,2),xy(:,3));

            offsetPinterp=kriging(S,xy(ix_points,2),xy(ix_points,3),atm_tropo_O,xy(:,2),xy(:,3));
        end


    end
%     save test
    clear  atm_tropo_O atm_tropo_S ixnon_points
    %     figure();
    %     %     offsetPinterp=atm_tropo_O;
    %     %     scalePinterp=atm_tropo_S;
    %     subplot(121);scatter(xy(~isnan(offsetPinterp),2),xy(~isnan(offsetPinterp),3),[],offsetPinterp(~isnan(offsetPinterp)));colorbar;title('Offset factor map');a=caxis;
    %     atmp=colormap;atmp=[1,1,1;atmp];colormap(atmp);a=[a(1)-abs(a(1))*0.3,a(2)-abs(a(2))*0.1];caxis(a)
    %
    %     subplot(122);scatter(xy(~isnan(scalePinterp),2),xy(~isnan(scalePinterp),3),[],offsetPinterp(~isnan(scalePinterp)));title('Scale factor map');colorbar;a=caxis;
    %     atmp=colormap;atmp=[1,1,1;atmp];colormap(atmp);a=[a(1)-abs(a(1))*0.3,a(2)-abs(a(2))*0.1];caxis(a)
    %
    %     saveas(gcf,['aps_patches/off_scale_before_fielter',num2str(k,'%6d')],'fig');
    %     print(['aps_patches/off_scale_before_fielter',num2str(k,'%6d')],gcf,'-djpeg','-r150')
    close all
    % clear xy
    %% filter
    %     [scalePinterp,offsetPinterp]=gaussian_int(xy(:,2:3),[scalePinterp,offsetPinterp],[8 8],0.2,4);
    %     %     offsetPinterp=gaussian_int(xy(:,2:3),offsetPinterp,[8 8],0.2,4);
    % %     save test
        figure(21);
        %     offsetPinterp=atm_tropo_O;
        %     scalePinterp=atm_tropo_S;
        subplot(121);scatter(xy(~isnan(offsetPinterp),2),xy(~isnan(offsetPinterp),3),[],offsetPinterp(~isnan(offsetPinterp)));colorbar;title('Offset factor map');a=caxis;
        %     atmp=colormap;atmp=[1,1,1;atmp];colormap(atmp);a=[a(1)-abs(a(1))*0.3,a(2)-abs(a(2))*0.1];caxis(a)
        subplot(122);scatter(xy(~isnan(scalePinterp),2),xy(~isnan(scalePinterp),3),[],scalePinterp(~isnan(scalePinterp)));title('Scale factor map');colorbar;a=caxis;
        %     atmp=colormap;atmp=[1,1,1;atmp];colormap(atmp);a=[a(1)-abs(a(1))*0.3,a(2)-abs(a(2))*0.1];caxis(a)
        %     fprintf('pause,please enter any keyboard')
        %     pause;
    
%         saveas(gcf,['aps_patches/off_scale_after_fielter',num2str(k,'%6d')],'fig');
        print(['aps_patches/off_scale_after_fielter',num2str(k,'%6d')],gcf,'-djpeg','-r150')
    % %     keyboard
    %     close all
    % computation of the delay   
%     clear xy

    ph_tropo_linear(:,k) = scalePinterp.*hgt+offsetPinterp;

    % set those pixels not used back to NaN
    ph_tropo_linear(isnan(phuw(:,k)),k)=NaN;
    %     imagesc( ph_tropo_linear(:,k));
    figure()
    scatter(lonlat(:,1),lonlat(:,2),[],phuw(:,k)-ph_tropo_linear(:,k));colorbar;title('Phuw-Tropo');colormap('jet');%caxis([-100 20])
    set(gcf,'PaperSize',[5,7],'PaperUnits','centimeters');box on;
    saveas(gcf,['aps_patches/phuw_cor_',num2str(k,'%6d')],'fig');
    print(['aps_patches/phuw_cor_',num2str(k,'%6d')],gcf,'-djpeg','-r150')
    close all
    figure()
    scatter(lonlat(ix_points,1),lonlat(ix_points,2),[],phuw(ix_points,k));colorbar;title('Phuw-defor');colormap('jet');%caxis([-100 20])
    set(gcf,'PaperSize',[5,7],'PaperUnits','centimeters');box on;
    saveas(gcf,['aps_patches/phuw_mask_',num2str(k,'%6d')],'fig');
    print(['aps_patches/phuw_mask_',num2str(k,'%6d')],gcf,'-djpeg','-r150')
    figure()
    close all
    scatter(lonlat(:,1),lonlat(:,2),[],ph_tropo_linear(:,k));colorbar;title('Tropo');colormap('jet');%caxis([-100 20])
    set(gcf,'PaperSize',[5,7],'PaperUnits','centimeters');box on;
    saveas(gcf,['aps_patches/aps_tropo_',num2str(k,'%6d')],'fig');
    print(['aps_patches/aps_tropo_',num2str(k,'%6d')],gcf,'-djpeg','-r150')
    close all
    %% statisitcs
    idx0=find(~isnan(phuw(ix_points,k)));
    figure();
    histogram(phuw(idx0,k)-ph_tropo_linear(idx0,k),100,'FaceColor',[0.67843 0.921568 1]);hold on;
    histogram(phuw(idx0,k),100,'FaceAlpha',0.5,'FaceColor','none','LineWidth',1.2); hold off;
    legend('Stra. Corr', 'Original');
    saveas(gcf,['aps_patches/statis_hist_',num2str(k,'%6d')],'fig');
    print(gcf,'-djpeg','-r150',['aps_patches/statis_hist_',num2str(k,'%6d'),'.jpg']);
    close all
    if k==1
        fid=fopen('aps_patches/before_after_correction_stats.txt','wt+');
        fprintf(fid,'      befor correction              after correction\n');
        fprintf(fid,'        mean    std                    mean    std  (rad)\n');
        fprintf(fid,'        %f         %f              %f         %f\n',mean(phuw(idx0,k)),std(phuw(idx0,k)),mean(phuw(idx0,k)-ph_tropo_linear(idx0,k)),std(phuw(idx0,k)-ph_tropo_linear(idx0,k)));
        %         fclose(fid);
    else

        %         fid=fopen('aps_patches/before_after_correction_stats.txt','at+');
        fprintf(fid,'        %f         %f              %f         %f\n',mean(phuw(idx0,k)),std(phuw(idx0,k)),mean(phuw(idx0,k)-ph_tropo_linear(idx0,k)),std(phuw(idx0,k)-ph_tropo_linear(idx0,k)));

    end
    clear idx0;
    close all
    if test_fig==1
        if k==1
            figure('name','Linear relation between phase and topography')
            if n_dates/n_fig_line<1
                n_rows = 1;
                n_columns = n_dates;
            else
                n_rows = ceil(n_dates/n_fig_line);
                n_columns = n_fig_line;
            end

        end
        subplot(n_rows,n_columns,k)
        plot(hgt(ix_points),phuw(ix_points,k),'g.')
        hold on
        plot(hgt(ixnon_points),phuw(ixnon_points,k),'k.')
        hold on
        plot(hgt,ph_tropo_linear(:,k),'r-','linewidth',2)
        xlim(hgt_range)
        xlabel('Height','fontsize',fontsize)
        ylabel('Phase','fontsize',fontsize)
        set(gca,'fontsize',fontsize)
    end
end
ph_tropo_linear_patches=ph_tropo_linear;
fclose(fid);
% aps_save(apsname,ph_tropo_linear_patches)
figure()
scatter(lonlat(ix_points,1),lonlat(ix_points,2),[],hgt(ix_points));colorbar;title('hgt');colormap('jet');%caxis([-100 20])
set(gcf,'PaperSize',[5,7],'PaperUnits','centimeters');box on;
saveas(gcf,['aps_patches/hgt'],'fig');
print(['aps_patches/hgt.jpg'],gcf,'-djpeg','-r150')
close all
save aps_patches/aps_Phuw_tropo.mat phuw ph_tropo_linear_patches
%% saving the data
% checking if this is StaMPS or not
% if strcmp(stamps_processed,'y')
% This is StaMPS
%     if strcmp(getparm('small_baseline_flag'),'y')
aps_save(apsname,ph_tropo_linear_patches)
phuw_cor=[lonlat phuw-ph_tropo_linear_patches];
save phuw_cor.txt -ascii phuw_cor
%     else
%         aps_save(apsname,ph_tropo_linear_patches)
%     end
% else
% This is not StaMPS
%     aps_save(apsname,ph_tropo_linear_patches)
%     phuw_cor=[lonlat phuw-ph_tropo_linear_patches];

%     save phuw_cor.txt -ascii phuw_cor
% end
end



