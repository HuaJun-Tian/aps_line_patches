function [scalePstep,offsetPstep,Xstep,Ystep,atm_tropo_S,atm_tropo_O]=strafication_correction_HJ(xy,rate,dem,poly_cell,centers,ts,class,index_class,lonlat,xy_range)
% % apply APS corrction in below
% % % %
% Author: first corrected by Wenyu Gong, 2021-2023; Modfied by Jun Hua,
% 2023.4
% e-mail:huajun@ies.cc.cn
% Release:1
% Release date: 30/04/2023
% option: class, hgt
%% prepare windows splitting
class_id=1;
if nargin < 7
    class=[];
    class_id=nan;

end

patches=numel(poly_cell);
scalePstep = nan(patches,1);
offsetPstep=nan(patches,1);
colorid=1;
con_each_line=ceil(sqrt(patches));
atm_tropo_S=nan(size(rate));
atm_tropo_O=nan(size(rate));
% step_x=nan(size(rate));
% step_y=nan(size(rate));
%% processing for the stratification delay
figure(22);
total_number=0;
total_pha=[];
total_ll=[];
total_dem=[];
for i = 1:patches
    %
    %     keyboard
    poly_single=poly_cell{i};

    if isnan(class_id)
        ix_temp = inpolygon(xy(:,2),xy(:,3),poly_single(:,2),poly_single(:,3));
        totalpix=sum(ix_temp);% total pixels in the patch
        obs_phavec=rate(ix_temp);
        demh=dem(ix_temp);
        %     step_x=xy(ix_temp,2);
        %     step_y=xy(ix_temp,3);
        Xstep(i,1)=centers(i,1);
        Ystep(i,1)=centers(i,2);

        idxNaN = isnan(obs_phavec);
        subplot(con_each_line,con_each_line,i);
        plot(demh(:),obs_phavec(:),'x','color',[1/colorid,0,1/colorid]);hold on
        colorid=colorid+1;

        % estimate the linear regression for each patch
        stdr = std(obs_phavec (~idxNaN));
        meanr=mean(obs_phavec (~idxNaN));
        idxrmv=find(obs_phavec > (meanr)+3*stdr | obs_phavec<(meanr)-3*stdr );
        obs_phavec(idxrmv)=NaN;
        idxNaN = isnan(obs_phavec);
        if  (totalpix-numel(idxNaN) <3) || 1.0*numel(idxNaN)/totalpix>0.6
            disp('Too less obs in the patch')
            scalePstep(i,1)=NaN;
            offsetPstep(i,1)=NaN;%         r_corr= p(1)*de +p(2);
            continue
        end
        %% process the a linear regression
        p=polyfit(demh(~idxNaN),obs_phavec(~idxNaN),1);
        scalePstep(i,1)=p(1);
        offsetPstep(i,1)=p(2);%         r_corr= p(1)*de +p(2);
        y_est = p(1)*demh + p(2);
        atm_tropo_S(ix_temp,1)=p(1);
        atm_tropo_O(ix_temp,1)=p(2);
        plot(demh,y_est,'x','color','red');
        hold off

        TSS=sum((obs_phavec(~idxNaN)-mean(obs_phavec(~idxNaN))).^2);
        RSS=sum((obs_phavec(~idxNaN)-y_est(~idxNaN)).^2);
        R2=1-RSS/TSS;

        clear idx idxNaN
    else
        %         keyboard
        ix_temp=index_class{i}(:,1);

        obs_phavec=rate(ix_temp);
        xy_vec=xy(ix_temp,2:3);
        demh=dem(ix_temp);
        total_pha=[total_pha;obs_phavec];
        total_ll=[total_ll;lonlat(ix_temp,1:2)];
        total_dem=[total_dem;demh];
        if (size(obs_phavec,1)~=size(demh,1))

            keyboard
        end
        totalpix=size(demh,1);
        total_number=total_number+totalpix;
        %         ix_temp = index_class{i}(:,1);
        Xstep(i,1)=centers(i,1);
        Ystep(i,1)=centers(i,2);
        %         keyboard
        rate_temp=rate(ix_temp);
        dem_temp=single(dem(ix_temp));
        if (norm(obs_phavec-rate_temp)==0)



            disp('using ISODATA method');

            if (numel(xy_range)>1)
                idx_temp=find(xy_vec(:,1)>=xy_range(1,1) & xy_vec(:,1)<=xy_range(2,1) & ...
                    xy_vec(:,2)>=xy_range(1,2) & xy_vec(:,2)<=xy_range(2,2));
                idxNaN = unique([find(isnan(obs_phavec)); idx_temp]);
                clear idx_temp
            else
                idxNaN = isnan(obs_phavec);

            end
            obs_phavec(idxNaN)=nan;
            subplot(con_each_line,con_each_line,i);
            plot(demh(:),obs_phavec(:),'x','color',[1/colorid,0,1/colorid]);hold on
            colorindex_classid=colorid+1;

            % estimate the linear regression for each patch
            stdr = std(obs_phavec (~idxNaN));
            meanr=mean(obs_phavec (~idxNaN));
            idxrmv=find(obs_phavec > (meanr)+3*stdr | obs_phavec<(meanr)-3*stdr );
            obs_phavec(idxrmv)=NaN;
            idxNaN = isnan(obs_phavec);
            if  totalpix<3 || 1.0*numel(idxrmv)/totalpix>0.8
                disp('Too less obs in the patch')
                scalePstep(i,1)=NaN;
                offsetPstep(i,1)=NaN;%         r_corr= p(1)*de +p(2);
                continue
            end
            %% process the a linear regression
%             tol=1.0e-15;
%             maxit=500;
%             B_coe=[demh(~idxNaN),ones(size(demh(~idxNaN),1),1)];
%             obs_par=obs_phavec(~idxNaN);
% 
%             nbb = sparse(B_coe'*B_coe);
%             d_part=sparse(B_coe'*(obs_par));
%             step_uli.type='crout';
%             step_uli.droptol=tol;
%             step_uli.milu='row';
%             [l1,u1]=ilu(nbb,step_uli);
%             p=bicg(nbb,d_part,tol,maxit,l1,u1);
%             clear B_coe obs_par nbb d_part step_uli l1 u1
            p=polyfit(demh(~idxNaN),obs_phavec(~idxNaN),1);
            scalePstep(i,1)=p(1);
            offsetPstep(i,1)=p(2);%         r_corr= p(1)*de +p(2);
            y_est = p(1)*demh + p(2);
            atm_tropo_S(ix_temp,1)=p(1);
            atm_tropo_O(ix_temp,1)=p(2);
            plot(demh,y_est,'x','color','red');
            hold off

            TSS=sum((obs_phavec(~idxNaN)-mean(obs_phavec(~idxNaN))).^2);
            RSS=sum((obs_phavec(~idxNaN)-y_est(~idxNaN)).^2);
            R2=1-RSS/TSS;

            clear idx idxNaN rate_temp   dem_temp 

        end


    end

end
if ~exist('aps_patches','dir')
    mkdir('aps_patches')
end
saveas(gcf,['aps_patches/dem_cor_',num2str(ts,'%6d')],'fig');
print(['aps_patches/dem_cor_',num2str(ts,'%6d')],gcf,'-djpeg','-r150')
fprintf('total number: %f\n',total_number);
close all
% keyboard
% figure();hold on;
% scatter(total_ll(:,1),total_ll(:,2),[],total_pha(:,1));colorbar;title(['originnal']);colormap('jet');
% % scatter(temp_mean(:,1),temp_mean(:,2),"filled",'^','Color','red');
% % text(temp_mean(:,1),temp_mean(:,2),string_patch,'Color','red');
% hold off;
% figure();hold on;
% scatter(lonlat(:,1),lonlat(:,2),[],atm_tropo_S);colorbar;title(['originnal']);colormap('jet');
% % scatter(temp_mean(:,1),temp_mean(:,2),"filled",'^','Color','red');
% % text(temp_mean(:,1),temp_mean(:,2),string_patch,'Color','red');
% hold off;
% fprintf('pause,please enter any keyboard\n')
% pause;
% close all
end
