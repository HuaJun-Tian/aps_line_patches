function V = fitVariogram_HJ(insar_data)
%% the programe was used to calculate the std of displacements 
%insar_data arrayx3[lon,lat,dis]
% author: HuaJun, 2022/6/7
% corrected by HuaJun . 2023/5/15
% add V structruct to preare for kring
% Extract subset from rectangular area or after masking
subset = insar_data(:,3);
llon = insar_data(:,1);
llat = insar_data(:,2);
% central point as reference point

% Display subregion from selection
% figure('Position', [1, 1, 1200, 1000]);
% subplot(2,3,1)
% scatter(llon(:),llat(:),[],subset(:),'.')
% colormap('jet')
% subset_max=max(abs(min(subset(:))),abs(max(subset(:))));
% caxis([-subset_max subset_max])
% axis xy
% axis equal
% axis tight
% title('Selected region, NON-DETRENDED')
% xlabel('Longitude (degrees)')
% ylabel('Latitude (degrees)')
% colorbar

%% Remove linear trend from subregion
xy=[llon llat];


A = [xy ones([length(xy) 1])];

coeff = lscov(A,subset);
deramped = subset - A*coeff;

%% Display trend and subregion after removal of trend
subplot(2,3,2)
scatter(llon(:),llat(:),[],A(:,:)*coeff,'.')
colormap('jet')
A_max=max(abs(min(A(:,:)*coeff)),abs(max(A(:,:)*coeff)));
caxis([-A_max A_max])
axis xy
axis equal
axis tight
title('Selected region, ESTIMATED TREND')
% xlabel('Longitude (degrees)')
% ylabel('Latitude (degrees)')
colorbar

subplot(2,3,3)
scatter(llon(:),llat(:),[],deramped(:),'.')
colormap('jet')

der_max=max(abs(min(deramped(:))),abs(max(deramped(:))));
caxis([-der_max der_max])
axis xy
axis equal
axis tight
title('Selected region, DETRENDED')
% xlabel('Longitude (degrees)')
% ylabel('Latitude (degrees)')
colorbar

%% Calculate and display variogram before plane removal
subplot(2,3,4)
variog = variogram(xy,double(subset),'plotit',true,'subsample',3000);
title('Semi-variogram, NON-DETRENDED')

% Calculate and display variogram after detrending
variogDtrnd = variogram(xy,double(deramped),'plotit',false,'subsample',3000,'nrbins',30);

%% Fit exponential function to experimental variogram and display
subplot(2,3,5)
[a,c,n,V] = variogramfit(variogDtrnd.distance,variogDtrnd.val,20000,1e-04,variogDtrnd.num, 'model', 'exponential', 'nugget', 1);
title('Semi-variogram and fit, DETRENDED')
thresh=c-n;
h =subplot(2,3,6);
set(h,'visible','off')
text(0.1,1.0,'Fitted exponential semi-variogram parameters:','FontSize',14)
text(0.1,0.8,['Sill:  ', num2str(c)],'FontSize',14)
text(0.1,0.6,['Range:  ', num2str(a)],'FontSize',14)
text(0.1,0.4,['Nugget:  ', num2str(n)],'FontSize',14)

% Print variogram exponential fit parameters to screen
disp(['Sill:  ',num2str(c)])
disp(['Range:  ',num2str(a)])
disp(['Nugget:  ',num2str(n)])


