clc,clear,close all % 清理命令区、清理工作区、关闭显示图形
% feature jit off % 加速代码运行
origin=textread('20211202_20220908.quad.unw.txt');
result=origin;
%%figure befor adf
figure();
scatter(origin(:,1),origin(:,2),[],origin(:,3),'filled','s');%caxis([-0.3 0.3])
axis equal; axis tight;
cbar = colorbar; 
colormap('jet')
xlabel('longititude (degree)','FontSize', 14)
ylabel('lattitude (degree)','FontSize', 14)
% subset via origin
% polyMask=impoly;
% pos=getPosition(polyMask);
% in = inpolygon(origin(:,1),origin(:,2),pos(:,1),pos(:,2));
% % ixSubset = find(in == 0);
% sub1=origin(find(in == 0),:);
% sub2=origin(find(in == 1),:);
% figure();
% scatter(sub1(:,1),sub1(:,2),[],sub1(:,3),'filled','s');caxis([-0.3 0.3])
% axis equal; axis tight;
% cbar = colorbar; 
% colormap('jet')
% xlabel('longititude (degree)','FontSize', 14)
% ylabel('lattitude (degree)','FontSize', 14)
% figure();
% scatter(sub2(:,1),sub2(:,2),[],sub2(:,3),'filled','s');caxis([-0.3 0.3])
% axis equal; axis tight;
% cbar = colorbar; 
% colormap('jet')
% xlabel('longititude (degree)','FontSize', 14)
% ylabel('lattitude (degree)','FontSize', 14)
% % % adf
interval=0.1; % 经纬度间隔
% im1 = adaptsmooth_filter_near( sub2,[10,10],interval ); % near 应用自适应平滑滤波
% im2 = adaptsmooth_filter_far( sub1,[2,2],interval ); % near 应用自适应平滑滤波

% restore real shear
% result(find(in == 0),:)=im2;
% result(find(in == 1),:)=im1;
% im=result;
im = adaptsmooth_filter_near( origin,[8,8],interval ); % near 应用自适应平滑滤波

%% figure
% shear
figure('color',[1,1,1])
scatter(im(:,1),im(:,2),[],im(:,3),'filled','s');%caxis([-0.3 0.3])
axis equal; axis tight;
cbar = colorbar; 
colormap('jet')
xlabel('longititude (degree)','FontSize', 14)
ylabel('lattitude (degree)','FontSize', 14)
% dilatation
figure('color',[1,1,1])
scatter(im(:,1),im(:,2),[],im(:,3),'filled','s');caxis([-0.3 0.3])
axis equal; axis tight;
cbar = colorbar; 
colormap('jet')
xlabel('longititude (degree)','FontSize', 14)
ylabel('lattitude (degree)','FontSize', 14)
% rotation
figure('color',[1,1,1])
scatter(im(:,1),im(:,2),[],im(:,3),'filled','s');caxis([-0.15 0.15])
axis equal; axis tight;
cbar = colorbar; 
colormap('jet')
xlabel('longititude (degree)','FontSize', 14)
ylabel('lattitude (degree)','FontSize', 14)
save adf_data
%% write and save 
file_name='result/20211202_20220908.quad.unw.adf.txt';
fid=fopen(file_name,'wt+');
fprintf(fid,' %9.4f %9.4f %9.4f\n',im(:,1:3)');
% %     fprintf(fid,'%9.4f%10.4f% 9.4f% 9.4f% 9.4f% 9.4f\n',los1{i}(:,1:6)');
% fprintf(fid,' %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n',im(:,1:10)');
% 
% fclose(fid);


%% function 
function Z = adaptsmooth_filter_near(X,mn,interval)
% 函数对输入的二维图像矩阵进行自适应平滑滤波
% input：
% X：输入的二维图像矩阵
% m：lon的滤波模板
% n：lat的滤波末班
% output:
% Z：输出对m x n的二维图像矩阵的运算结果
if nargin < 2
    m = 3; % 滤波模板尺寸
    n = 3;
else
    m = mn(1);
    n = mn(2); % 模板大小
end
[n1,n2] =size(X);
Z=zeros(n1,n2)*nan;
% mask coffecient
for i=1:n1
%     save test
    H_subset=find( X(:,1)>=X(i,1)-m*interval & X(:,1)<=X(i,1)+m*interval & X(:,2)>=X(i,2)-n*interval & X(:,2)<=X(i,2)+2*interval);
    H=X(H_subset~=i,:);
    Glon=0.5*(H(:,1)-X(i,1));
    Glat=0.5*(H(:,2)-X(i,2));
    d=sqrt( Glon.^2 + Glat.^2);
    Coff=exp(-d./2.0);
    Coff_kuozhan=repmat(Coff,[1,1]);
    if numel(H_subset) > 1
        Z(i,1:2)=X(i,1:2);
%         Z(i,3)=sum(X(i,3).*Coff)/(sum(Coff));
%         Z(i,4)=sum(X(i,4).*Coff)/(sum(Coff));
%         Z(i,5)=sum(X(i,5).*Coff)/(sum(Coff));
%         Z(i,6)=sum(X(i,6).*Coff)/(sum(Coff));
%         Z(i,7)=sum(X(i,7).*Coff)/(sum(Coff));
%         Z(i,8)=sum(X(i,8).*Coff)/(sum(Coff));
%         Z(i,9)=sum(X(i,9).*Coff)/(sum(Coff));
%         Z(i,10)=sum(X(i,10).*Coff)/(sum(Coff));
        Z(i,3)=sum(X(i,3).*Coff_kuozhan)./(sum(Coff_kuozhan));
    end
    
    
end


end



function Z = adaptsmooth_filter_far(X,mn,interval)
% 函数对输入的二维图像矩阵进行自适应平滑滤波
% input：
% X：输入的二维图像矩阵
% m：lon的滤波模板
% n：lat的滤波末班
% output:
% Z：输出对m x n的二维图像矩阵的运算结果
if nargin < 2
    m = 3; % 滤波模板尺寸
    n = 3;
else
    m = mn(1);
    n = mn(2); % 模板大小
end
[n1,n2] =size(X);
Z=zeros(n1,n2)*nan;
% mask coffecient
for i=1:n1
%     save test
    H_subset=find( X(:,1)>=X(i,1)-m*interval & X(:,1)<=X(i,1)+m*interval & X(:,2)>=X(i,2)-n*interval & X(:,2)<=X(i,2)+2*interval);
    H=X(H_subset,:);
    Glon=0.5*(H(:,1)-X(i,1));
    Glat=0.5*(H(:,2)-X(i,2));
    d=sqrt( Glon.^2 + Glat.^2);
%     Coff=exp(-d./2.0);
    if numel(H_subset) > 1
        Z(i,1:2)=X(i,1:2);
%         Z(i,3)=mean(X(i,3));
%         Z(i,4)=mode(X(i,4));
%         Z(i,5)=mode(X(i,5));
%         Z(i,6)=mode(X(i,6));
%         Z(i,7)=mode(X(i,7));
%         Z(i,8)=mode(X(i,8));
%         Z(i,9)=mode(X(i,9));
%         Z(i,10)=mode(X(i,10));
        Z(i,3)=median(H(:,3));
%         X(i,:)=Z(i,:);
    end
    
    
end


end
