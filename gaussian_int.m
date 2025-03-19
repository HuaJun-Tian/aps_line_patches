function [Z1,Z2] = gaussian_int(xy,phase,mn,interval,sigma)
% filter and interp based on Gaussia
%
%
% input：
% xy: coordinate n*2
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
[n1,n2] =size(xy);
Z1=zeros(n1,1)*nan;
Z2=zeros(n1,1)*nan;
p=parpool(24);
% mask coffecient
parfor i=1:n1
    %     save test
%     clear pha_single dis G H_subset;
    H_subset=find( xy(:,1)>=xy(i,1)-m*interval & xy(:,1)<=xy(i,1)+m*interval & xy(:,2)>=xy(i,2)-n*interval & xy(:,2)<=xy(i,2)+2*interval);
    H=xy(H_subset,:);
    pha_single=phase(H_subset,:);
    H_demean=H-mean(H);
    dis=sum(H_demean.^2,2);
    G=exp((-1.0)*dis/(2*sigma*sigma))/(2*pi*(sigma^2));
    G=G/sum(G(:));
    if numel(H_subset) > 1

        Z1(i,1)=sum(pha_single(:,1).*G);
        Z2(i,1)=sum(pha_single(:,2).*G);
        %         X(i,:)=Z(i,:);
    end
        


end
delete(p);
% save G_ G -mat
end