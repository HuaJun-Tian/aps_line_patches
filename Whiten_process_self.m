function whitening=Whiten_process_self(x)
%% Whiten_process.m
% author Hua Jun. 2023/5/12
x=double(x);
[row,~]=size(x);
meanx=mean(x); %求矩阵每列的平均值
% tempx=repmat(meanx,row,1); %堆叠矩阵
x=x-repmat(meanx,row,1);   
c=x'*x/row; %求矩阵x的协方差矩阵
% keyboard
[F,V]=eigs(c); %F的列向量为对应特征向量,V为最大特征值对角阵
epsilion=1e-8;
d=diag(V)+epsilion;
v=diag(d.^(-0.5));
% Y_whiten=x*(F*v);
Y_whiten=F*v*F'*x';
whitening=Y_whiten'; 
cor=whitening'*whitening/row;
%% datainput.m

figure()
scatter(whitening(:,1),whitening(:,2),[],whitening(:,3)); colorbar;colormap('jet');title('whitening');

saveas(gcf,'aps_patches/whitening.png')

fprintf('convariance after whitening: \n');
disp(cor);
close all
% keyboard
% reference_std=F*v*F*[0 0 1]
% fprintf('reference std is %f\n',reference_std);
end
%% selfpca.m
function pca=selfpca(x)
[row,~]=size(x);
c=cov(x); %求矩阵x的协方差矩阵
[F,~]=eigs(c); %F的列向量为对应特征向量
meanx=mean(x); %求矩阵每列的平均值
tempx=repmat(meanx,row,1); %堆叠矩阵
Y_pca=(x-tempx)*F;
pca=Y_pca(:,1:2); %取第1、2列
end

