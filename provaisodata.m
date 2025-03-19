clear all
close all
clc
load('\dades.mat')
%������������X��Y
% X  1x267 double array
% Y  1x267 double array

ON=15;  %һ���������������� 
OC=7;  %����������֮�����С����.
OS=12;  %һ������������������ֲ��ı�׼��.
k=8;   %�����õ��ľ�����.
L=2;   %һ�ε��������п��Ժϲ��ľ������ĵ�������.
I=60;  %��������Ĵ���.
NO=1;  %����Ĳ������Զ��ش��γɲ�û���κβ�����Ҫ��
min=50;%��ÿһ�����ĵ���С����Ӧ������κ�ʱ������Ҫɾ������߶����ۡ�
 
[centro, Xcluster, Ycluster, A, clustering]=isodata(X, Y, k, L, I, ON, OC, OS, NO, min);
clc;
fprintf('�۳������Ŀ: %d',A);

 
% ��������ɫ
colr=zeros(A,3);
for i=1:A
    colr(i,:)=rand(1,3);
end;

%��ʾ�������Ϣ
figure;
hold on;
for i=1:A,
    n=find(clustering==i);
    p=plot(X(n), Y(n),'.');
    set(p,'Color',colr(i,:));
    title(A);
end;
 
clc;
fprintf('�۳������Ŀ: %d',A);

%ɾ����ʱ������
%clear n;clear i;clear p;clear colr;clear NO;