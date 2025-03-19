function [centers,polys,class_ori,index_class]=ISODATA(x,K,theta_N,theta_S,theta_c,L,I,others_val,lonlat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%input parameters%%%%%%
% x : data including index
% K : 预期的聚类中心数
% theta_N : 每一聚类中心中最少的样本数，少于此数就不作为一个独立的聚类
% theta_S ：一个聚类中样本距离分布的标准差
% theta_c : 两聚类中心之间的最小距离，如小于此数，两个聚类进行合并
% L : 在一次迭代运算中可以和并的聚类中心的最多对数
% I ：迭代运算的次数序号
% Jun Hua, 2023.5.11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% step1
n = size(x,1);
%normalize
ori=double(x); %original data
x=double(x);
% ori=[(1:size(x,1))',ori];
% coor=[x(:,1);x(:,2)];
% coor_nor=normalize(coor,'norm');
% dis_temp=normalize(x(:,3),'norm');
% clear x;
% x=[coor_nor(1:n),coor_nor(n+1:2*n),dis_temp];
% clear coor_nor coor dis_temp;
x=Whiten_process_self(x(:,2:end));
% return
x=[ori(:,1),x];
index_sample=[(1:size(x,1))']; % index every sample
N_c = K;
mean_pat = cell(K,1);
num_temp=floor(n/K);
for i=1:K
    mean_pat{i} = x(i+(i-1)*num_temp,2:end);
end
ite = 1;
while ite<I
    flag = 1;
    while flag
        %% step2
        index_class=cell(size(mean_pat));
        class = cell(size(mean_pat));
        %     for i=1:n
        %         num = Belong2(x(i,:),mean_pat);
        %         class{num} =  [class{num};x(i,:)];
        %     end
        num = Belong(x(:,2:end),mean_pat);
        %     class{num} =  [x(i,:)];
        for i=1:N_c
            idx_temp=find(num==i);
            if (~isnan(idx_temp))
                class{i} =  [x(idx_temp,:)];
                %                 keyboard
                index_class{i}=[index_sample(idx_temp,:),others_val(idx_temp,:)];
            end

            clear idx_temp;
        end
        %         keyboard;
        clear num;
        %% step3
        flag_N=1;
        while flag_N
            for i=1:N_c
                %                     keyboard;
                size_i = size(class{i},1);
                if size_i<theta_N
                    %                 keyboard;
                    class_i = class{i};
                    mean_pat = DeleteRow(mean_pat,i);
                    class = DeleteRow(class,i);
                    index_class=DeleteRow(index_class,i);
                    N_c = N_c-1;
                    num = Belong(class_i(:,2:end),mean_pat);
                    for iii=1:N_c
                        idx_temp=find(num==iii);
                        if (~isnan(idx_temp))
                            class{iii} =  [class{iii};x(idx_temp,:)];
                            index_class{iii}=[index_class{iii};[index_sample(idx_temp,:),others_val(idx_temp,:)]];
                        end
                        clear idx_temp;
                    end
                    clear num;
                    flag_N=1;
                    break;
                else
                    flag_N=0;

                    %         clear idx_temp;
                    %     end
                    %           for j=1:size_i
                    %             class_ij = class_i(j,:);%the j'th row of class{i}
                    %             num = Belong2(class_ij,mean_pat);
                    %             class{num} = [class{num};class_ij];
                    %           end
                end
            end

        end

        %% step4
        for i=1:N_c
            if ~isempty(mean_pat{i})
                mean_pat{i} = sum(class{i}(:,2:end))./size(class{i},1);
            end
        end
        %% step5
        Dis = zeros(N_c,1);
        for i=1:N_c
            if ~isempty(class{i})
                [N_i,~] =size(class{i});
                %             tmp = bsxfun(@minus,class{i},mean_pat{i});
                tmp=class{i}(:,2:end)-repmat(mean_pat{i},N_i,1);
                Dis(i)=sum(sqrt(sum(tmp.^2,2)))/N_i;
                %             Dis(i) = sum(arrayfun(@(x)norm(tmp(x,:)),1:N_i))/N_i;
                clear tmp;
            end
        end
        %% step6
        D = 0;
        for i=1:N_c
            if ~isempty(class{i})
                N_i =size(class{i},1);
                D = D + N_i*Dis(i);
            end
        end
        D = D/n;
        %% step7
        flag = 0;
        if ite == I
            theta_c = 0;
            flag = 0;
        elseif ~(N_c > K/2)
            flag = 1;
        elseif mod(ite,2)==0 || ~(N_c<2*K)
            flag = 0;
        elseif mod(ite,2)==1
            flag = 1;
        end
        %% 分裂处理
        %% step8
        if flag
            flag = 0;
            delta = cell(N_c,1);
            for i=1:N_c
                if ~isempty(class{i})
                    N_i =size(class{i},1);
                    tmp=class{i}(:,2:end)-repmat(mean_pat{i},N_i,1);

                    %                  tmp = bsxfun(@minus,class{i},mean_pat{i});
                    delta{i} = arrayfun(@(x)norm(tmp(:,x)),1:size(tmp,2))/N_i;
                end
            end

            %% step9
            delta_max = cell(N_c,1);
            for i=1:N_c
                if ~isempty(class{i})
                    max_i = max(delta{i});
                    sub = find(delta{i}==max_i,1);
                    delta_max{i} = [max_i,sub];
                end
            end
            %% step10
            for i=1:N_c
                if delta_max{i}(1) > theta_S
                    N_i =size(class{i},1);
                    con1 = (Dis(i)>D && N_i>2*(theta_N + 1));
                    con2 = ~(N_c>K/2);
                    if con1 || con2
                        %%%%这里分裂%%%%%
                        flag = 1;%一旦发生分裂，那么分裂一次后就返回第二步；若没发生分裂，则直接进入合并处理步
                        lamda = 0.5;
                        max_sub = delta_max{i}(2);
                        mean_pat{i}(max_sub) = mean_pat{i}(max_sub) + lamda * delta_max{i}(1);
                        addOneMean =  mean_pat{i};
                        addOneMean(max_sub) = addOneMean(max_sub) - lamda * delta_max{i}(1);
                        mean_pat = [mean_pat;addOneMean];
                        N_c = N_c+1;
                        break;
                    end
                end
            end

        end

    end
    %% 合并处理
    if L
        %% step11
        Distance = zeros(N_c,N_c);
        for i=1:N_c-1
            for j=i:N_c
                Distance(i,j) = norm(mean_pat{i}-mean_pat{j});
            end
        end
        %         keyboard
        %% step12
        index = find(Distance<=theta_c & Distance>0);
%         keyboard
        keepIndex = [Distance(index),index];
        try
            [~, index] = sort(keepIndex(:,1));
        catch
            index
            Distance
            ite=ite+1;
            continue;
            error('please adjust parameter')
        end
        if size(index,1) > L
            index = index(1:L,:);
        end
        %% step13
        iscombin=zeros(N_c,1);
        num_class=N_c;
        if size(index,1) ~= 0
            combin_num=0;
            %             keyboard
            for id=1:size(index,1)
                [m_i m_j]= seq2idx(keepIndex(index(id),2),N_c);
                % flag mark
                if (~iscombin(m_i) && ~iscombin(m_j))
                    %%%%%这里合并%%%%%
                    combin_num=combin_num+1;
                    dele_index(combin_num)=m_j;
                    N_mi = size(class{m_i},1);
                    N_mj = size(class{m_j},1);
                    mean_pat{m_i}=(N_mi*mean_pat{m_i} + N_mj*mean_pat{m_j})/(N_mi+N_mj);
                    %                     mean_pat_new{m_i}=[m_i m_j mean_new];
                    class{m_i} = [class{m_i};class{m_j}];
                    index_class{m_i}=[index_class{m_i};index_class{m_j}];
                    iscombin(m_i)=1;
                    iscombin(m_j)=1;
                    num_class=num_class-1;
                end
                %                 mean_pat{m_i} = (N_mi*mean_pat{m_i} + N_mj*mean_pat{m_j})/(N_mi+N_mj);
                %                 mean_pat = DeleteRow(mean_pat,m_j);
                %                 class{m_i} = [class{m_i};class{m_j}];
                %                 class = DeleteRow(class,m_j);
            end
            %             keyboard;
            mean_pat(dele_index)=[];
            class(dele_index)=[];
            index_class(dele_index)=[];
            N_c=num_class;
            clear dele_index keepIndex  index num_class
            %             new_num=numel(mean_pat_new);
            % %             combin_num=0;
            %             for i=1:new_num
            %                 c
            %             end
        end
    end
    %% step14
    %     figure();hold on;
    %     for i=1:N_c
    %         scatter(class{i}(:,2),class{i}(:,3),[],[i/(N_c+2),1-i/(N_c+2),i/(N_c+2)]);colorbar;title(['iteration: ',num2str(ite)]);%colormap('jet');
    %
    %     end
    %     hold off;
    %     set(gcf,'PaperSize',[5,7],'PaperUnits','centimeters');box on;
    %     saveas(gcf,['aps_patches/iteration_',num2str(ite)],'fig');
    %     print(['aps_patches/iteration_',num2str(ite)],gcf,'-djpeg','-r150')
    ite=ite+1;
end
% figure();hold on;
temp_total=[];
temp_mean=[];
string_patch=[];
for i=1:N_c
    %     keyboard
    temp_total=[temp_total;[class{i}(:,2),class{i}(:,3),repmat(mean_pat{i}(3),size(class{i},1),1)]];
    temp_mean=[temp_mean;[mean_pat{i}(1),mean_pat{i}(2)]];
    %     keyboard
    string_patch=[string_patch;["patch"+num2str(i,'%4d')]];
    %     scatter(class{i}(:,2),class{i}(:,3),[],repmat(mean_pat{i}(3),size(class{i},1),1));colorbar;title(['iteration: ',num2str(ite)]);%colormap('jet');
    %     text(mean_pat{i}(1),mean_pat{i}(2),['patch',num2str(i)],'Color','red');
    %     keyboard
end
% hold off;
set(gcf,'PaperSize',[5,7],'PaperUnits','centimeters');box on;
saveas(gcf,['aps_patches/iteration_',num2str(ite)],'fig');
print(['aps_patches/iteration_',num2str(ite)],gcf,'-djpeg','-r150')
figure();hold on;
scatter(temp_total(:,1),temp_total(:,2),[],temp_total(:,3));colorbar;title(['iteration: ',num2str(ite),' The final patches']);colormap('jet');
scatter(temp_mean(:,1),temp_mean(:,2),"filled",'^','Color','red');
text(temp_mean(:,1),temp_mean(:,2),string_patch,'Color','red');
box on;
hold off;
set(gcf,'PaperSize',[5,7],'PaperUnits','centimeters');box on;
saveas(gcf,['aps_patches/iteration_',num2str(ite),'_The final patches'],'fig');
print(['aps_patches/iteration_',num2str(ite),'_The final patches'],gcf,'-djpeg','-r150');
clear temp_total temp_mean string_patch;
%    keyboard
class_ori=class;
centers=[];
polys={};
number_total=0;
temp_total=[];
temp_mean=[];
temp_total_ll=[];
temp_mean_ll=[];
string_patch=[];
index_temp=[];
% figure()
for  i=1:N_c
    number_total=size(class_ori{i},1)+number_total;

    ix_ori=index_class{i}(:,1);
    class_ori{i}=ori(ix_ori,:);
    K=convhull(class_ori{i}(:,2),class_ori{i}(:,3));
    centers=[centers;mean(class_ori{i}(:,2:3))];
    polys=[polys;class_ori{i}(K,:)];

    %     plot(class_ori{i}(K,2),class_ori{i}(K,3),'color',[0.8 0.8 0.8]);colorbar;hold on;
    %     patch(class_ori{i}(K,2),class_ori{i}(K,3),[i/(N_c+2),1-i/(N_c+2),i/(N_c+2)]);colorbar;
    %     temp_total=class_ori{i};
    %     keyboard
    temp_total=[temp_total;[class_ori{i}(:,2),class_ori{i}(:,3),repmat(mean(class_ori{i}(:,4)),size(class_ori{i},1),1)]];
    temp_mean=[temp_mean;centers(i,1:2)];
    temp_total_ll=[temp_total_ll;[lonlat(ix_ori,1),lonlat(ix_ori,2),repmat(mean(class_ori{i}(:,4)),size(class_ori{i},1),1)]];
    temp_mean_ll=[temp_mean_ll;mean(lonlat(ix_ori,:))];
    index_temp=[index_temp;index_class{i}(:,1)];
    %     string_patch=[string_patch;['patch',num2str(i,'%4d')]];
    clear ix_ori K;
end
%% judge if there are classfied correctly
% xy_inique=[unique(temp_total(:,1)),unique(temp_total(:,2))];
% if (size(xy_inique,1)~=number_total || min(index_temp)~=1 || max(index_temp)~=number_total || size(index_temp,1)~=number_total)
%
%     error('the program is something wrong !');
% else
%     fprintf('Yes, the patches is  classfied correctly\n');
% end
% hold off;
% keyboard
figure();hold on;
scatter(temp_total(:,1),temp_total(:,2),[],temp_total(:,3));colorbar;title(['iteration: ',num2str(ite),' final patches']);colormap('jet');
scatter(temp_mean(:,1),temp_mean(:,2),"filled",'^','Color','red');
% text(temp_mean(:,1),temp_mean(:,2),string_patch,'Color','red');
box on;
hold off;

fprintf('toall pixel:%d\n',number_total);
fprintf('maximun dis: %f, minmum dis: %f\n',max(max(Distance)),min(min(Distance(Distance>0))));
fprintf('Average dis: %f\n',D);
fprintf('maximum std:\n');
delta_max
saveas(gcf,'aps_patches/patches.png')
saveas(gcf,'aps_patches/patches.fig')
figure();hold on;
scatter(temp_total_ll(:,1),temp_total_ll(:,2),[],temp_total_ll(:,3));colorbar;title(['iteration: ',num2str(ite),' final patches']);colormap('jet');
scatter(temp_mean_ll(:,1),temp_mean_ll(:,2),"filled",'^','Color','red');

hold off;
clear temp_tota temp_mean string_patch xy_inique index_temp;
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function number = Belong(x_i,means_pat)
%     INF = 10000;
%     min_ini = inf;
[m_row,m_col]=size(x_i);
kk = numel(means_pat);
%     number = 1;
dis=zeros(m_row,kk)+inf;
for i=1:kk
    if ~isempty(means_pat{i}) & size(means_pat{i},2)==m_col
        %             dis(:,i)=norm(x_i-repmat(means_pat{i},m_row,1));
%                 keyboard;
        dis(:,i)=sqrt(sum((x_i - repmat(means_pat{i},m_row,1)).^2,2));
        a_temp=sqrt(sum((x_i - repmat(means_pat{i},m_row,1)).^2,2));
        %             if norm(x_i - means_pat{i}) < min
        %                 min = norm(x_i - means_pat{i});
        %                 number = i;
        %             end
    end
end
%     keyboard;
[min_dis,number]=min(dis,[],2,'omitnan');
end


function A_del = DeleteRow(A,r)
%     n = size(A,1);
A_del=A;
A_del(r)=[];
%     if r == 1
%         A_del = A(2:n,:);
%     elseif r == n
%         A_del = A(1:n-1,:);
%     else
%         A_del = [A(1:r-1,:);A(r+1:n,:)];
%     end
end


function [row col] = seq2idx(id,n)
if mod(id,n)==0
    row = n;
    col = id/n;
else
    row = mod(1.0*id,n);
    col = ceil(1.0*id/n);
end
end
