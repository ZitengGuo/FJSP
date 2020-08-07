clc;clear
%% 下载数据
% 加工数据包括加工时间，加工机器，机器数，各机器权重，工件数，各工件对应的工序数
load data operation_time operation_machine num_machine machine_weight num_job num_op

%% 基本参数
MAXGEN=200;             % 最大迭代次数
sizepop=201;            % 种群规模
e=0.5;                  % 目标值权重
N_size=30;              % 邻域解数量
S_size=15;              % 共享解数量
G=5;                    % 巡回次数
G1=20;                  % 竞争机制1参数
G2=10;                  % 竞争机制2参数
trace=zeros(2,MAXGEN);
chrom_best=[];

%% ===========================种群初始化============================
total_op_num=sum(num_op);
chroms=initialization(num_op,num_job,total_op_num,sizepop,operation_machine,operation_time);
[Z,~,~,~,~]=fitness(chroms,num_machine,e,num_job,num_op);
% 将最好的解划分为领飞鸟
[Z_leader,ind]=min(Z);
leader=chroms(ind,:); 
% 从chroms中移出领飞鸟，然后划分左右两个跟飞鸟种群
chroms(ind,:)=[];
Z(ind)=[];
sp=(sizepop-1)/2;
lefts=chroms(1:sp,:);
Z_left=Z(1:sp);
rights=chroms(sp+1:end,:);
Z_right=Z(sp+1:end);

%% ============================迭代过程=============================
for gen=1:MAXGEN
    fprintf('当前迭代次数：'),disp(gen)
    %% 巡回阶段
    for i=1:G
        %% 领飞鸟进化
        [leader,Z_leader,share,Z_share]=bird_evolution(leader,Z_leader,[],[],...
            N_size,S_size,total_op_num,num_machine,e,num_job,num_op);
        %% 跟飞鸟进化
        % 初始化左右队列的共享解集
        share_left=share;
        Z_share_left=Z_share;
        share_right=share;
        Z_share_right=Z_share;
        for j=1:sp
            % 左队列
            [lefts(j,:),Z_left(j),share_left,Z_share_left]=bird_evolution(lefts(j,:),Z_left(j),share_left,Z_share_left,...
                N_size-S_size,S_size,total_op_num,num_machine,e,num_job,num_op);
            % 右队列
            [rights(j,:),Z_right(j),share_right,Z_share_right]=bird_evolution(rights(j,:),Z_right(j),share_right,Z_share_right,...
                N_size-S_size,S_size,total_op_num,num_machine,e,num_job,num_op);
        end
        %% 竞争机制2：队间交叉
        % 随机产生G2对位置相同的个体
        ind=randperm(sp,G2);
        [rights(ind,:),Z_right(ind),lefts(ind,:),Z_left(ind)]= crossover(lefts(ind,:),rights(ind,:),...
            Z_left(ind),Z_right(ind),total_op_num,num_machine,e,num_job,num_op);
    end
    %% 竞争机制1：队内竞争
    ind=randperm(sp,G1);
    % 左队列
    [~,ind1]=sort(Z_left(ind));
    lefts(ind,:)=lefts(ind(ind1),:);
    Z_left(ind)=Z_left(ind(ind1));
    % 右队列
    [~,ind2]=sort(Z_right(ind));
    rights(ind,:)=rights(ind(ind2),:);
    Z_right(ind)=Z_right(ind(ind2));
    %% 领飞鸟替换
    if rand<0.5
        % 选择左队列首只跟飞鸟
        [leader,Z_leader,lefts,Z_left]=update_leader(leader,Z_leader,lefts,Z_left);
    else
        % 选择右队列首只跟飞鸟
        [leader,Z_leader,rights,Z_right]=update_leader(leader,Z_leader,rights,Z_right);
    end
    %% 记录相关数据
    % 记录每代的最优适应度与平均适应度
    Z=[Z_leader,Z_left,Z_right];
    [val,ind]=min(Z);
    trace(1,gen)=val;
    trace(2,gen)=mean(Z);
    % 更新全局最优适应度
    if gen==1 || MinVal>trace(1,gen)
        MinVal=trace(1,gen);
    end
end

%% ============================输出结果=============================
%% 输出最优适应度
fprintf('最优适应度：'),disp(MinVal)
%% 绘制最优适应度与平均适应度的迭代曲线图
figure(1)
plot(trace(1,:));
hold on;
plot(trace(2,:),'-.');grid;
legend('解的变化','种群均值的变化');
%% 绘制全局最优解的甘特图
[Z,~,~,machine_time,Pvals]=fitness(leader,num_machine,e,num_job,num_op);
Pval1=Pvals{1,1};
figure(2);
for i=1:total_op_num
    mText=leader(total_op_num+i);
    b=leader(i);
    x1=Pval1(1,i);
    x2=Pval1(2,i); 
    y1=mText-0.2;
    y2=mText;
    hold on; 
    fill([x1,x2,x2,x1],[y1,y1,y2,y2],[1-1/b,1/b,b/num_job]);
    text((x1+x2)/2,mText-0.1,num2str(b));
end