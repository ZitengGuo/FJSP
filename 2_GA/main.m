clc;clear
%% 下载数据
% 加工数据包括加工时间，加工机器，机器数，各机器权重，工件数，各工件对应的工序数
load data operation_time operation_machine num_machine machine_weight num_job num_op

%% 基本参数
MAXGEN = 200;               % 最大迭代次数
Ps = 0.8;                   % 选择率
Pc = 0.7;                   % 交叉率
Pm = 0.3;                   % 变异率
sizepop = 200;              % 个体数目
e = 0.5;                    % 目标值权重
trace = zeros(2,MAXGEN);

%% ===========================种群初始化============================
total_op_num=sum(num_op);
chroms=initialization(num_op,num_job,total_op_num,sizepop,operation_machine,operation_time);
[Z,~,~]=fitness(chroms,num_machine,e,num_job,num_op);

%% ============================迭代过程=============================
for gen=1:MAXGEN
    fprintf('当前迭代次数：'),disp(gen)
    % 轮盘赌选择
    chroms_new=selection(chroms,Z,Ps);
    % 交叉操作
    chroms_new=crossover(chroms_new,Pc,total_op_num,num_job,num_op);
    % 变异操作
    chroms_new=mutation(chroms_new,total_op_num,Pm,num_machine,e,num_job,num_op,operation_machine,operation_time);
    % 计算选择交叉变异后个体的适应度
    [Z_new,~,~]=fitness(chroms_new,num_machine,e,num_job,num_op);
    % 根据适应度在原种群和遗传操作后的种群中选出sizepop个更优个体
    [chroms,Z,chrom_best]=update_chroms(chroms,chroms_new,Z,Z_new,sizepop);
    % 记录每代的最优适应度与平均适应度
    trace(1, gen)=Z(1);       
    trace(2, gen)=mean(Z);  
    % 更新全局最优适应度
    if gen==1 || MinVal>trace(1,gen)
        MinVal=trace(1,gen);
    end
end

%% ============================输出结果=============================
%% 输出最优适应度
fprintf('最优适应度：'),disp(MinVal)
%% 描绘解的变化
figure(1)
plot(trace(1,:));
hold on;
plot(trace(2,:),'-.');grid;
legend('解的变化','种群均值的变化');
%% 显示最优解
[Z,machine_weight1,Pvals]=fitness(chrom_best,num_machine,e,num_job,num_op);
Pval1=Pvals{1,1};
figure(2);
for i=1:total_op_num
    mText=chrom_best(total_op_num+i);
    b=chrom_best(i);
    x1=Pval1(1,i);
    x2=Pval1(2,i); 
    y1=mText-0.2;
    y2=mText;
    hold on; 
    fill([x1,x2,x2,x1],[y1,y1,y2,y2],[1-1/b,1/b,b/num_job]);
    text((x1+x2)/2,mText-0.1,num2str(b));
end