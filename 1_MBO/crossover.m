%% 候鸟算法中的交叉函数与遗传算法的不同
%% 候鸟算法输入两个染色体种群，分别来自左右队列
%--------------------------------------------------------------------------
function [lefts,Z_left,rights,Z_right]= crossover(lefts,rights,Z_left,Z_right,total_op_num,num_machine,e,num_job,num_op)
chroms1=lefts;
chroms2=rights;
for i=1:size(chroms1,1)
    %% 面向工序码的交叉操作
    % 父代染色体
    parent1=lefts(i,:);
    parent2=rights(i,:);
    Job=randperm(num_job);
    % 将工件随机分成两个集合
    J1=Job(1:round(num_job/2));
    J2=Job(length(J1)+1:end);
    op_p1=[];
    op_p2=[];
    for j=1:length(J2)
        %找出父代中J2片段对应的位置
        op_p1=[op_p1,find(parent1(1:total_op_num)==J2(j))];
        op_p2=[op_p2,find(parent2(1:total_op_num)==J2(j))];
    end
    op_s1=sort(op_p1);
    op_s2=sort(op_p2);
    % 子代1交换J2片段的基因，机器码对应位置的基因，工时码对应位置的基因
    chroms1(i,op_s1)=parent2(op_s2);
    chroms1(i,total_op_num+op_s1)=parent2(total_op_num+op_s2);
    chroms1(i,total_op_num*2+op_s1)=parent2(total_op_num*2+op_s2);
    % 子代2同理
    chroms2(i,op_s2)=parent1(op_s1);
    chroms2(i,total_op_num+op_s2)=parent1(total_op_num+op_s1);
    chroms2(i,total_op_num*2+op_s2)=parent1(total_op_num*2+op_s1);
    
    %% 面向机器码的交叉操作
    parent1=chroms1(i,:);
    parent2=chroms2(i,:);
    % 随机产生与染色体长度相等的0,1序列
    rand0_1=randi([0,1],1,total_op_num);
    for n=1:num_job
        ind_0=find(rand0_1(num_op(n)*(n-1)+1:num_op(n)*n)==0);
        if ~isempty(ind_0)
            ind1=find(parent1(1:total_op_num)==n);
            ind2=find(parent2(1:total_op_num)==n);
            chroms1(i,total_op_num+ind1(ind_0))=parent2(total_op_num+ind2(ind_0));
            chroms1(i,total_op_num*2+ind1(ind_0))=parent2(total_op_num*2+ind2(ind_0));
            chroms2(i,total_op_num+ind2(ind_0))=parent1(total_op_num+ind1(ind_0));
            chroms2(i,total_op_num*2+ind2(ind_0))=parent1(total_op_num*2+ind1(ind_0));
        end
    end
end
%% 判断个体是否可以更新
[Z1,~,~,~,~]=fitness(chroms1,num_machine,e,num_job,num_op);
[Z2,~,~,~,~]=fitness(chroms2,num_machine,e,num_job,num_op);
lefts(Z1<Z_left,:)=chroms1(Z1<Z_left,:);
Z_left(Z1<Z_left)=Z1(Z1<Z_left);
rights(Z2<Z_right,:)=chroms2(Z2<Z_right,:);
Z_right(Z2<Z_right)=Z2(Z2<Z_right);