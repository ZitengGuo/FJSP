function [ chroms_new] = crossover(chroms,Pc,total_op_num,num_job,num_op)
size_chrom=size(chroms,1);  % 染色体数
chroms_new=chroms;
%% 面向工序码的交叉操作
for i=1:2:size_chrom-1 
    if Pc>rand
        % 父代染色体
        parent1=chroms(i,:);
        parent2=chroms(i+1,:);
        Job=randperm(num_job);
        % 将工件随机分成两个集合
        J1=Job(1:round(num_job/2));
        J2=Job(length(J1)+1:end);
        % 子代染色体
        child1=parent1;
        child2=parent2;
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
        child1(op_s1)=parent2(op_s2);
        child1(total_op_num+op_s1)=parent2(total_op_num+op_s2);
        child1(total_op_num*2+op_s1)=parent2(total_op_num*2+op_s2);
        % 子代2同理
        child2(op_s2)=parent1(op_s1);
        child2(total_op_num+op_s2)=parent1(total_op_num+op_s1);
        child2(total_op_num*2+op_s2)=parent1(total_op_num*2+op_s1);
        chroms_new(i,:)=child1;
        chroms_new(i+1,:)=child2;
    end
end
%% 面向机器码的交叉操作
for k=1:2:size_chrom-1
    if Pc>rand
        parent1=chroms_new(k,:);
        parent2=chroms_new(k+1,:);
        child1=parent1;
        child2=parent2;
        % 随机产生与染色体长度相等的0,1序列
        rand0_1=randi([0,1],1,total_op_num);
        for n=1:num_job
            ind_0=find(rand0_1(num_op(n)*(n-1)+1:num_op(n)*n)==0);
            if ~isempty(ind_0)
                temp1=find(parent1(1:total_op_num)==n);
                temp2=find(parent2(1:total_op_num)==n);
                child1(total_op_num+temp1(ind_0))=parent2(total_op_num+temp2(ind_0));
                child2(total_op_num+temp2(ind_0))=parent1(total_op_num+temp1(ind_0));
                child1(total_op_num*2+temp1(ind_0))=parent2(total_op_num*2+temp2(ind_0));
                child2(total_op_num*2+temp2(ind_0))=parent1(total_op_num*2+temp1(ind_0));
            end
        end
        chroms_new(k,:)=child1;
        chroms_new(k+1,:)=child2;
    end
end