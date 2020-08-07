function [Z,machine_weight,pvals] = fitness(chroms,num_machine,e,num_job,num_op)
sizepop=size(chroms,1);
pvals=cell(1,sizepop);
Z1=zeros(1,sizepop);
Z2=Z1;
total_op_num=sum(num_op);  % 总工序数
for k=1:sizepop
    chrom=chroms(k,:);
    machine=zeros(1,num_machine);  % 记录各机器变化时间
    job=zeros(1,num_job);  % 记录各工件变化时间
    machine_time=zeros(1,num_machine);  % 计算各机器的实际加工时间
    pval=zeros(2,total_op_num);  % 记录各工序开始和结束时间
    for i=1:total_op_num
        % 机器时间大于工件时间
        if machine(chrom(total_op_num+i))>=job(chrom(i))
            pval(1,i)=machine(chrom(total_op_num+i));  % 记录工件开始时间
            machine(chrom(total_op_num+i))=machine(chrom(total_op_num+i))+chrom(total_op_num*2+i);
            job(chrom(i))=machine(chrom(total_op_num+i));
            pval(2,i)=machine(chrom(total_op_num+i));  % 记录工件结束时间
            % 机器时间小于工件时间
        else
            pval(1,i)=job(chrom(i));
            job(chrom(i))=job(chrom(i))+chrom(total_op_num*2+i);
            machine(chrom(total_op_num+i))=job(chrom(i));
            pval(2,i)=job(chrom(i));
        end
        machine_time(chrom(total_op_num+i))=machine_time(chrom(total_op_num+i))+chrom(total_op_num*2+i);
    end
    Z1(k)=max(machine);  % 最大机器时间值，对应makespan
    % machine_weight=machine_time/sum(machine_time);  % 计算各机器的负荷
    machine_weight=machine_time;
    Z2(k)=max(machine_weight)-min(machine_weight);
    pvals{k}=pval;
end
% min_makespan=min(Z1);%所有染色体的makespan最优值
% max_makespan=max(Z1);
% min_weight=min(Z2);%负载最优值
% max_weight=max(Z2);
% Z=e*((Z1-min_makespan)./(max_makespan-min_makespan))+(1-e)*((Z2-min_weight)./(max_weight-min_weight));%计算适应度
Z=e*Z1+(1-e)*Z2;

