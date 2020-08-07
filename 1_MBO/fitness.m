function [Z,makespan,machine_load,machine_weight,pvals] = fitness(chroms,num_machine,e,num_job,num_op)
sizepop=size(chroms,1);
pvals=cell(1,sizepop);
makespan=zeros(1,sizepop);
machine_load=makespan;
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
    makespan(k)=max(machine);
    machine_weight=machine_time;
    machine_load(k)=max(machine_weight)-min(machine_weight);
    pvals{k}=pval;
end
Z=e*makespan+(1-e)*machine_load;