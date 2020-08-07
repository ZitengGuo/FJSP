function chroms = mutation(chroms,total_op_num,Pm,num_machine,e,num_job,num_op,operation_machine,operation_time)
%% 针对工序码，设计变邻域变异
for i=1:size(chroms,1)
    if Pm>rand
        chrom=chroms(i,:);
        chrom_best=[];
        Z_best=[];
        tau_min=2;
        tau_max=4;
        for tau=tau_min:tau_max
            % 选择tau个不同位置、代表不同工件的基因
            job=randperm(num_job,tau);  % 选择的工件
            job_ind=zeros(1,tau);  % 工件的位置
            for j=1:tau
                temp=find(chrom(1:total_op_num)==job(j));
                jj=randperm(num_op(job(tau)),1);
                job_ind(j)=temp(jj);
            end
            % tau个位置的全排列（索引）
            ind=perms(1:tau);
            % 评价当前邻域中最优个体
            Z=zeros(1,length(ind));
            chrom_neigh=[];  % 邻域解
            for k=1:length(ind)
                chrom1=chrom;
                chrom1(job_ind)=job(ind(k,:));
                chrom1(total_op_num+job_ind)=chrom(total_op_num+job_ind(ind(k,:)));
                chrom1(total_op_num*2+job_ind)=chrom(total_op_num*2+job_ind(ind(k,:)));
                chrom_neigh=[chrom_neigh;chrom1];
                [Z(k),~,~] = fitness(chrom1,num_machine,e,num_job,num_op);
            end
            % 更新当前解
            [val,ii]=min(Z);
            Z_best=[Z_best,val];
            chrom_best=[chrom_best;chrom_neigh(ii,:)];
        end
        [~,ii]=min(Z_best);
        chroms(i,:)=chrom_best(ii,:);
    end
end
%% 针对机器码，设计单点变异,即替换可选机器
for i=1:size(chroms,1)
    if Pm>rand
        chrom=chroms(i,:);
        ind=randperm(total_op_num,1);  % 随机选择的位置
        job=chrom(ind);  % 所选工件
        job_ind=find(chrom(1:total_op_num)==job);
        op=find(job_ind==ind);  % 所选位置的工序
        machines=operation_machine{job}{op};  % 对应的加工机器集
        times=operation_time{job}{op}; %  对应的加工时间集
        if length(machines)>1
            ii=randperm(length(machines)-1,1);
            ind1=find(machines==chrom(total_op_num+ind));
            machines(ind1)=[];
            times(ind1)=[];
            chrom(total_op_num+ind)=machines(ii);
            chrom(total_op_num*2+ind)=times(ii);
        end
        chroms(i,:)=chrom;
    end
end