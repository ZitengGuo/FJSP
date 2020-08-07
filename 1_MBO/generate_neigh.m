function chroms=generate_neigh(chrom,neigh_size,total_op_num,num_job,num_op)
chroms=[];
flag=0;  % 用于判断是否跳出主循环
tau=randi([2,4]);  % 随机产生tau值
for i=1:neigh_size
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
    for j=1:length(ind)
        chrom1=chrom;
        chrom1(job_ind)=job(ind(j,:));
        chrom1(total_op_num+job_ind)=chrom(total_op_num+job_ind(ind(j,:)));
        chrom1(total_op_num*2+job_ind)=chrom(total_op_num*2+job_ind(ind(j,:)));
        chroms=[chroms;chrom1];
        % 只需要产生neigh_size个邻域解
        if size(chroms,1)>=neigh_size
            flag=1;
            break;
        end
    end
    % 达到所需邻域解数量则跳出主循环
    if flag==1
        break;
    end
end