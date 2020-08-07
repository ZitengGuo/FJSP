%% 种群初始化
function chroms=initialization(num_op,num_job,total_op_num,sizepop,operation_machine,operation_time)
% 随机生成与工序数对应的工序
a=1;
for j=1:num_job
    b=num_op(j);
    temp(1,a:(a+b-1))=j;
    a=a+b;
end   
% 染色体第一层为工序号，第二层为对应加工机器，第三层为对应加工时间
chroms=zeros(sizepop,total_op_num*3);
for i=1:sizepop
    % 随机生成工序码
    rt=randperm(total_op_num);
    chroms(i,1:total_op_num)=temp(rt);
    % 随机生成机器码，并记录对应的工时
    for j=1:num_job
        job=find(chroms(i,1:total_op_num)==j);
        for k=1:length(job)
            machines=operation_machine{j}{k};  % 工件j工序k可以选择的机器集
            times=operation_time{j}{k};  % 相应的工时集
            ind=randperm(length(machines),1);  % 随机选择一个机器
            chroms(i,total_op_num+job(k))=machines(ind);
            chroms(i,total_op_num*2+job(k))=times(ind);
        end
    end
end