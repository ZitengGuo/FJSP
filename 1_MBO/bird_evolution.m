function [chrom,Z_chrom,share,Z_share]=bird_evolution(chrom,Z_chrom,share,Z_share,neigh_size,share_size,total_op_num,num_machine,e,num_job,num_op)
% 产生邻域解并计算其适应度
neigh=generate_neigh(chrom,neigh_size,total_op_num,num_job,num_op);
[Z_neigh,~,~,~,~] = fitness(neigh,num_machine,e,num_job,num_op);
% 产生的邻域解、共享解与当前解合并
neigh=[neigh;share;chrom];
Z_neigh=[Z_neigh,Z_share,Z_chrom];
% 选择最优解，更新当前个体
[~,ind]=sort(Z_neigh);
chrom=neigh(ind(1),:);
Z_chrom=Z_neigh(ind(1));
% 更新共享邻域解集
share=neigh(2:share_size+1,:);
Z_share=Z_neigh(2:share_size+1);