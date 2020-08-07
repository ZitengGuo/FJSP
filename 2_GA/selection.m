function [chrom1] = selection(chrom,Z,Ps)

size_chrom=round(size(chrom,1)*Ps);  % 被选择进行交叉变异的染色体数
Z=Z+0.00000001;
Z1=1./Z; % 取倒数
P=Z1./sum(Z1);
Q=cumsum(P,2);
chrom1=chrom(1:size_chrom,:);
for i=1:size_chrom
    temp=find(Q>rand,1);
    chrom1(i,:)=chrom(temp,:);
end