function loading_order = merit_order(X,gendata)
% conventional merit order ตามความถูกของ fuel consumption
% col1 คือ capacity, col2 คือ FOR, col3 คือ operating, col4 คือ fix O&M

tmp = X>0;
tmp_gen = X(tmp);
gendata_tmp = gendata(tmp,:);
gendata_comp =[];

for i=1:size(gendata_tmp,1)
    gen_tmp = [];
    numgen = tmp_gen(i); %
    for j=1:numgen
        gen_tmp = vertcat(gen_tmp,gendata_tmp(i,:));
    end
    gendata_comp = vertcat(gendata_comp,gen_tmp);
end

[a,IX] = sort(gendata_comp,'ascend');      % นำมาเรียงลำดับจากค่าเชื้อเพลิงถูกจากน้อยไปมาก
loading_order=gendata_comp(IX(:,4),:);
end
