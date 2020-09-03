function [newnest,new_fitness,new_ZVal]=CalcFitness(current_nest,levy_nest,fitness,ZVal)
% Evaluating all new solutions
% โดย format คำตอบจอว nest จะอยู่ในรูปของ virtual mapping

new_fitness = fitness;
new_ZVal = ZVal;
newnest = current_nest;

[fnew,Znew]=CalObj_rev01(levy_nest,fitness,ZVal);
for i=1:size(current_nest,1),
    if fnew(i)<=fitness(i),
       new_fitness(i)=fnew(i);
       new_ZVal(i)=Znew(i);
       newnest(i,:)=levy_nest(i,:);
    end
end
 
end

% [fitness,ZVal]= CalObj(feas_row,fitness,ZVal);