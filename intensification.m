function newnest=intensification(nest,Kbest)
% Levy flights
global LowerBound UpperBound T VMP CumCap CumOilCap CumLNGCap CumCoalCap CumNuclearCap Mmax Mmin

bestnest= nest(Kbest,:);
nest(Kbest,:)=[];
%CumCap(Kbest,:)=[];
%CumOilCap(Kbest,:)=[]; 
%CumLNGCap(Kbest,:)=[];
%CumCoalCap(Kbest,:)=[];
%CumNuclearCap(Kbest,:)=[];

% Levy exponent and coefficient
% For details, see equation (2.21), Page 16 (chapter 2) of the book
% X. S. Yang, Nature-Inspired Metaheuristic Algorithms, 2nd Edition, Luniver Press, (2010).

beta=3/2;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);

n=size(nest,1);
CumCap = zeros(n,1);
CumOilCap = zeros(n,1);
CumLNGCap = zeros(n,1);
CumCoalCap = zeros(n,1);
CumNuclearCap = zeros(n,1);
%
minStep = LowerBound(1)-5450;
maxStep = UpperBound(1)-5450;

feas_trans_step = find(VMP(:,2)>=minStep & VMP(:,2)<=maxStep);
eggs = nest(:,1);                                               % eggs จะเป็นคำตอบ original ก่อนทำการ levy walk

u=randn(size(eggs))*sigma;
v=randn(size(eggs));
alpha = 10; 
step_length = round(u./abs(v).^(1/beta)*alpha);

feas_eggs = eggs + step_length;
chkStep = (min(feas_trans_step) <= feas_eggs) & (feas_eggs) <= max(feas_trans_step);
[r,c]= find(chkStep<1);

% ดูว่ามีคำตอบ infeasible solution หรือไม่
if size(r,1)>0, BadStepBit = 1; else BadStepBit = 0; end

while BadStepBit==1
    u=randn(size(r))*sigma;
    v=randn(size(r)); 
    step_length = round(u./abs(v).^(1/beta)*alpha);
    
    for row = 1:size(r,1)
        feas_eggs(r(row),1) = eggs(r(row),1) + step_length(row);
    end  
    
    chkStep = (min(feas_trans_step) <= feas_eggs) & (feas_eggs) <= max(feas_trans_step);
    [r,c]= find(chkStep<1);
    if size(r,1)>0, BadStepBit = 1; else BadStepBit = 0; end 
end  

%Check Fuel-Mix Constraint

CumCap(:,1) = VMP(feas_eggs,2)+ 5450;
CumOilCap(:,1) = VMP(feas_eggs,3).*ones(n,1)*200 + 550;                                              % เพราะ existing Oil มีอยู่ 550
CumLNGCap(:,1) = VMP(feas_eggs,4).*ones(n,1)*450 + 1400;                                             % เพราะ existing LNG มีอยู่ 1400
CumCoalCap(:,1) = VMP(feas_eggs,5).*ones(n,1)*500 + 1500;                                            % เพราะ existing Coal มีอยู่ 1500   
CumNuclearCap(:,1) = (VMP(feas_eggs,6).*ones(n,1)*1000 + 2000)+(VMP(feas_eggs,7).*ones(n,1)*700);     % เพราะ existing PWR มีอยู่ 2000

chkOil  = CumOilCap ./ CumCap <= Mmax(1);
chkLNG  = CumLNGCap ./ CumCap <= Mmax(2);
chkCoal = (Mmin(3) <= CumCoalCap ./ CumCap) & (CumCoalCap ./ CumCap) <= Mmax(3);
chkNuclear = (Mmin(4) <= CumNuclearCap ./ CumCap) & (CumNuclearCap ./ CumCap <= Mmax(4));

chkFuel = chkOil & chkLNG & chkCoal & chkNuclear;
[r,c]= find(chkFuel<1);                                                    % หาตำแหน่งที่ infeasible ในเงื่อนไข fuel diversification
if size(r,1)>0, FuelMixBit = 1; else FuelMixBit = 0; end                   % FuelMixBit = 1 กรณีมี infeasible solution ต้องแก้ bit ดังกล่าว

while FuelMixBit==1
    for row = 1:size(r,1)
        feas_eggs(r(row),1) = randi([min(feas_trans_step) max(feas_trans_step)],1,1);
    end 
    
    %Check Fuel-Mix Constraint
    CumCap(:,1) = VMP(feas_eggs,2)+ 5450;
    CumOilCap(:,1) = VMP(feas_eggs,3).*ones(n,1)*200 + 550;                                              % เพราะ existing Oil มีอยู่ 550
    CumLNGCap(:,1) = VMP(feas_eggs,4).*ones(n,1)*450 + 1400;                                             % เพราะ existing LNG มีอยู่ 1400
    CumCoalCap(:,1) = VMP(feas_eggs,5).*ones(n,1)*500 + 1500;                                            % เพราะ existing Coal มีอยู่ 1500   
    CumNuclearCap(:,1) = (VMP(feas_eggs,6).*ones(n,1)*1000 + 2000)+(VMP(feas_eggs,7).*ones(n,1)*700);     % เพราะ existing PWR มีอยู่ 2000

    chkOil  = CumOilCap ./ CumCap <= Mmax(1);
    chkLNG  = CumLNGCap ./ CumCap <= Mmax(2);
    chkCoal = (Mmin(3) <= CumCoalCap ./ CumCap) & (CumCoalCap ./ CumCap) <= Mmax(3);
    chkNuclear = (Mmin(4) <= CumNuclearCap ./ CumCap) & (CumNuclearCap ./ CumCap <= Mmax(4));

    chkFuel = chkOil & chkLNG & chkCoal & chkNuclear;
    [r,c]= find(chkFuel<1); 
    if size(r,1)>0, FuelMixBit = 1; else FuelMixBit = 0; end                   % FuelMixBit = 1 กรณีมี infeasible solution ต้องแก้ bit ดังกล่าว
end

for t=2:T  
   
eggs = nest(:,t);
for i=1:size(feas_eggs)
    minStep = LowerBound(t)-CumCap(i,t-1);
    maxStep = UpperBound(t)-CumCap(i,t-1);

    trans_feas_step = find(VMP(:,2)>=minStep & VMP(:,2)<=maxStep);
    %min(trans_feas_step)
    %max(trans_feas_step)
    u=randn()*sigma;
    v=randn();
    step_length = round(u./abs(v).^(1/beta)*alpha);
    feas_eggs(i,t) = eggs(i) + step_length;
    chkStep = (min(trans_feas_step) <= feas_eggs(i,t)) & (feas_eggs(i,t)) <= max(trans_feas_step);
    
    [r,c]= find(chkStep<1);
    if size(r,1)>0, BadStepBit = 1; else BadStepBit = 0; end
    
    while BadStepBit==1
    %u=randn()*sigma;
    %v=randn();
    %step_length = round(u./abs(v).^(1/beta)*alpha);
    %feas_eggs(i,t) = eggs(i) + step_length
    feas_eggs(i,t) = randi([min(trans_feas_step) max(trans_feas_step)],1,1);    
    chkStep = (min(trans_feas_step) <= feas_eggs(i,t)) & (feas_eggs(i,t)) <= max(trans_feas_step);
    [r,c]= find(chkStep<1);
    if size(r,1)>0, BadStepBit = 1; else BadStepBit = 0; end     
    end         
end  

%Check Fuel-Mix Constraint
CumCap(:,t) = VMP(feas_eggs(:,t),2)+ CumCap(:,t-1);
CumOilCap(:,t)  = VMP(feas_eggs(:,t),3).*ones(n,1)*200 + CumOilCap(:,t-1);   % column vector ของ cumulative oil cpacity
CumLNGCap(:,t)  = VMP(feas_eggs(:,t),4).*ones(n,1)*450 + CumLNGCap(:,t-1);   % column vector ของ cumulative LNG cpacity
CumCoalCap(:,t) = VMP(feas_eggs(:,t),5).*ones(n,1)*500 + CumCoalCap(:,t-1);   % column vector ของ cumulative Coal cpacity  
CumNuclearCap(:,t) = VMP(feas_eggs(:,t),6).*ones(n,1)*1000 + VMP(feas_eggs(:,t),7).*ones(n,1)*700 + CumNuclearCap(:,t-1);     % column vector ของ cumulative Nuclear cpacity

chkOil  = CumOilCap(:,t) ./ CumCap(:,t) <= Mmax(1);
chkLNG  = CumLNGCap(:,t) ./ CumCap(:,t) <= Mmax(2);
chkCoal = (Mmin(3) <= CumCoalCap(:,t) ./ CumCap(:,t)) & (CumCoalCap(:,t) ./ CumCap(:,t)) <= Mmax(3);
chkNuclear = (Mmin(4) <= CumNuclearCap(:,t) ./ CumCap(:,t)) & (CumNuclearCap(:,t) ./ CumCap(:,t) <= Mmax(4));

chkFuel = chkOil & chkLNG & chkCoal & chkNuclear;
[r,c]= find(chkFuel<1);                                                    % หาตำแหน่งที่ infeasible ในเงื่อนไข fuel diversification
if size(r,1)>0, FuelMixBit = 1; else FuelMixBit = 0; end                   % FuelMixBit = 1 กรณีมี infeasible solution ต้องแก้ bit ดังกล่าว

while FuelMixBit==1  
    for row = 1:size(r,1)
        minStep = LowerBound(t)-CumCap(r(row),t-1);
        maxStep = UpperBound(t)-CumCap(r(row),t-1);
        trans_feas_step = find(VMP(:,2)>=minStep & VMP(:,2)<=maxStep);
        %min(trans_feas_step)
        %max(trans_feas_step)
        feas_eggs(r(row),t) = randi([min(trans_feas_step) max(trans_feas_step)],1,1);
    end 

    %Check Fuel-Mix Constraint
    CumCap(:,t) = VMP(feas_eggs(:,t),2)+ CumCap(:,t-1);
    CumOilCap(:,t)  = VMP(feas_eggs(:,t),3).*ones(n,1)*200 + CumOilCap(:,t-1);   % column vector ของ cumulative oil cpacity
    CumLNGCap(:,t)  = VMP(feas_eggs(:,t),4).*ones(n,1)*450 + CumLNGCap(:,t-1);   % column vector ของ cumulative LNG cpacity
    CumCoalCap(:,t) = VMP(feas_eggs(:,t),5).*ones(n,1)*500 + CumCoalCap(:,t-1);   % column vector ของ cumulative Coal cpacity  
    CumNuclearCap(:,t) = VMP(feas_eggs(:,t),6).*ones(n,1)*1000 + VMP(feas_eggs(:,t),7).*ones(n,1)*700 + CumNuclearCap(:,t-1);     % column vector ของ cumulative Nuclear cpacity

    chkOil  = CumOilCap(:,t) ./ CumCap(:,t) <= Mmax(1);
    chkLNG  = CumLNGCap(:,t) ./ CumCap(:,t) <= Mmax(2);
    chkCoal = (Mmin(3) <= CumCoalCap(:,t) ./ CumCap(:,t)) & (CumCoalCap(:,t) ./ CumCap(:,t)) <= Mmax(3);
    chkNuclear = (Mmin(4) <= CumNuclearCap(:,t) ./ CumCap(:,t)) & (CumNuclearCap(:,t) ./ CumCap(:,t) <= Mmax(4));

    chkFuel = chkOil & chkLNG & chkCoal & chkNuclear;
    [r,c]= find(chkFuel<1);                                                    % หาตำแหน่งที่ infeasible ในเงื่อนไข fuel diversification
    if size(r,1)>0, FuelMixBit = 1; else FuelMixBit = 0; end                   % FuelMixBit = 1 กรณีมี infeasible solution ต้องแก้ bit ดังกล่าว
end

end

newnest = [bestnest; feas_eggs;];



