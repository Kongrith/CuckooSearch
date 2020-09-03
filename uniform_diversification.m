function new_nest= uniform_diversification(nest,Kbest,Pa)
bestnest= nest(Kbest,:);
nest(Kbest,:)=[];

global VMP LowerBound UpperBound Mmin Mmax T
% A fraction of worse nests are discovered with a probability pa

% Discovered or not -- a status vector
K=rand(size(nest,1),1)<= Pa;
discovered_nest = nest(K,:);
residue_nest = nest(~K,:);             % เป็นรังนกที่ไม่ได้ทำการ randomization
n = size(discovered_nest,1);

CumCap = zeros(n,1);
CumOilCap = zeros(n,1);
CumLNGCap = zeros(n,1);
CumCoalCap = zeros(n,1);
CumNuclearCap = zeros(n,1);
num_discover = size(find(K>0),1);
if num_discover > 0, DiscoveredBit = 1; else DiscoveredBit = 0; end
if  DiscoveredBit ==1
minStep = LowerBound(1)-5450;
maxStep = UpperBound(1)-5450;
       
feas_random = find(VMP(:,2)>=minStep & VMP(:,2)<=maxStep);  %เฉพาะ time stage ที่1 ที่จะ common กัน
feas_eggs = randi([min(feas_random) max(feas_random)],n,1);

%Check Fuel-Mix Constraint
CumCap(:,1) = VMP(feas_eggs,2) + 5450;
CumOilCap(:,1) = VMP(feas_eggs,3).*ones(n,1)*200 + 550;                                              % เพราะ existing Oil มีอยู่ 550
CumLNGCap(:,1) = VMP(feas_eggs,4).*ones(n,1)*450 + 1400;                                             % เพราะ existing LNG มีอยู่ 1400
CumCoalCap(:,1) = VMP(feas_eggs,5).*ones(n,1)*500 + 1500;                                            % เพราะ existing Coal มีอยู่ 1500   
CumNuclearCap(:,1) = (VMP(feas_eggs,6).*ones(n,1)*1000 + 2000)+(VMP(feas_eggs,7).*ones(n,1)*700);    % เพราะ existing PWR มีอยู่ 2000
    
chkOil  = CumOilCap ./ CumCap <= Mmax(1);
chkLNG  = CumLNGCap ./ CumCap <= Mmax(2);
chkCoal = (Mmin(3) <= CumCoalCap ./ CumCap) & (CumCoalCap ./ CumCap) <= Mmax(3);
chkNuclear = (Mmin(4) <= CumNuclearCap ./ CumCap) & (CumNuclearCap ./ CumCap <= Mmax(4));
    
chkFuel = chkOil & chkLNG & chkCoal & chkNuclear;
[r,c]= find(chkFuel<1);                                                    % หาตำแหน่งที่ infeasible ในเงื่อนไข fuel diversification
if size(r,1)>0, FuelMixBit = 1; else FuelMixBit = 0; end                   % FuelMixBit = 1 กรณีมี infeasible solution ต้องแก้ bit ดังกล่าว

while FuelMixBit==1
    for row = 1:size(r,1)
        feas_eggs(r(row),1) = randi([min(feas_random) max(feas_random)],1,1);
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
%eggs = nest(:,t);
for i=1:size(feas_eggs)
    minStep = LowerBound(t)-CumCap(i,t-1);
    maxStep = UpperBound(t)-CumCap(i,t-1);
    
    feas_random = find(VMP(:,2)>=minStep & VMP(:,2)<=maxStep);  %เฉพาะ time stage ที่1 ที่จะ common กัน
    feas_eggs(i,t) = randi([min(feas_random) max(feas_random)],1,1);      
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
    %feas_eggs
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

else
    feas_eggs =[];
end

new_nest = [bestnest; residue_nest; feas_eggs];


