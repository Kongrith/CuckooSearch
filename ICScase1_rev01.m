% ===============================================================  
% Doctoral seminar D3 ใช้ Cuckoo Search Optimization เพื่อวางแผน 6 ปี
% สำหรับจูนพารามิเตอร์ในระบบตัวอย่าง (เฉลย: วิธี DP = 5006.19)
% ===============================================================

%% Part0: problem initialization
clc, clear, close all

% rev01 นำต้นแบบไฟล์มาจาก case1_rev05_7

format short;

%help case1_rev1.m

% เป็นการสร้าง cuckoo search สำหรับ 3 time stages
% จำนวน combination ทั้งหมด = 7,077,888,000 
% global solution =  คำตอบ 4541.7 ล้าน dollarUS, เซตคำตอบ คือ U=[2 1 1; 4 2 0; 1 1 0;0 0 0; 1 1 1;] 
% disp('searching solution ... it may take a few minutes.');


tic
% Virtual Mapping Procedure  ใช้เวลา 0.127574 seconds ต่อ 1 ฟอง
global VMP
load VMP.mat;

% ข้อมูลทางเทคนิคของโรงไฟฟ้าที่นำเสนอเข้าแผน
global gendata CI Life VOM FOM 

% ข้อมูลทางเทคนิคของโรงไฟฟ้าที่นำเสนอเข้าแผน
% เรียกข้อมูลโรงไฟฟ้าที่ต่อในระบบ
%                | No |  Cap   |  FOR  |  VOM      |    FOM       | 
%                |    |  (MW)  |       |  ล้าน$/GWh | ล้าน$/MW-year |
existingdata =   [ 1     200     0.07     0.024         5.4 ;      % Oil#1 (Heavy oil)
                   2     200     0.068    0.027         5.4 ;      % Oil#2 (Heavy oil)
                   3     150     0.06     0.03          3.78 ;     % Oil#3 (Heavy oil)
                   4      50     0.03     0.043         2.712 ;    % LNG GT#1 (LNG)            
                   5      50     0.03     0.043         2.712 ;    % LNG GT#2 (LNG)           
                   6      50     0.03     0.043         2.712 ;    % LNG GT#3 (LNG)
                   7     400     0.10     0.038         7.824 ;    % LNG CC#1 (LNG) 
                   8     400     0.10     0.04          7.824 ;    % LNG CC#2 (LNG)            
                   9     450     0.11     0.035        10.8  ;     % LNG CC#3 (LNG) 
                  10     250     0.15     0.023        19.95 ;     % Coal#1 (Anthracite)
                  11     250     0.15     0.023        19.95 ;     % Coal#2 (Anthracite)
                  12     500     0.09     0.019        16.86 ;     % Coal#3 (Bituminous)
                  13     500     0.085    0.015        16.86 ;     % Coal#4 (Bituminous)
                  14    1000     0.09     0.005        59.28 ;     % Nuclear#1 (PWR)  
                  15    1000     0.088    0.005        55.2  ;];   % Nuclear#2 (PWR)
    
% เรียกข้อมูลโรงไฟฟ้าที่นำเสนอเข้าแผนในระบบ
%                | No | Cap  |  FOR  |  VOM      |    FOM       | Investment  |  Life  |
%                |    | (MW) |       |  ล้าน$/GWh | ล้าน$/MW-year |  ล้าน$/MW    |  year  |
candidatedata = [  16    200   0.070     0.021        5.28           162.5        25;
                   17    450   0.100     0.035        4.86           225.0        20; 
                   18    500   0.095     0.014        16.5           531.25       25;
                   19   1000   0.090     0.004        55.2          1625.0        25;
                   20    700   0.070     0.003        46.2          1225.0        25;];

gendata = vertcat(existingdata,candidatedata(:,1:5));            
                        
CI = candidatedata(:,6);                  % capital Investment (หน่วย: ล้าน$/MW)ที่นำเสนอเข้าแผนในระบบ
Life = candidatedata(:,7);                % Useful life (หน่วย: ปี )ที่นำเสนอเข้าแผนในระบบ

VOM = vertcat(existingdata(:,4),candidatedata(:,4));
FOM = vertcat(existingdata(:,5),candidatedata(:,5));

%% Part1: studied assumption
% หน่วยเงิน = ล้าน$, พลังงาน=GWh 
global T tt TT 
T=3;, k=5;, t0=2;                         % study time, ประเภทโรงไฟฟ้า candidate, เวลา reference ที่คิด discounting
inteval = 2;                              % ระยะปีระหว่าง time stage

tt= 0:inteval:T*inteval;                  % เวคเตอร์เวลา
TT = t0+inteval*T;                        % เวลาที่นับถอยหลังในการคิดมูลค่าซาก

% สมมติฐาน parameter ของปัญหาที่สนใจ
global penalty epsilon d CEENS Mmin Mmax
penalty = 10^6;
CEENS = 0.05;                                               % cost outage (หน่วย: ล้าน dollar per GWh)
epsilon = 0.01;                                             % LOLP criteria (หน่วย: 1วันต่อปี หรือ 0.27% หรือ 0.0027  )
d = 0.085;                                                  % discount_rate (หน่วย: 8.5% per year)
Mmin = [0 0 0.2 0.3]'; Mmax = [0.3 0.4 0.6 0.6]';           % หมายเหตุ: 1=oil-fired, 2=LNG-fired, 3=coal-fired, 4=nuclear
peak =[7000 9000 10000 12000 13000 14000 15000 17000 18000 20000 22000 24000];         % peak demand data (หน่วย: MW) ไม่นับปีอ้างอิง
Rmin = 0.15; Rmax = 0.6;                                    % min and max reserve margin (หน่วย: %)
existing = 5450;

% ส่วนของการคำนวณ EF0 ไว้ล่วงหน้า
global X_ref gcd N EF0  
X_ref(:,1) = ones(15,1);                                % existing_plant
gcd = min(gendata(:,2));                                % หา load block ที่เล็กที่สุด

loaddata = TrapezoidLoad(5000,gcd);                     % col1 MW, col2 freq, col3 indivifual prob., col4 cum. prob., col5 duration 
N = ceil(5000/gcd);                                     % หาจำนวน block สูงสุงของ load
EF0 = trapezoidhist2EF(loaddata,N,gcd);                 % คำนวน energy function ด้วยโหลดแบบ histogram สมมติฐานเป็นสี่เหลี่ยมคางหมู
     
GenUnit = merit_order(X_ref(:,1),gendata);
[E,LOLP(1),EENS(1)] = EEF(EF0,GenUnit,gcd,N,8760);      % หา E,LOLP,EENS ด้วยวิธี equivalent energy function(EEF)

EF0 = cell(1,T);                                        % กำหนดให้ EF0 เป็น cell array ไว้เก็บคำตอบคอลัมละ time stage
for t = 1:T
    loaddata = TrapezoidLoad(peak(t),gcd);
    N = ceil(peak(t)/gcd);
    EF0{1,t}= trapezoidhist2EF(loaddata,N,gcd);
end

% ส่วนของ parameter cuckoo search
n=2;, Pa=0.1;                                          % จำนวนรัง และ Discovery rate

ZVal=zeros(n,1);
fitness=ones(n,1)*10^10;                                % ถ้าเป็น minimization problems ให้ใส่ค่า fitness เยอะๆ 

% Initial Intelligent (เอา criteria power demand มาใช้ในการสร้างเลข)
global LowerBound UpperBound
LowerBound  =(1+Rmin)*peak(1:T);
UpperBound=(1+Rmax)*peak(1:T);

L1 = LowerBound(1)-existing;
U1 = UpperBound(1)-existing;

feas_set = find(VMP(:,2)>=L1 & VMP(:,2)<=U1);
feas_row = randi([min(feas_set) max(feas_set)],n,1)

%Check Fuel-Mix Constraint
global CumCap CumOilCap CumLNGCap CumCoalCap CumNuclearCap
CumCap = VMP(feas_row,2)+ existing;
CumOilCap = zeros(n,1);
CumLNGCap = zeros(n,1);
CumCoalCap = zeros(n,1);
CumNuclearCap = zeros(n,1);

CumOilCap = VMP(feas_row,3).*ones(n,1)*200 + 550;                                              % เพราะ existing Oil มีอยู่ 550
CumLNGCap = VMP(feas_row,4).*ones(n,1)*450 + 1400;                                             % เพราะ existing LNG มีอยู่ 1400
CumCoalCap = VMP(feas_row,5).*ones(n,1)*500 + 1500;                                            % เพราะ existing Coal มีอยู่ 1500   
CumNuclearCap = (VMP(feas_row,6).*ones(n,1)*1000 + 2000)+(VMP(feas_row,7).*ones(n,1)*700);     % เพราะ existing PWR มีอยู่ 2000

chkOil  = CumOilCap ./ CumCap <= Mmax(1);
chkLNG  = CumLNGCap ./ CumCap <= Mmax(2);
chkCoal = (Mmin(3) <= CumCoalCap ./ CumCap) & (CumCoalCap ./ CumCap) <= Mmax(3);
chkNuclear = (Mmin(4) <= CumNuclearCap ./ CumCap) & (CumNuclearCap ./ CumCap <= Mmax(4));

chkFuel = chkOil & chkLNG & chkCoal & chkNuclear;
[r,c]= find(chkFuel<1);                                                    % หาตำแหน่งที่ infeasible ในเงื่อนไข fuel diversification

% ดูว่ามีคำตอบ infeasible solution หรือไม่
if size(r,1)>0, FuelMixBit = 1; else FuelMixBit = 0; end                   % FuelMixBit = 1 กรณีมี infeasible solution ต้องแก้ bit ดังกล่าว
while FuelMixBit==1
    feas_row(r) = randi([min(feas_set) max(feas_set)],size(r,1),1);

    CumCap = VMP(feas_row,2)+ existing;
    CumOilCap = VMP(feas_row,3).*ones(n,1)*200 + 550;                                              % เพราะ existing Oil มีอยู่ 550
    CumLNGCap = VMP(feas_row,4).*ones(n,1)*450 + 1400;                                             % เพราะ existing LNG มีอยู่ 1400
    CumCoalCap = VMP(feas_row,5).*ones(n,1)*500 + 1500;                                            % เพราะ existing Coal มีอยู่ 1500   
    CumNuclearCap = (VMP(feas_row,6).*ones(n,1)*1000 + 2000)+(VMP(feas_row,7).*ones(n,1)*700);     % เพราะ existing PWR มีอยู่ 2000

    chkOil  = CumOilCap ./ CumCap <= Mmax(1);
    chkLNG  = CumLNGCap ./ CumCap <= Mmax(2);
    chkCoal = (Mmin(3) <= CumCoalCap ./ CumCap) & (CumCoalCap ./ CumCap) <= Mmax(3);
    chkNuclear = (Mmin(4) <= CumNuclearCap ./ CumCap) & (CumNuclearCap ./ CumCap <= Mmax(4));

    chkFuel = chkOil & chkLNG & chkCoal & chkNuclear;   
    [r,c]= find(chkFuel<1);
    if size(r,1)>0, FuelMixBit = 1; else FuelMixBit = 0; end 
end

for t=2:T                                 % 0.019289 seconds ของ 2 รัง และ 0.022629 seconds.สำหรับ 25 รัง
%IIICumCap = existing + VMP(feas_row,2);
minRequire = LowerBound(t)-CumCap(:,t-1);
maxRequire = UpperBound(t)-CumCap(:,t-1);

trans_feas_row = zeros(n,1);
    for i=1:n
        trans_feas_set = find(VMP(:,2)>=minRequire(i) & VMP(:,2)<=maxRequire(i));
        trans_feas_row(i,1) = randi([min(trans_feas_set) max(trans_feas_set)],1,1);
    end

%Check Fuel-Mix Constraint
CumCap(:,t) = VMP(trans_feas_row,2)+ CumCap(:,t-1);
CumOilCap(:,t)  = VMP(trans_feas_row,3).*ones(n,1)*200 + CumOilCap(:,t-1);   % column vector ของ cumulative oil cpacity
CumLNGCap(:,t)  = VMP(trans_feas_row,4).*ones(n,1)*450 + CumLNGCap(:,t-1);   % column vector ของ cumulative LNG cpacity
CumCoalCap(:,t) = VMP(trans_feas_row,5).*ones(n,1)*500 + CumCoalCap(:,t-1);   % column vector ของ cumulative Coal cpacity  
CumNuclearCap(:,t) = VMP(trans_feas_row,6).*ones(n,1)*1000 + VMP(trans_feas_row,7).*ones(n,1)*700 + CumNuclearCap(:,t-1);     % column vector ของ cumulative Nuclear cpacity

chkOil  = CumOilCap ./ CumCap <= Mmax(1);
chkLNG  = CumLNGCap ./ CumCap <= Mmax(2);
chkCoal = (Mmin(3) <= CumCoalCap ./ CumCap) & (CumCoalCap ./ CumCap) <= Mmax(3);
chkNuclear = (Mmin(4) <= CumNuclearCap ./ CumCap) & (CumNuclearCap ./ CumCap <= Mmax(4));

chkFuel = chkOil & chkLNG & chkCoal & chkNuclear;
[r,c]= find(chkFuel<1);     
        
% ลอกคำตอบมาเก็บไว้
feas_row(:,t) = trans_feas_row;

% ดูว่ามีคำตอบ infeasible solution หรือไม่
if size(r,1)>0, FuelMixBit = 1; else FuelMixBit = 0; end                   % FuelMixBit = 1 กรณีมี infeasible solution ต้องแก้ bit ดังกล่าว
while FuelMixBit==1
    trans_feas_row = zeros(size(r,1),1);
    for row = 1:size(r,1)
        trans_feas_set = find(VMP(:,2)>=minRequire(r(row)) & VMP(:,2)<=maxRequire(r(row)));
        %trans_feas_row(row,1) = randi([min(trans_feas_set) max(trans_feas_set)],1,1)
        feas_row(r(row),t) = randi([min(trans_feas_set) max(trans_feas_set)],1,1);
    end
    %feas_row(:,t)= trans_feas_row(:,1);
    %trans_feas_set = find(VMP(:,2)>=minRequire(r) & VMP(:,2)<=maxRequire(r))
    %trans_feas_row(r) = randi([min(trans_feas_set) max(trans_feas_set)],size(r,1),1);
    CumCap(:,t) = VMP(feas_row(:,t),2)+ CumCap(:,t-1);
    CumOilCap(:,t)  = VMP(feas_row(:,t),3).*ones(n,1)*200 + CumOilCap(:,t-1);   % column vector ของ cumulative oil cpacity
    CumLNGCap(:,t)  = VMP(feas_row(:,t),4).*ones(n,1)*450 + CumLNGCap(:,t-1);   % column vector ของ cumulative LNG cpacity
    CumCoalCap(:,t) = VMP(feas_row(:,t),5).*ones(n,1)*500 + CumCoalCap(:,t-1);   % column vector ของ cumulative Coal cpacity  
    CumNuclearCap(:,t) = VMP(feas_row(:,t),6).*ones(n,1)*1000 + VMP(feas_row(:,t),7).*ones(n,1)*700 + CumNuclearCap(:,t-1);     % column vector ของ cumulative Nuclear cpacity

    chkOil  = CumOilCap ./ CumCap <= Mmax(1);
    chkLNG  = CumLNGCap ./ CumCap <= Mmax(2);
    chkCoal = (Mmin(3) <= CumCoalCap ./ CumCap) & (CumCoalCap ./ CumCap) <= Mmax(3);
    chkNuclear = (Mmin(4) <= CumNuclearCap ./ CumCap) & (CumNuclearCap ./ CumCap <= Mmax(4));
    chkFuel = chkOil & chkLNG & chkCoal & chkNuclear;   
    [r,c]= find(chkFuel<1);   
    if size(r,1)>0, FuelMixBit = 1; else FuelMixBit = 0; end     
end

end
nest = feas_row;
[fitness,ZVal]= CalObj_rev01(feas_row,fitness,ZVal);

[fmin,Kbest]=min(fitness);                               % หาค่าตำตอบที่มี fmin จาก vector คำตอบ
bestnest= nest(Kbest,:);                                 % หาเวคเตอร์ตำตอบที่มี fmin น้อยที่สุด
iter = 0;  seed=0;
tagetBit = 0;
%Zvector = zeros(maxiter,1);
maxiter= 1;
%% starting optimization
while (iter<=maxiter) & (tagetBit == 0)     
    % intensification/exploitation via levy flight
      new_nest=intensification(nest,Kbest);   
      [nest,fitness,ZVal]=CalcFitness(nest,new_nest,fitness,ZVal);   
      seed=seed+n;                                      % Update search seeds again                                
    % diversification/exploration via randomization
      new_nest=uniform_diversification(nest,Kbest,Pa);
      [nest,fitness,ZVal]=CalcFitness(nest,new_nest,fitness,ZVal);    
      %seed = seed+n;                                    % Update search seeds again
    % update bestnest    
      [fmin,Kbest]=min(fitness);                        % หาค่าตำตอบที่มี fmin จาก vector คำตอบ
      bestnest=nest(Kbest,:);                           % หาตำตอบที่มี fmin น้อยที่สุด
      Zmin = ZVal(Kbest);
      %[Zmin,I] = min(ZVal);                             % หาตำตอบที่มี Zmin น้อยที่สุด
      %if fmin < 6646 & (tagetBit == 0)
          %iter
          %tagetBit = 1;
      %end    
      iter = iter+1;                                    % Update the counter
      %Zvector(iter)=Zmin;
end %% End of iterations
toc 
Zmin

if Zmin == 4740.042073
    disp('global solution found')
end
%fmin
%bestnest
%seed
%Zvector
