% ===============================================================  
% Doctoral seminar D3 �� Cuckoo Search Optimization �����ҧἹ 6 ��
% ����Ѻ�ٹ������������к�������ҧ (���: �Ը� DP = 5006.19)
% ===============================================================

%% Part0: problem initialization
clc, clear, close all

% rev01 �ӵ�Ẻ����Ҩҡ case1_rev05_7

format short;

%help case1_rev1.m

% �繡�����ҧ cuckoo search ����Ѻ 3 time stages
% �ӹǹ combination ������ = 7,077,888,000 
% global solution =  �ӵͺ 4541.7 ��ҹ dollarUS, ૵�ӵͺ ��� U=[2 1 1; 4 2 0; 1 1 0;0 0 0; 1 1 1;] 
% disp('searching solution ... it may take a few minutes.');


tic
% Virtual Mapping Procedure  ������ 0.127574 seconds ��� 1 �ͧ
global VMP
load VMP.mat;

% �����ŷҧ෤�Ԥ�ͧ�ç俿�ҷ����ʹ����Ἱ
global gendata CI Life VOM FOM 

% �����ŷҧ෤�Ԥ�ͧ�ç俿�ҷ����ʹ����Ἱ
% ���¡�������ç俿�ҷ������к�
%                | No |  Cap   |  FOR  |  VOM      |    FOM       | 
%                |    |  (MW)  |       |  ��ҹ$/GWh | ��ҹ$/MW-year |
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
    
% ���¡�������ç俿�ҷ����ʹ����Ἱ��к�
%                | No | Cap  |  FOR  |  VOM      |    FOM       | Investment  |  Life  |
%                |    | (MW) |       |  ��ҹ$/GWh | ��ҹ$/MW-year |  ��ҹ$/MW    |  year  |
candidatedata = [  16    200   0.070     0.021        5.28           162.5        25;
                   17    450   0.100     0.035        4.86           225.0        20; 
                   18    500   0.095     0.014        16.5           531.25       25;
                   19   1000   0.090     0.004        55.2          1625.0        25;
                   20    700   0.070     0.003        46.2          1225.0        25;];

gendata = vertcat(existingdata,candidatedata(:,1:5));            
                        
CI = candidatedata(:,6);                  % capital Investment (˹���: ��ҹ$/MW)�����ʹ����Ἱ��к�
Life = candidatedata(:,7);                % Useful life (˹���: �� )�����ʹ����Ἱ��к�

VOM = vertcat(existingdata(:,4),candidatedata(:,4));
FOM = vertcat(existingdata(:,5),candidatedata(:,5));

%% Part1: studied assumption
% ˹����Թ = ��ҹ$, ��ѧ�ҹ=GWh 
global T tt TT 
T=3;, k=5;, t0=2;                         % study time, �������ç俿�� candidate, ���� reference ���Դ discounting
inteval = 2;                              % ���л������ҧ time stage

tt= 0:inteval:T*inteval;                  % �Ǥ��������
TT = t0+inteval*T;                        % ���ҷ��Ѻ�����ѧ㹡�äԴ��Ť�ҫҡ

% ����԰ҹ parameter �ͧ�ѭ�ҷ��ʹ�
global penalty epsilon d CEENS Mmin Mmax
penalty = 10^6;
CEENS = 0.05;                                               % cost outage (˹���: ��ҹ dollar per GWh)
epsilon = 0.01;                                             % LOLP criteria (˹���: 1�ѹ��ͻ� ���� 0.27% ���� 0.0027  )
d = 0.085;                                                  % discount_rate (˹���: 8.5% per year)
Mmin = [0 0 0.2 0.3]'; Mmax = [0.3 0.4 0.6 0.6]';           % �����˵�: 1=oil-fired, 2=LNG-fired, 3=coal-fired, 4=nuclear
peak =[7000 9000 10000 12000 13000 14000 15000 17000 18000 20000 22000 24000];         % peak demand data (˹���: MW) ���Ѻ����ҧ�ԧ
Rmin = 0.15; Rmax = 0.6;                                    % min and max reserve margin (˹���: %)
existing = 5450;

% ��ǹ�ͧ��äӹǳ EF0 �����ǧ˹��
global X_ref gcd N EF0  
X_ref(:,1) = ones(15,1);                                % existing_plant
gcd = min(gendata(:,2));                                % �� load block �����硷���ش

loaddata = TrapezoidLoad(5000,gcd);                     % col1 MW, col2 freq, col3 indivifual prob., col4 cum. prob., col5 duration 
N = ceil(5000/gcd);                                     % �Ҩӹǹ block �٧�ا�ͧ load
EF0 = trapezoidhist2EF(loaddata,N,gcd);                 % �ӹǹ energy function ������ŴẺ histogram ����԰ҹ�������������ҧ���
     
GenUnit = merit_order(X_ref(:,1),gendata);
[E,LOLP(1),EENS(1)] = EEF(EF0,GenUnit,gcd,N,8760);      % �� E,LOLP,EENS �����Ը� equivalent energy function(EEF)

EF0 = cell(1,T);                                        % ��˹���� EF0 �� cell array ����纤ӵͺ������� time stage
for t = 1:T
    loaddata = TrapezoidLoad(peak(t),gcd);
    N = ceil(peak(t)/gcd);
    EF0{1,t}= trapezoidhist2EF(loaddata,N,gcd);
end

% ��ǹ�ͧ parameter cuckoo search
n=2;, Pa=0.1;                                          % �ӹǹ�ѧ ��� Discovery rate

ZVal=zeros(n,1);
fitness=ones(n,1)*10^10;                                % ����� minimization problems �������� fitness ����� 

% Initial Intelligent (��� criteria power demand ����㹡�����ҧ�Ţ)
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

CumOilCap = VMP(feas_row,3).*ones(n,1)*200 + 550;                                              % ���� existing Oil ������ 550
CumLNGCap = VMP(feas_row,4).*ones(n,1)*450 + 1400;                                             % ���� existing LNG ������ 1400
CumCoalCap = VMP(feas_row,5).*ones(n,1)*500 + 1500;                                            % ���� existing Coal ������ 1500   
CumNuclearCap = (VMP(feas_row,6).*ones(n,1)*1000 + 2000)+(VMP(feas_row,7).*ones(n,1)*700);     % ���� existing PWR ������ 2000

chkOil  = CumOilCap ./ CumCap <= Mmax(1);
chkLNG  = CumLNGCap ./ CumCap <= Mmax(2);
chkCoal = (Mmin(3) <= CumCoalCap ./ CumCap) & (CumCoalCap ./ CumCap) <= Mmax(3);
chkNuclear = (Mmin(4) <= CumNuclearCap ./ CumCap) & (CumNuclearCap ./ CumCap <= Mmax(4));

chkFuel = chkOil & chkLNG & chkCoal & chkNuclear;
[r,c]= find(chkFuel<1);                                                    % �ҵ��˹觷�� infeasible ����͹� fuel diversification

% ������դӵͺ infeasible solution �������
if size(r,1)>0, FuelMixBit = 1; else FuelMixBit = 0; end                   % FuelMixBit = 1 �ó��� infeasible solution ��ͧ�� bit �ѧ�����
while FuelMixBit==1
    feas_row(r) = randi([min(feas_set) max(feas_set)],size(r,1),1);

    CumCap = VMP(feas_row,2)+ existing;
    CumOilCap = VMP(feas_row,3).*ones(n,1)*200 + 550;                                              % ���� existing Oil ������ 550
    CumLNGCap = VMP(feas_row,4).*ones(n,1)*450 + 1400;                                             % ���� existing LNG ������ 1400
    CumCoalCap = VMP(feas_row,5).*ones(n,1)*500 + 1500;                                            % ���� existing Coal ������ 1500   
    CumNuclearCap = (VMP(feas_row,6).*ones(n,1)*1000 + 2000)+(VMP(feas_row,7).*ones(n,1)*700);     % ���� existing PWR ������ 2000

    chkOil  = CumOilCap ./ CumCap <= Mmax(1);
    chkLNG  = CumLNGCap ./ CumCap <= Mmax(2);
    chkCoal = (Mmin(3) <= CumCoalCap ./ CumCap) & (CumCoalCap ./ CumCap) <= Mmax(3);
    chkNuclear = (Mmin(4) <= CumNuclearCap ./ CumCap) & (CumNuclearCap ./ CumCap <= Mmax(4));

    chkFuel = chkOil & chkLNG & chkCoal & chkNuclear;   
    [r,c]= find(chkFuel<1);
    if size(r,1)>0, FuelMixBit = 1; else FuelMixBit = 0; end 
end

for t=2:T                                 % 0.019289 seconds �ͧ 2 �ѧ ��� 0.022629 seconds.����Ѻ 25 �ѧ
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
CumOilCap(:,t)  = VMP(trans_feas_row,3).*ones(n,1)*200 + CumOilCap(:,t-1);   % column vector �ͧ cumulative oil cpacity
CumLNGCap(:,t)  = VMP(trans_feas_row,4).*ones(n,1)*450 + CumLNGCap(:,t-1);   % column vector �ͧ cumulative LNG cpacity
CumCoalCap(:,t) = VMP(trans_feas_row,5).*ones(n,1)*500 + CumCoalCap(:,t-1);   % column vector �ͧ cumulative Coal cpacity  
CumNuclearCap(:,t) = VMP(trans_feas_row,6).*ones(n,1)*1000 + VMP(trans_feas_row,7).*ones(n,1)*700 + CumNuclearCap(:,t-1);     % column vector �ͧ cumulative Nuclear cpacity

chkOil  = CumOilCap ./ CumCap <= Mmax(1);
chkLNG  = CumLNGCap ./ CumCap <= Mmax(2);
chkCoal = (Mmin(3) <= CumCoalCap ./ CumCap) & (CumCoalCap ./ CumCap) <= Mmax(3);
chkNuclear = (Mmin(4) <= CumNuclearCap ./ CumCap) & (CumNuclearCap ./ CumCap <= Mmax(4));

chkFuel = chkOil & chkLNG & chkCoal & chkNuclear;
[r,c]= find(chkFuel<1);     
        
% �͡�ӵͺ�������
feas_row(:,t) = trans_feas_row;

% ������դӵͺ infeasible solution �������
if size(r,1)>0, FuelMixBit = 1; else FuelMixBit = 0; end                   % FuelMixBit = 1 �ó��� infeasible solution ��ͧ�� bit �ѧ�����
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
    CumOilCap(:,t)  = VMP(feas_row(:,t),3).*ones(n,1)*200 + CumOilCap(:,t-1);   % column vector �ͧ cumulative oil cpacity
    CumLNGCap(:,t)  = VMP(feas_row(:,t),4).*ones(n,1)*450 + CumLNGCap(:,t-1);   % column vector �ͧ cumulative LNG cpacity
    CumCoalCap(:,t) = VMP(feas_row(:,t),5).*ones(n,1)*500 + CumCoalCap(:,t-1);   % column vector �ͧ cumulative Coal cpacity  
    CumNuclearCap(:,t) = VMP(feas_row(:,t),6).*ones(n,1)*1000 + VMP(feas_row(:,t),7).*ones(n,1)*700 + CumNuclearCap(:,t-1);     % column vector �ͧ cumulative Nuclear cpacity

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

[fmin,Kbest]=min(fitness);                               % �Ҥ�ҵӵͺ����� fmin �ҡ vector �ӵͺ
bestnest= nest(Kbest,:);                                 % ���Ǥ����ӵͺ����� fmin ���·���ش
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
      [fmin,Kbest]=min(fitness);                        % �Ҥ�ҵӵͺ����� fmin �ҡ vector �ӵͺ
      bestnest=nest(Kbest,:);                           % �ҵӵͺ����� fmin ���·���ش
      Zmin = ZVal(Kbest);
      %[Zmin,I] = min(ZVal);                             % �ҵӵͺ����� Zmin ���·���ش
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
