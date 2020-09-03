function [NewFitness,NewZ] = CalObj_rev01(feas_row,fitness,ZVal)
% 
% เวลาในการทำงานของ function 0.692367 วินาที ต่อ 25 รัง 7 time stage
% เพิ่ม outage cost ในสมการ Z

global T tt TT VOM FOM  CI Life d  EF0 X_ref gendata N gcd VMP penalty epsilon CEENS

% เริ่มต้น
I = zeros(1,T); M = zeros(1,T); LOLP = zeros(1,T); EENS = zeros(1,T);
 
%NewZ = zeros(size(nest,1),1); 
%NewFitness = zeros(size(nest,1),1);

NewZ = zeros(size(feas_row,1),1); 
NewFitness = zeros(size(feas_row,1),1);

for j=1:size(feas_row,1)                %% loop นอกสุด    %0.005
%U = RecallAns(nest,j,T,k)
feas_row
Strings1 = VMP(feas_row(j,1:T),3);
Strings2 = VMP(feas_row(j,1:T),4);
Strings3 = VMP(feas_row(j,1:T),5);
Strings4 = VMP(feas_row(j,1:T),6);
Strings5 = VMP(feas_row(j,1:T),7);

U = [Strings1 Strings2 Strings3 Strings4 Strings5]'; % จะได้คำตอบเป็น U = [nest1; nest2; .. nestสุดท้าย;]

for t = 1:T                          %% loop ใน 
    if t == 1                        % หา state variables
        X(1:15,t)= X_ref;
        X(16:20,t)= U(:,t);
    else
        X(1:15,t) = X(1:15,t-1);                  
        X(16:20,t) = X(16:20,t-1) + U(:,t);
    end
    
    GenUnit = merit_order(X(:,t),gendata);   % ทำการเรียงลำดับโรงไฟฟ้าด้วย merit order %col1 คือ capacity, col2 คือ FOR, col3 คือ operating, col4 คือ fix O&M
    [E,LOLP(1,t),EENS(1,t)] = EEF(EF0{1,t},GenUnit,gcd,N,8760);
           
    I(1,t) = sum((CI.*U(:,t))/(1+d)^tt(t+1));
    M(1,t) = ((sum(FOM.*X(:,t)) + sum(VOM(GenUnit(:,1)).*(E/1000)))/(1+d)^(tt(t+1)+0.5)) + ...
             ((sum(FOM.*X(:,t)) + sum(VOM(GenUnit(:,1)).*(E/1000)))/(1+d)^(tt(t+1)+0.5+1));
    O(1,t) = (((EENS(1,t)*CEENS)/(1+d)^(tt(t+1)+0.5)) + ((EENS(1,t)*CEENS)/(1+d)^(tt(t+1)+0.5+1)))/1000;
    TotalGen(t)=sum(GenUnit(:,2));
        
end

% คำนวนค่า spent life ของโรงไฟฟ้าสำหรับการคำนวนค่าซากของโรงไฟฟ้าจากเวคเตอร U
row = 1; col = 1; maxrow = size(U(:,1:T),1); maxcol = size(U(:,1:T),2); IDLife = zeros(sum(sum(U)),1);  % sum(U) คือ บวกลงมาในแนวคอลัม
for RowIndex = 1:maxrow
    for ColIndex = 1:maxcol
        numrepeat = U(RowIndex,ColIndex);        
        for r = 1:numrepeat
            count = 1;           
            for fillcol = col+1:maxcol+1                       % col+1 เพราะเราต้องนับอายุ 1 time stage = 2 ปี  ส่วน mxacol+1 เพราะ salvage ต้องนับที่ end of studied period
                SpentLife(row,fillcol)= 2 + 2*(count-1);       % อนุกรม 2 4 6 8 10
                count = count+1;
            end
            IDLife(row,1) = RowIndex;
            row = row + 1;
        end
        col = col+1;
    end
    col = col - maxcol;
end

P = CI(IDLife); L = 0.157669630022891*CI(IDLife);
S = sum(P-(P-L)./Life(IDLife).*SpentLife(:,end))/(1+d)^TT;

%NewZ(j) = sum(I)+ sum(M)- S;                                    % evaluate objective function
NewZ(j) = sum(I)+ sum(M) + sum(O) - S;                                    % evaluate objective function

% สมการเงื่อนไข LOLP
chkLOLP = LOLP>epsilon ;                                        % ตรวจสอบ LOLP criteria: ได้ 0 คือ ไม่เกิน criteria
TotalPenalty = sum(chkLOLP.* O * penalty);

NewFitness(j) = NewZ(j) + TotalPenalty ;                          %PenaltyFunction(LOLP);
SpentLife =[];, RemainLife =[];, SFF=[];, SCAF=[];              % ล้างค่า เพื่อไม่ให้เกิด error
end

%function [U] = RecallAns(nest,num,T,k)
% ดึงคำตอบจาก nest มาแสดงแบบ [เซตคำตอบของโรงไฟ้า1; เซตคำตอบของโรงไฟ้า2, ..., เซตคำตอบของโรงไฟ้าk]
%ans=[];, setans=[];
%for i=1:k
    %numseq = 1+(i-1)*T;              % กระบวนการ map คำตอบจะเป็นแบบอนุกรมเลขคณิต 
    %ans = nest(num,numseq:numseq+T-1);
    %setans = vertcat(setans,ans);
    %ans =[];
%end
%U = setans;
%function [TotalPenalty] = PenaltyFunction(LOLP)
%global penalty epsilon
%TotalPenalty=0;

% สมการเงื่อนไข LOLP
%chkLOLP = LOLP>epsilon                                        % ตรวจสอบ LOLP criteria: ได้ 0 คือ ไม่เกิน criteria
%TotalPenalty = TotalPenalty + sum(chkLOLP*penalty);





