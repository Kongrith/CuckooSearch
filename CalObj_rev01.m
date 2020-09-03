function [NewFitness,NewZ] = CalObj_rev01(feas_row,fitness,ZVal)
% 
% ����㹡�÷ӧҹ�ͧ function 0.692367 �Թҷ� ��� 25 �ѧ 7 time stage
% ���� outage cost ������ Z

global T tt TT VOM FOM  CI Life d  EF0 X_ref gendata N gcd VMP penalty epsilon CEENS

% �������
I = zeros(1,T); M = zeros(1,T); LOLP = zeros(1,T); EENS = zeros(1,T);
 
%NewZ = zeros(size(nest,1),1); 
%NewFitness = zeros(size(nest,1),1);

NewZ = zeros(size(feas_row,1),1); 
NewFitness = zeros(size(feas_row,1),1);

for j=1:size(feas_row,1)                %% loop �͡�ش    %0.005
%U = RecallAns(nest,j,T,k)
feas_row
Strings1 = VMP(feas_row(j,1:T),3);
Strings2 = VMP(feas_row(j,1:T),4);
Strings3 = VMP(feas_row(j,1:T),5);
Strings4 = VMP(feas_row(j,1:T),6);
Strings5 = VMP(feas_row(j,1:T),7);

U = [Strings1 Strings2 Strings3 Strings4 Strings5]'; % ����ӵͺ�� U = [nest1; nest2; .. nest�ش����;]

for t = 1:T                          %% loop � 
    if t == 1                        % �� state variables
        X(1:15,t)= X_ref;
        X(16:20,t)= U(:,t);
    else
        X(1:15,t) = X(1:15,t-1);                  
        X(16:20,t) = X(16:20,t-1) + U(:,t);
    end
    
    GenUnit = merit_order(X(:,t),gendata);   % �ӡ�����§�ӴѺ�ç俿�Ҵ��� merit order %col1 ��� capacity, col2 ��� FOR, col3 ��� operating, col4 ��� fix O&M
    [E,LOLP(1,t),EENS(1,t)] = EEF(EF0{1,t},GenUnit,gcd,N,8760);
           
    I(1,t) = sum((CI.*U(:,t))/(1+d)^tt(t+1));
    M(1,t) = ((sum(FOM.*X(:,t)) + sum(VOM(GenUnit(:,1)).*(E/1000)))/(1+d)^(tt(t+1)+0.5)) + ...
             ((sum(FOM.*X(:,t)) + sum(VOM(GenUnit(:,1)).*(E/1000)))/(1+d)^(tt(t+1)+0.5+1));
    O(1,t) = (((EENS(1,t)*CEENS)/(1+d)^(tt(t+1)+0.5)) + ((EENS(1,t)*CEENS)/(1+d)^(tt(t+1)+0.5+1)))/1000;
    TotalGen(t)=sum(GenUnit(:,2));
        
end

% �ӹǹ��� spent life �ͧ�ç俿������Ѻ��äӹǹ��ҫҡ�ͧ�ç俿�Ҩҡ�Ǥ��� U
row = 1; col = 1; maxrow = size(U(:,1:T),1); maxcol = size(U(:,1:T),2); IDLife = zeros(sum(sum(U)),1);  % sum(U) ��� �ǡŧ����Ǥ����
for RowIndex = 1:maxrow
    for ColIndex = 1:maxcol
        numrepeat = U(RowIndex,ColIndex);        
        for r = 1:numrepeat
            count = 1;           
            for fillcol = col+1:maxcol+1                       % col+1 ������ҵ�ͧ�Ѻ���� 1 time stage = 2 ��  ��ǹ mxacol+1 ���� salvage ��ͧ�Ѻ��� end of studied period
                SpentLife(row,fillcol)= 2 + 2*(count-1);       % ͹ء�� 2 4 6 8 10
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

% ��������͹� LOLP
chkLOLP = LOLP>epsilon ;                                        % ��Ǩ�ͺ LOLP criteria: �� 0 ��� ����Թ criteria
TotalPenalty = sum(chkLOLP.* O * penalty);

NewFitness(j) = NewZ(j) + TotalPenalty ;                          %PenaltyFunction(LOLP);
SpentLife =[];, RemainLife =[];, SFF=[];, SCAF=[];              % ��ҧ��� �����������Դ error
end

%function [U] = RecallAns(nest,num,T,k)
% �֧�ӵͺ�ҡ nest ���ʴ�Ẻ [૵�ӵͺ�ͧ�ç���1; ૵�ӵͺ�ͧ�ç���2, ..., ૵�ӵͺ�ͧ�ç���k]
%ans=[];, setans=[];
%for i=1:k
    %numseq = 1+(i-1)*T;              % ��кǹ��� map �ӵͺ����Ẻ͹ء���Ţ��Ե 
    %ans = nest(num,numseq:numseq+T-1);
    %setans = vertcat(setans,ans);
    %ans =[];
%end
%U = setans;
%function [TotalPenalty] = PenaltyFunction(LOLP)
%global penalty epsilon
%TotalPenalty=0;

% ��������͹� LOLP
%chkLOLP = LOLP>epsilon                                        % ��Ǩ�ͺ LOLP criteria: �� 0 ��� ����Թ criteria
%TotalPenalty = TotalPenalty + sum(chkLOLP*penalty);





