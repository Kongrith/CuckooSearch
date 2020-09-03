function [Eg,LOLP,EENS] = EEF(EF0,GenUnit,gcd,N,time)
% �ӹǹ�� LOLP ��� EENS �����Ը� equivalent energy function method (EEF)
% ��͸Ժ��: N ��� block �٧�ا�ͧ load , time ��� period of study , EF0 ��� energy function

numGen = size(GenUnit(:,1),1);       % �ӹǹ���駷���ͧ�� convolution
Jn = sum(GenUnit(:,2))/gcd;          % �Ҩӹǹ block ����ͧ shift �ҡ gen unit     

Eg=zeros(numGen,1);                  % ���ҧ vector ����Ѻ�纤�� energy ��� gen ��Ե��
Ed=zeros(numGen,1);                  % ���ҧ vector ����Ѻ�纤�� expected energy not served
EF = zeros(Jn+N,numGen);             % ���ҧ ����ԡ������Ѻ energy function
EF_tmp = zeros(Jn+N,numGen);         % ���ҧ ����ԡ������Ѻ energy function ��� shift � K

EF0 = vertcat(EF0,zeros(Jn,1));      % �������������� ����ԡ� EF0 (energy function)
totalE = sum(EF0);                   % ��鹷�����ҿ�ͧ Load ������
EENS = 0; index = 0;

for i = 1:numGen
   q = GenUnit(i,3);                 % ��� FOR �ͧ Generator unit
   p = 1-q;                          % ��� availability �ͧ Generator unit
   K = GenUnit(i,2)/gcd;             % �ӹǹ block ������ʹ�
   index = index + K;                %
   
   if i==1                           % �ӡ�� shift energy function 价��� K ���ͤ
       EF_tmp(K,i)=EF0(1);
       for j=1:size(EF0,1)
           EF_tmp(j+K,i)=EF0(j);
       end        
   else
       for j=1:size(EF,1)
           EF_tmp(j+K,i)=EF(j,i-1);    
       end
   end  
   
   if i==1                           % convolution
       for j=index:size(EF0,1)
           EF(j,i)= p*EF0(j,i) + q*EF_tmp(j,i);
       end       
   else
       for j=index:size(EF,1)
           EF(j,i)= p*EF(j,i-1) + q*EF_tmp(j,i);    
       end
   end   
      
end

index = 0;
for i=1:numGen                       % �ӹǹ�Ҥ�� Energy ��� EENS ��� gen ���е��
   q = GenUnit(i,3);          
   p = 1-q;                   
   K = GenUnit(i,2)/gcd;
   
   if i==1
   Eg(i) = sum(EF0(1:K))*p;     
   Ed(i) = totalE-Eg(i);   
   index = index + K;
   else
   Eg(i) = sum(EF(index+1:index+K,i-1))*p;     
   Ed(i) = Ed(i-1)-Eg(i);
   index = index + K;
   end
   
end

LOLP= ((EF(Jn,numGen)+EF(Jn+1,numGen)))/(2*time*gcd);
EENS=Ed(numGen);

end

