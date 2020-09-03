function [Eg,LOLP,EENS] = EEF(EF0,GenUnit,gcd,N,time)
% คำนวนหา LOLP และ EENS ด้วยวิธี equivalent energy function method (EEF)
% คำอธิบาย: N คือ block สูงสุงของ load , time คือ period of study , EF0 คือ energy function

numGen = size(GenUnit(:,1),1);       % จำนวนครั้งที่ต้องทำ convolution
Jn = sum(GenUnit(:,2))/gcd;          % หาจำนวน block ที่ต้อง shift จาก gen unit     

Eg=zeros(numGen,1);                  % สร้าง vector สำหรับเก็บค่า energy ที่ gen ผลิตได้
Ed=zeros(numGen,1);                  % สร้าง vector สำหรับเก็บค่า expected energy not served
EF = zeros(Jn+N,numGen);             % สร้าง เมตริกซ์สำหรับ energy function
EF_tmp = zeros(Jn+N,numGen);         % สร้าง เมตริกซ์สำหรับ energy function ที่ shift ไป K

EF0 = vertcat(EF0,zeros(Jn,1));      % เติมค่าเปล่าใส่ใน เมตริกซ EF0 (energy function)
totalE = sum(EF0);                   % พื้นที่ใต้กราฟของ Load ทั้งหมด
EENS = 0; index = 0;

for i = 1:numGen
   q = GenUnit(i,3);                 % ค่า FOR ของ Generator unit
   p = 1-q;                          % ค่า availability ของ Generator unit
   K = GenUnit(i,2)/gcd;             % จำนวน block ที่เราสนใจ
   index = index + K;                %
   
   if i==1                           % ทำการ shift energy function ไปทีละ K บล๊อค
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
for i=1:numGen                       % คำนวนหาค่า Energy และ EENS ที่ gen แต่ละตัว
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

