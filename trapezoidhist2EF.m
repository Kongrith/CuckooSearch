function [EF] = trapezoidhist2EF(x,N,gcd)
% calculate energy function ด้วยโหลดแบบ histogram สมมติฐานเป็นสี่เหลี่ยมคางหมู
Load = x(:,1); duration = x(:,5);
CumArea = zeros(N,1); E = zeros(N,1);
for k=1:N
    J=N-k+1;
    baseline=gcd*(J-1);
    I_tmp = find(Load>baseline);      % เพื่อหา index ของ load block ทีมีค่ามากกว่า step load ที่สนใจ
    maxrow = size(I_tmp,1);           % หาจำนวนแถวของ block load 
    BlockLoadArea = zeros(maxrow,1);   
    for i=1:maxrow                    
    % hint การหาพื้นที่จะใช้วิธีคือใช้จุดที่ต่ำกว่า1จุดมาช่วยในการหาพื้นที่ โดยจะทำการ check ว่าจุด ต่ำกว่า criteria หรือเปล่า
        i1 = I_tmp(i);  i2 = i1-1;                      
        if i1 ~= 1                                 % ทำ poka-yoke โดยจะไม่คิดจุดที่ 1
            x1 = Load(i1);     x2 = Load(i2);
            y1 = duration(i1); y2 = duration(i2);
            criteria = Load(i1-1);        
            if baseline > criteria          
                % กรณีนี้จะคิด พท.สี่เหลี่ยมคางหมูแบบไม่เต็ม
                %disp('partial-trapeziod area')
                unknown = y1+(y1-y2)*(baseline-x1)/(x1-x2);
                b = abs(unknown-y1);
                h = abs(x1-baseline);
                BlockLoadArea(i) = 0.5*b*h; 
            else
                % baseline <= criteria กรณีนี้จะคิด พท.สี่เหลี่ยมคางหมู
                %disp('trapezoid area')
                BlockLoadArea(i) = 0.5*abs(y1-y2)*((x1-baseline)+(x2-baseline));
            end
        end
    end    
    CumArea(J) = sum(BlockLoadArea);        % พื้นที่ใต้กราฟสะสมแต่ละ energy function ใน block ที่ J    
    if k==1
        E(J)= CumArea(J);
    else          
        E(J)= CumArea(J)-CumArea(J+1);      % หาพื้นที่ใต้กราฟแต่ละ energy function ใน block ที่ J
    end              
end    

EF=E;
totalE = sum(EF); 

end

