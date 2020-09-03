function [LDC] = TrapezoidLoad(peak,StepLoad)
% Simulate Trapezoid Load Profile base load at 50% peak load
% col1 MW, col2 freq, col3 indi prob., col4 cum. prob., col5 duration

i=1; T=8760;
maxrow=(0.5*peak/StepLoad)+1;                       % base load at 50% peak load
LDC=zeros(maxrow,5);

while i <= maxrow
    LDC(i,1)= (peak/2)+(i-1)*StepLoad;
    LDC(i,5)= (LDC(i,1)-peak)*8760/(-0.5*peak);     % Duration    
    LDC(i,4)=  LDC(i,5)/T;                          % Cumulative probability
    i=i+1;
end
 
for i=1:maxrow
    k = maxrow-i+1;
    if i == 1
        LDC(k,3)= 0;                                % individual probability
        LDC(k,2)= LDC(k,3)*T;                       % frequency
    else
        LDC(k,3)= LDC(k,4)-LDC(k+1,4);
        LDC(k,2)= LDC(k,3)*T;
    end        
end    

end

