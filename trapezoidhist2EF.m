function [EF] = trapezoidhist2EF(x,N,gcd)
% calculate energy function ������ŴẺ histogram ����԰ҹ�������������ҧ���
Load = x(:,1); duration = x(:,5);
CumArea = zeros(N,1); E = zeros(N,1);
for k=1:N
    J=N-k+1;
    baseline=gcd*(J-1);
    I_tmp = find(Load>baseline);      % ������ index �ͧ load block ���դ���ҡ���� step load ���ʹ�
    maxrow = size(I_tmp,1);           % �Ҩӹǹ�Ǣͧ block load 
    BlockLoadArea = zeros(maxrow,1);   
    for i=1:maxrow                    
    % hint ����Ҿ�鹷������Ըդ����ش����ӡ���1�ش�Ҫ���㹡���Ҿ�鹷�� �¨зӡ�� check ��Ҩش ��ӡ��� criteria ��������
        i1 = I_tmp(i);  i2 = i1-1;                      
        if i1 ~= 1                                 % �� poka-yoke �¨����Դ�ش��� 1
            x1 = Load(i1);     x2 = Load(i2);
            y1 = duration(i1); y2 = duration(i2);
            criteria = Load(i1-1);        
            if baseline > criteria          
                % �óչ��ФԴ ��.�����������ҧ���Ẻ������
                %disp('partial-trapeziod area')
                unknown = y1+(y1-y2)*(baseline-x1)/(x1-x2);
                b = abs(unknown-y1);
                h = abs(x1-baseline);
                BlockLoadArea(i) = 0.5*b*h; 
            else
                % baseline <= criteria �óչ��ФԴ ��.�����������ҧ���
                %disp('trapezoid area')
                BlockLoadArea(i) = 0.5*abs(y1-y2)*((x1-baseline)+(x2-baseline));
            end
        end
    end    
    CumArea(J) = sum(BlockLoadArea);        % ��鹷�����ҿ�������� energy function � block ��� J    
    if k==1
        E(J)= CumArea(J);
    else          
        E(J)= CumArea(J)-CumArea(J+1);      % �Ҿ�鹷�����ҿ���� energy function � block ��� J
    end              
end    

EF=E;
totalE = sum(EF); 

end

