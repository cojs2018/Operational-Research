function opt=OperationalResearch(A,xb,z) % A, xb, z are matrices
    m=1;
    A2=A;
    while m~=0
        A3=A2;
        B=zeros(length(xb));
        for i=1:length(xb)
            for j=1:length(xb)
                B(j,i)=A2(j,length(A2)-length(xb)+i);
            end
        end
        c=zeros(length(xb),1);
        for i=1:length(xb)
            c(i,1)=z(1,length(A2)-length(xb)+i);
        end
        S1=size(A2);
        Z=zeros(1,length(A2)-length(B));
        P=zeros(S1(2),length(A2)-length(B));
        for i=1:length(A2)-length(B)
            Z(1,i)=z(1,i);
            for j=1:S1(2)
                P(j,i)=A2(j,i);
            end
        end
        CB=(c'/B)*P-Z;
        m=min(CB);                          %largest negative value
        E=1;
        while CB(i)~=m
            E=E+1;                %entering variable
        end
        for i=1:length(xb)
            P2=A2(i,E);
        end
        BP=Bi*P2;
        ch=zeros(1,S1(2));
        for i=1:S1(2)
            ch(i)=xb(i)*BP(i);
        end
        for i=1:length(xb)         %smallest positive value
            if ch(i)==-abs(ch(i))
                ch(i)=100;
            else
                ch(i)=abs(ch(i));
            end
        end
        m2=min(ch);
        L=1;
        while ch(i)~=m2
            L=L+1;             %leaving variable
        end
        for i=1:length(xb)       %swapping columns
            A2(i,L)=A3(i,E);
            A2(i,E)=A3(i,L);
        end
    end
    disp('the optimised solution is found')
    opt=zeros(length(xb));
    for i=1:length(xb)
        for j=1:length(xb)
            opt(j,i)=A3(j,length(A3)-length(xb)+i);
        end
    end
    disp(opt)
end
A=[1,0,1,1;-1,1,1,-1;2,-1,-3,0;0,1,1,0;1,0,1,1;1,3,0,1;1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1]
xb=[18,8,36,23]
z=[2;3;4;1;8;1;0;0;0;0]
OperationalReasearch(A,xb,z)
    
        
        
            