clc 
clear

%TEST PG19 OF NOTES 
C = [4;1];
A = [1,2;
    3,1;
    4,3];
b= [4;3;6]  ;
equality = [1,0,-1];
Max = -1;


%TEST WITH ONLY <= ON PG 30 of NOTES
% C = [5;4]
% A = [6,4; 
%     1,2;
%    -1,1;
%     0,1]
% b = [24;6;1;2]
% equality = [1;1;1;1]
% Max= 1



%MAIN PROBLEM
C = [2;3;4;1;8;1];
 A =[1,-1,2,0,1,1
   0,1,-1,1,0,3;
   1,1,-3,1,1,0;
    1,-1,0,0,1,1];
     
 b = [18;8;36;23];
 equality = [0,1,1,1];
 Max = 1;

 %pg26 infeasible
%  C = [3;2]
%  A = [2,1;
%       3,4]
%  b = [2;12]
%  equality = [1, -1]
%  Max = 1
 
 %pg23 degeneracy
%  C = [3;9]
%  A = [1,4;
%       1,2]
%  b = [8;4]
%  equality = [1, 1]
%  Max = 1

 
 %pg23 Alternate optima
%  C = [2;4]
%  A = [1,2;
%       1,1]
%  b = [5;4]
%  equality = [1, 1]
%  Max = 1
 
 %pg24 unbounded
%  C = [2;1]
%  A = [1,-1;
%       2,0]
%  b = [10;40]
%  equality = [1, 1]
%  Max = 1
%  
RevSimMeth(C,A,b,equality,Max);


clf

%Graphs
X = 0.1:0.25:30;
Y = [];

subplot(4,2,1)
hold on
for i = 1:length(X)
   tmp = RevSimMeth(C,A,[X(i);8  ;36 ;23],equality,Max);
   tmp = tmp(7);
   Y(i) = tmp;
end 
plot(X,Y)
plot([18 18],[0 232],'black')
plot([0 18],[232 232],'black')
title("Graph of how the optimal solution changes with resource 1")
ylabel("Optimal solution")
xlabel("Resouce 1")


subplot(4,2,2)
hold on
for j = 1:6
for i = 1:length(X)
   tmp = RevSimMeth(C,A,[X(i);8  ;36 ;23],equality,Max);
   tmp = tmp(j);
   Y(i) = tmp;
end 
plot(X,Y)
end
plot([18 18],[0 40],'black')
legend("x1","x2","x3","x4","x5","x6", 'Location','northeastoutside')
title("Graph of how the values of x change with resource 1")
ylabel("x")
xlabel("Resouce 1")


X = 0.1:0.25:30;
Y = [];

subplot(4,2,3)
hold on
for i = 1:length(X)
   tmp = RevSimMeth(C,A,[18; X(i);36 ;23],equality,Max);
   tmp = tmp(7);
   Y(i) = tmp;
end 
plot(X,Y)
plot([8 8],[0 232],'black')
plot([0 8],[232 232],'black')
title("Graph of how the optimal solution changes with resource 2")
ylabel("Optimal solution")
xlabel("Resouce 2")


subplot(4,2,4)
hold on
for j = 1:6
for i = 1:length(X)
   tmp = RevSimMeth(C,A,[18; X(i)  ;36 ;23],equality,Max);
   tmp = tmp(j);
   Y(i) = tmp;
end 
plot(X,Y)
end
plot([8 8],[0 60],'black')
legend("x1","x2","x3","x4","x5","x6", 'Location','northeastoutside')
title("Graph of how the values of x change with resource 2")
ylabel("x")
xlabel("Resouce 2")


X = 0.1:0.25:50;
Y = [];

subplot(4,2,5)
hold on
for i = 1:length(X)
   tmp = RevSimMeth(C,A,[18; 8; X(i) ;23],equality,Max);
   tmp = tmp(7);
   Y(i) = tmp;
end 
plot(X,Y)
plot([36 36],[200 232],'black')
plot([0 36],[232 232],'black')
title("Graph of how the optimal solution changes with resource 3")
ylabel("Optimal solution")
xlabel("Resouce 3")


subplot(4,2,6)
hold on
for j = 1:6
for i = 1:length(X)
   tmp = RevSimMeth(C,A,[18; 8  ;X(i) ;23],equality,Max);
   tmp = tmp(j);
   Y(i) = tmp;
end 
plot(X,Y)
end
plot([36 36],[0 30],'black')
legend("x1","x2","x3","x4","x5","x6", 'Location','northeastoutside')
title("Graph of how the values of x change with resource 3")
ylabel("x")
xlabel("Resouce 3")

X = 0.1:0.25:50;
Y = [];

subplot(4,2,7)
hold on
for i = 1:length(X)
   tmp = RevSimMeth(C,A,[18; 8; 36 ; X(i)],equality,Max);
   tmp = tmp(7);
   Y(i) = tmp;
end 
plot(X,Y)
plot([23 23],[200 232],'black')
plot([0 23],[232 232],'black')
title("Graph of how the optimal solution changes with resource 4")
ylabel("Optimal solution")
xlabel("Resouce 4")


subplot(4,2,8)
hold on
for j = 1:6
for i = 1:length(X)
   tmp = RevSimMeth(C,A,[18; 8  ;36 ; X(i)],equality,Max);
   tmp = tmp(j);
   Y(i) = tmp;
end 
plot(X,Y)
end
plot([36 36],[0 30],'black')
legend("x1","x2","x3","x4","x5","x6", 'Location','northeastoutside')
title("Graph of how the values of x change with resource 4")
ylabel("x")
xlabel("Resouce 4")


hold off 

RevSimMeth(C,A,b,equality,Max)


function  solution = RevSimMeth(C, A,b,equality,Max)
clc


[m, ~] = size(A);


%Findinng number of slack and artificail variables
sla = sum(equality == 1) + sum(equality == -1);
art = sum(equality == 0) + sum(equality == -1);

Slack = zeros(m,sla);
Artificial = zeros(m,art);

C2 = [];
C3 = [];

%Adding slack/surpul variablies to c = 0 
C2 = zeros(1,sla);


%Addiing artifictial Variables to C if needed = +-M
if art > 0
    M = -1 * Max * 100000;
    for i = 1:art     
        C3(i) = M;        
    end 
end

C = vertcat(C, C2.', C3.');

if Max == -1 
    C = -1 * C;
end

%Adding Slack, Surplus and Artificial Variables to A 
countS = 1;
countA = 1;
for i = 1:m
    if equality(i) == -1
        Slack(i,countS) = -1;
        countS = countS + 1;
        Artificial(i,countA) = 1;
        countA = countA + 1 ;
    end   
    
    if equality(i) == 0
        Artificial(i,countA) = 1;
        countA = countA + 1  ;
    end
    
    if equality(i) == 1
        Slack(i,countS) = 1;
        countS = countS + 1;
    end
end 

Slack;
Artificial;

%Creating a table of A
SimTab1 = horzcat(A,Slack, Artificial);


%function for jth column of matrix A 
function Pj = P(j)       
     Pj = SimTab1(:,j);
end


[~, q] = size(SimTab1);
%finding basic and non basic variables
countBas = 1;
countNonBas = 1;

for j = 1:q
    if sum(P(j)) == 1 && nnz(P(j)) == 1
        BasVar(countBas) = j;
        countBas = countBas + 1;
    else 
        NonBasVar(countNonBas) = j;
        countNonBas = countNonBas + 1;
    end
end

BasVar

NonBasVar


%Starting B
for i = 1:length(BasVar)
    B(:,i) = P(BasVar(i));
end 
B;


%If there are artificial variables in starting table 
%Satisfying the zero coefficiant rule for basic variables



SimTab2 = horzcat(vertcat(C.',SimTab1),vertcat([0],b))
[s,t] = size(SimTab2);

if art > 0 
     
    for m = 1:art
        n = find(Artificial(:,m) == 1);
        for i = 1:t
            SimTab2(1,i) = SimTab2(1,i) - Max * SimTab2(n+1,i) * M;
        end     
    end
          
    SimTab2
    C = SimTab2(1,[1:t-1]) ;
    Msol = SimTab2(1,t);
end


CounterOpt = 1;
CounterFea = 1;


%OPTIMIALITY 
    function Optimality(B, C_B)
    "===OPTIMALITY==="
    
    
    C_B * inv(B);
    enter = [];
        for  j = 1:length(NonBasVar) 
            enter(j) = C_B * inv(B) * P(NonBasVar(j)) - C(NonBasVar(j));          
        end
    enter  
        %Testing if optimised
        if all(enter >= 0)
            Optimised = 1
            return
        else Optimised = 0;
        end  
    
    enterV = min(enter(enter<=0));
    enterP = NonBasVar(find(enter == enterV))
    
    
    if length(enterP) == 1
        CounterOpt = 1 ;
    else 
        enterP = enterP(CounterOpt)
        CounterOpt = CounterOpt + 1;
    end 
    
    
    end


%Feasibility
    function Feasibility(B, enterP)
        "===FEASIBITY==="
        X_B;
        inv(B) * P(enterP);
        exit = (X_B ./ (inv(B) * P(enterP)))
        
        %Stops because all exit <  0
        exitV = min(exit(exit>=0));
        
        
        exitP = BasVar(find(exit == exitV)) 
        

        
        if length(exitP) == 1
            CounterFea = 1 ;
        else 
            exitP = exitP(CounterFea)
        CounterFea = CounterFea + 1;
        end
        
    end



%Next Iteration 
    function [B,C_B] = NextIteration(B, exitP, enterP)
        "===NEXTITERATION==="
        %Swapping Variables
        NonBasVar(NonBasVar == enterP) = exitP;
        BasVar(BasVar == exitP) = enterP;
        
        %Next C_B 
        for i = 1:length(BasVar)
            C_B(i) = C(BasVar(i));
        end 
        C_B;
        
        %Next B Matrix
        for i = 1:length(BasVar)
            B(:,i) = P(BasVar(i));
        end 
        B
        
        
        %Next X_B 
        X_B = inv(B) * b;
        
        %Next Z 
        if art > 0 
            Z =  Max * -1 *(Msol - (C_B * X_B));
        else 
            Z = C_B * X_B;
        end 
   
    end 

"=====STARTING REVISED SIMPLEX METHOD========"
BasVar;
NonBasVar;
X_B = b;
C_B = zeros(1,length(BasVar));

Optimality(B, C_B)

counter = 0 ;
while Optimised == 0 && counter < 10
Feasibility(B, enterP)
[B,C_B] = NextIteration(B, exitP, enterP);
Optimality(B, C_B)
counter = counter + 1;
end 


"===SOLUTION==="
BasVar ;
C_B;
X_B;


solution = zeros(1,length(NonBasVar));
for i = find(BasVar <= length(NonBasVar))
solution(BasVar(i)) =  X_B(i);
end


Z;
solution(length(NonBasVar)+ 1) = Z


"From online calculator Optimal Solution: z = 232; x1 = 0, x2 = 8, x3 = 0, x4 = 0, x5 = 26, x6 = 0"


end 


