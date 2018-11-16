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
%C = [5;4]
%A = [6,4; 
%    1,2;
%   -1,1;
%     0,1]
%b = [24;6;1;2]
%equality = [1;1;1;1]
%Max= 1



%MAIN PROBLEM
C = [2;3;4;1;8;1]
A =[1,-1,2,0,1,1
  0,1,-1,1,0,3;
  1,1,-3,1,1,0;
   1,-1,0,0,1,1]
    
b = [18;8;36;23]
equality = [0,1,1,1]
Max = 1


RevSimMeth(C,A,b,equality,Max);


function RevSimMeth(C, A,b,equality,Max)
clc



[m, n] = size(A);


%Findinng number of slack and artificail variables
sla = sum(equality == 1) + sum(equality == -1)
art = sum(equality == 0) + sum(equality == -1)

Slack = zeros(m,sla);
Artificial = zeros(m,art);

C2 = [];
C3 = [];

%Adding slack/surpul variablies to c = 0 
C2 = zeros(1,sla);


%Addiing artifictial Variables to C if needed = +-M
if art > 0
    M = -1 * Max * 1000;
    for i = 1:art     
        C3(i) = M;        
    end 
end

C = vertcat(C, C2.', C3.')

if Max == -1 
    C = -1 * C
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

Slack
Artificial

%Creating a table of A
SimTab1 = horzcat(A,Slack, Artificial)

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
B


%If there are artificial variables in starting table 
%Satisfying the zero coefficiant rule for basic variables



SimTab2 = horzcat(vertcat(C.',SimTab1),vertcat([0],b))
[s,t] = size(SimTab2);

if art > 0 
     
    for m = 1:art
        n = find(Artificial(:,m) == 1)
        for i = 1:t
            SimTab2(1,i) = SimTab2(1,i) - Max * SimTab2(n+1,i) * M;
        end     
    end
          
    SimTab2
    C = SimTab2(1,[1:t-1]) 
    Msol = SimTab2(1,t)
end





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
    
    enterV = min(enter(enter<0))
    enterP = NonBasVar(find(enter == enterV))
    
             
    end


%Feasibility
    function Feasibility(B, enterP)
        "===FEASIBITY==="
        X_B
        inv(B) * P(enterP)
        exit = (X_B ./ (inv(B) * P(enterP)))
        
        %Stops because all exit <  0
        exitV = min(exit(exit>0))
        
        
        exitP = BasVar(find(exit == exitV)) 
        
        
        
    end



%Next Iteration 
    function [B,C_B] = NextIteration(B, exitP, enterP)
        "===NEXTITERATION==="
        %Swapping Variables
        NonBasVar(NonBasVar == enterP) = exitP
        BasVar(BasVar == exitP) = enterP
        
        %Next C_B 
        for i = 1:length(BasVar)
            C_B(i) = C(BasVar(i));
        end 
        C_B
        
        %Next B Matrix
        for i = 1:length(BasVar)
            B(:,i) = P(BasVar(i));
        end 
        B
        
        
        %Next X_B 
        X_B = inv(B) * b
        
        %Next Z 
        if art > 0 
            Z =  Max * -1 *(Msol - (C_B * X_B))
        else 
            Z = C_B * X_B
        end 
   
    end 

"=====STARTING REVISED SIMPLEX METHOD========"
BasVar
NonBasVar
X_B = b
C_B = zeros(1,length(BasVar))

Optimality(B, C_B)

counter = 0 ;
while Optimised == 0 && counter < 10
Feasibility(B, enterP)
[B,C_B] = NextIteration(B, exitP, enterP)
Optimality(B, C_B)
counter = counter + 1;
end 


"===SOLUTION==="

BasVar 
C_B
X_B
Z




end 
