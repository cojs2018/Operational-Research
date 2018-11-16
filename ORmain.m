function ORmain
    c = [2;3;4;1;8;1;0;0;0;0];
    A = [1,-1,2,0,1,1,1,0,0,0;
        0,1,-1,1,0,3,0,1,0,0;
        1,1,-3,1,1,0,0,0,1,0;
        1,-1,0,0,1,1,0,0,0,1];
    
    b = [18;8;36;23];
    
    str = "max";
    
    [n, m] = size(A); %Find the size of A
    
    slack = zeros(n);
    for i=(n-1):-1:1
        slack(n-i) = m-i; %Get indecies of the slack variables
    end
    
    SOL = 0;
    
    c_B = c((m-n+1):m, :); %Cut c to the lenth of b
    B = A(:, (m-n+1):m); %Cut A to the length of b
    
    x = revised_simplex(c, A, b, c_B, B, slack, str, SOL);
end

function x = revised_simplex(c, A, b, c_B, B, slack, str, SOL)
    if str == "max"
        %Function solves the equation z = c^T * x subject to A*x = b.
        [n, m] = size(A); %Find the size of A

        %Work out the entering P variable (optimal)
        P_enter = 0;
        for j = 1:(m-n)
            if -1*c(j,:) < P_enter_coeff 
                P_enter = j; %entering variable
            end
        end

        %Works out leaving P variable (feasible)
        L = zeros(n);
        V = A(:, P_enter)/B;
        for i = 1:n
           L(i) = b(i) / V(i);
        end
        P_leave = 0;
        for j = 1:n
           if (L(j) == min(L)) 
               P_leave = j;
           end
        end

        c_B(P_leave,:) = c(P_enter,:);
        B(:, P_leave) = A(:, P_enter);

        x_B = b/B; 

        z = c_B' * x_B; %Work out value of z = c_B^T * x_B

        if z >= SOL
            x = revised_simplex(c, A, b, c_B, B, slack, str, z);
        else
            x = revised_simplex(c, A, b, c_B, B, slack, str, SOL);
        end
    elseif str == "min"
        %Function solves the equation z = c^T * x subject to A*x = b.
        [n, m] = size(A); %Find the size of A

        %Work out the entering P variable (optimal)
        P_enter = 0;
        for j = 1:(m-n)
            if c(j,:) > P_enter_coeff 
                P_enter = j; %entering variable
            end
        end

        %Works out leaving P variable (feasible)
        L = zeros(n);
        V = A(:, P_enter)/B;
        for i = 1:n
           L(i) = b(i) / V(i);
        end
        P_leave = 0;
        for j = 1:n
           if (L(j) == min(L)) 
               P_leave = j;
           end
        end

        c_B(P_leave,:) = c(P_enter,:);
        B(:, P_leave) = A(:, P_enter);

        x_B = b/B; 

        z = c_B' * x_B; %Work out value of z = c_B^T * x_B

        if z <= SOL
            x = revised_simplex(c, A, b, c_B, B, slack, str, z);
        else
            x = revised_simplex(c, A, b, c_B, B, slack, str, SOL);
        end
    else
        error("Recall function with str set to either 'min' or 'max'"); %Return error
    end
end

function P = get_P(A)
    P = A;
end

function simplex_out(z)
    
end