function ORmain
    c = [2;3;4;1;8;1;0;0;0;0];
    A = [1,-1,2,0,1,1,1,0,0,0;
        0,1,-1,1,0,3,0,1,0,0;
        1,1,-3,1,1,0,0,0,1,0;
        1,-1,0,0,1,1,0,0,0,1];
    
    b = [18;8;36;23];
    
    str = "max";
    
    x = revised_simplex(c, A, b, str);
end

function x = revised_simplex(c, A, b, str)
%Function solves the equation z = c^T * x subject to A*x = b.
    n = length(b); %Let n be the length of b
    
    c_B = c(1:n, :); %Cut c to the lenth of b
    B = A(:, 1:n); %Cut A to the length of b
    Binv = inv(B); %Find B^-1
    
    x_B = Binv * b; %Find basic values of x subject to x_B = B^-1 * b
    
    z = c_B' * x_B; %Work out value of z = c_B^T * x_B
    BinvA = Binv * A; %Work out B^-1 * A
    d = (c_B' * BinvA) - c'; %Calculate (c_B^T)*(B^-1 * A) - c^T
    
    for i = 1:n
        
    end
    
    %Output value in simplex table via new function
    simplex_out(z);
        
    x = 0; %ToDo: Find way of returning data
    
    if str == "max"
        
    elseif str == "min"
        x = 0; %ToDo: Find way of minimising a problem
    else
        error("Recall function with str set to either 'min' or 'max'"); %Return error
    end
end

function simplex_out(z)
    
end