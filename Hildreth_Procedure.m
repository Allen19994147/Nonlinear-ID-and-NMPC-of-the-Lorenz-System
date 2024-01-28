%%% This function uses Hildreth Procedure to find Constrained Optimization
%%% Problem
function eta = Hildreth_Procedure(E,F,M,gamma)

% Dimension check:
if(size(M,1)~=size(gamma,1))
    error("Dimensions of M and gamma are inconsistent!");
end
if(size(M,2)~=size(E,2) || size(M,2)~=size(F,1))
    error("Dimensions of decision variables are inconsistent!..." + ...
        "Check E, M and F.");
end
if(size(E,1)~=size(E,2))
    error('E should be a square matrix!')
end




eta=-E\F;
violation_num=0;
for j=1:size(M,1) % each row of M, #constraint
    if (M(j,:)*eta>gamma(j)) % check # violated constraints
        violation_num = violation_num + 1;
        break
    end
end
if (violation_num==0)
    return
else
    H = M*(E\M');
    K = M*(E\F) + gamma;
    [n,m] = size(K);
    lambda = zeros(n,m);
    for i = 1:100 % Enlarge i if lambda converge slowly
        lambda_p = lambda;
        for j=1:n
            w = H(j,:)*lambda-H(j,j)*lambda(j,1);
            w = w + K(j,1);
            lambda(j,1) = max(0,-w/H(j,j));
        end
        al=(lambda-lambda_p)'*(lambda-lambda_p); % if converged
        if (al<10e-8)
            break
        end
    end

    eta= -E\F -E\M'*lambda;
end

