%%% A library of candidate functions for SINDy
%%% You may extend this library to fit your use.

% Note that the state, X, and input, U are decoupled in this library,
% you may couple them if necessary.
function Theta = Candidate_Library(X,U)
[~,n] = size(X); % #state
[~,m] = size(U); % #input

% Build library
Theta = [X,U];

% 2nd order
for i=1:n
    for j=i:n
        Theta = [Theta X(:,i).*X(:,j)];
    end
end
for i=1:m
    for j=i:m
        Theta = [Theta U(:,i).*U(:,j)];
    end
end
% 3rd order
for i=1:n
    for j=i:n
        for k=j:n
            Theta = [Theta X(:,i).*X(:,j).*X(:,k)];
        end
    end
end
for i=1:m
    for j=i:m
        for k=j:m
            Theta = [Theta U(:,i).*U(:,j).*U(:,k)];
        end
    end
end
% 4rd order
for i=1:n
    for j=i:n
        for k=j:n
            for ll=k:n
                Theta = [Theta X(:,i).*X(:,j).*X(:,k).*X(:,ll)];
            end
        end
    end
end
for i=1:m
    for j=i:m
        for k=j:m
            for ll=k:m
                Theta = [Theta U(:,i).*U(:,j).*U(:,k).*U(:,ll)];
            end
        end
    end
end

% sin & cos, freq up to w
for w = 1:5
    for i=1:n
        Theta = [Theta sin(w.*X(:,i)) cos(w.*X(:,i))];
    end
end
for w = 1:5
    for i=1:m
        Theta = [Theta sin(w.*U(:,i)) cos(w.*U(:,i))];
    end
end

