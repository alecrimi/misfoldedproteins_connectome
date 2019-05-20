
list = dir('*.csv');

for kk = 1 : length(list)
list(kk).name
    A = csvread(list(kk).name);

    D = diag(sum(A));
    L = D - A ;
    beta = 0.005;

    X0 = zeros(size(A));
    % 8,13,97,182 alpha syn
    % 24,25,30,32, 28.35 for abeta, 28 and 35 are parahippocampal?

    X0(24,:) = A(24,:);
    X0(25,:) = A(25,:);
    X0(30,:) = A(30,:);
    X0(32,:) = A(32,:);
    X0(28,:) = A(28,:);
    X0(35,:) = A(35,:);

    X1 = exp(-L*0.0005).*X0;
 
 for k = 1 : 10
    
    [r,c] = find(X1);
    for jj = 1 : length(c)
        X0(c(jj),:) = A(c(jj),:);
    end
    
    X1 =  exp(-L*0.0005*k).*X0 + X1;

 end

 csvwrite(['sim_' list(kk).name],sum(X1));
end 
