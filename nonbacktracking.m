function [B Bmarg] = nonbacktracking(s,t)
%compute the nonbacktracking matrix from edges (s,t)
m = length(s);
n = max(max(s),max(t));
B = sparse(2*m,2*m);
sa = [s t];
ta = [t s];
Bmarg = sparse(n,2*m);
for i = 1:2*m
    for j = 1:2*m
        if (ta(i) == sa(j)) && (ta(j) ~= sa(i))
            B(i,j) = 1;
            Bmarg(sa(i),i) = 1;
        end
    end 
end

B;
end

