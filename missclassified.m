function [ind, err] = missclassified(ogroup, q,r,cf)
%Just to count the errors

b = 1:ogroup;
act = find(b == r);
oth = find(b ~= q);
indices = 10*oth+act*ones(size(oth));
ind = r;
err = 0;
for i = 1:size(oth,2)
    j = find(cf(:,1) == indices(i)) ;
    err = err + cf(j,2);
end
