function [A,s,t] = SBM(ngroup, nvertices, groups, Q)
%generate a graph based on the SBM parameters Q, and the group membership
%of the nodes (groups). nvertices is the number of vertices, ngroups is the
%number of groups

Q = Q/nvertices;
groupvertices = zeros(ngroup, nvertices);

for igroup = 1:ngroup
    i = 1;
    for ivert = 1:nvertices
        if groups(ivert) == igroup
            groupvertices(igroup,i) = ivert;
            i = i+1;
        end
    end
end
groupvertices;
s = [];
t = [];
A = sparse(nvertices,nvertices);

for ia = 1:ngroup
    for ib = ia:ngroup
        ga = groupvertices(ia,groupvertices(ia,:)>0);
        gb = groupvertices(ib,groupvertices(ib,:)>0);
        if ia ~= ib
            for na = 1:length(ga)
                for nb = 1:length(gb)
                    p = Q(ia,ib);
                    y = random('bino',1,p);
                    if y == 1
                        a = ga(na);
                        b = gb(nb);
                        s = [s a];
                        t = [t b];
                        A(a,b) = 1;
                        A(b,a) = 1;
                    end
                end
            end
        
        else
            for na = 1:length(ga)
                for nb = (na+1):length(gb)
                    p = Q(ia,ib);
                    y = random('bino',1,p);
                    if y == 1
                        a = ga(na);
                        b = gb(nb);
                        s = [s a];
                        t = [t b];
                        A(a,b) = 1;
                        A(b,a) = 1;
                    end
                end
            end
        end 
    
    end
end

%plot(graph(A))

end

