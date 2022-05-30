function y = clsfy(G,mrkr)
% Clasify the nodes into two groups based on the variable mrkr
% 0 - separate high degree ands low degree nodes 
% 1 - groups membership - based on node neighbours


A = adjacency(G);
D = diag(sum(A));
%         L = D-A;
%         C = diag(D);
%         E = diag(C.^-0.5);
%         E(find(E==inf)) = 0;
%         F = full(E);
%         %normalized graph laplacian
%         L_norm = F*full(L)*F;
ogroup = 2;
nvertices = length(A);
%equivalent nonbacktracking matrix (easier to compute)
Bdash = [zeros(nvertices) D - speye(nvertices); -1*speye(nvertices) A];

%[B1,Bmarg] = nonbacktracking(s,t);                                          %to compute the actual nonbacktracking matrix
%B = Bmarg*B1;                                                               % marginalize

%with Bdash

[Vdash,~] = eigs(Bdash,ogroup,'lm');                                         % compute eigenvalues
%Vdash = Bmarg*Vdash;

Vdash = Vdash(nvertices+1:2*nvertices,:);                                    % select half the eigenvalues

% to plot the eigenvalues
% Beig = eigs(B1,200);
%  plot(real(Beig),imag(Beig),'*')


if mrkr == 0
    % separate based on degree
    y = kmeans([real(Vdash(:,1:ogroup)), degree(G,1:numnodes(G))'],ogroup);
else
    % separate based on group membership only
    y = kmeans(real(Vdash(:,1:ogroup)),ogroup);
end

%             end

end
