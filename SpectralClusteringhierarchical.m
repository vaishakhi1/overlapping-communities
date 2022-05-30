%Implement a hierarchical 
clc
clear

% enter the number of vertices, groups. For now, the code is hard coded to support 2 groups with overlapping nodes.

prompt = 'Enter the number of vertices in the graph:';
beep
nvertices = input(prompt);
ngroup = 2;
ogroup = ngroup^2 - 1;                                                      % List of all the possible overlapping groups

%To generalize the code, change from here in future versions
smatrix = 0:0.1:1;                                                          % fraction of nodes that overlap
xmax = nvertices/10;                                                        % average degree of the nonoverlapping graph - the lower this value is set, the sparser the graph is
cmatrix = 0:xmax/10:xmax;                                                   % Range of separation (c_{in} - c_{out})
errormatrix = [];
% 2 loops, one for s, one for x
results_temp = zeros(length(smatrix),length(cmatrix));
results = results_temp;
niteration = 5;                                                             % Number of iterations
tic
bound = ksbound(smatrix,xmax)                                               % compute the KS bound
for iiter = 1:niteration                                                    % iteration loop
    for is = 1:length(smatrix)                                              % overlap loop
        overlap = smatrix(is);                                              % fraction of ovrelap
        
        %assign groups to the nodes
        %uniformly at random with some nodes overlapping - change soon
        p_group = [0.5 0.5];%(1/ngroup)*ones(1,ngroup);
        %groups = random('unid', ngroup,1,nvertices);
        for ix = 1:length(cmatrix)                                          % separation loop
            
            
            groups = random('bino',1,p_group(1),1,nvertices)+1;             % assign groups based on p_group
            
            % separate the overlapping nodes - set aside about s/2 of the
            % nodes in each group to be a part of the other group 
            for i = 1:ngroup
                ga1 = find(groups == i);
                ovlap = random('bino',1,overlap,length(ga1),1);
                ga1 = ovlap'.*ga1;
                ovg = ga1(find(ga1>0));
                groups(ovg) = 3;
            end
            
            
            x = cmatrix(ix);                                                % separation
            cin = (xmax/2)+x/2;                                             % in-degree
            cout = (xmax/2)-x/2;                                            % out-degree
            
            % Affinity matrix of the non-overlapping graph
            Q_simple = cout*ones(ngroup, ngroup) + (cin - cout)*eye(ngroup);
            
            if overlap== 0
                Q = Q_simple;
                ogroup = ngroup;
            else
                % set up affinity matrix for overlap sets, separated into disjoint sets
                Z = [1 0; 0 1; 1 1];   % this is hardcoded in - how to change?
                ogroup = ngroup^2-1;
                % affinity matrix of the overlapping graph
                Q = Z*Q_simple*Z';
            end
            
            
            % A check to ensure the affinity matrix is computed correctly
            % for the the number of groups
            check = sum(size(Q) == [ogroup, ogroup]);
            if (check ~= 2)
                disp('the affinity matrix is wrong')
                beep
                return;
            end
            
            %generate the SBM
            [A,s,t] = SBM(ogroup, nvertices, groups,Q);
            G = graph(A);
            
            
            % classification
            if (ogroup == 2)
                
                % classify without degree, since the average degrees of all
                % nodes should be similar
                g2 = clsfy(G,1);
                
                if (sum((g2==groups')) <= nvertices/2)
                    results_temp(is, ix) = nvertices - sum((g2==groups'));
                else
                    results_temp(is, ix) = sum((g2==groups'));
                end
                
            else
                
                % First classify the nodes based on the degree of the nodes
                g2 = clsfy(G,0);
                
                % identify the overlapping nodes - by looking at the set of
                % nodes with higher average degree
                
                nodelist1 = find(g2 == 1);
                mean1  = mean(degree(G,nodelist1));
                
                nodelist2 = find(g2 == 2);
                mean2  = mean(degree(G,nodelist2));
                
                % check the group assignment
                c2 = 10*g2+groups';
                cf = returncount(c2);
                                
                [g_hier, min_mean] = min([mean1,mean2]);                    %separate low degree nodes
                edgelist = find(g2 == min_mean);
                G1 = subgraph(G,edgelist);                                  % separate the nonoverlapping nodes
                g3 = clsfy(G1,1);                                           % classify the low degree nodes
                
                % some adjustments to assign groups correctly
                if min_mean == 1
                    [~ , max_mean] = max([mean1 mean2]);
                    g_max = find(g2 == 2);
                    g2(g_max) = 1;
                end
                g2(edgelist) = g3+1;
                y = g2;
                
                % check errors
                c = 10*y+groups';
                cf = returncount(c);
                groupassign = [1:ogroup;zeros(1,ogroup)];
                used = [];
                b = 1:ogroup;
                err = [];
                
                % compute the maximum overlap and compare against the
                % original group membership
                while (isempty(c) == 0)
                    a = mode(c);
                    q = fix(a/10);
                    r = mod(a,10);
                    [i2, err1] = missclassified(ogroup, q,r,cf);
                    groupassign(2,q) = r;
                    used = [used r];
                    ind1 = find(mod(c,10) == r);
                    c(ind1) = [];
                    ind1 = find(fix(c/10) == q);
                    c(ind1) = [];
                    err = [err;[i2, err1]];
                end
                
                err = sortrows(err,1)
                a = find(groupassign(2,:)==0);
                k = ismember([1 2 3],used);
                k = [1 1 1] - k;
                k1 = find((k.*[1 2 3])>0);
                for i = 1:length(k1)
                    groupassign(2,a(i)) = k1(i);
                end
                groupassign;
                groupassign(1,:)*10+groupassign(2,:);
                y_reassigned = zeros(size(y));
                for i = 1:ogroup
                    ind1 = find(y==i);
                    y_reassigned(ind1) = groupassign(2,i);
                end
                results_temp(is, ix) = sum((y_reassigned==groups'));
                %             if(size(err,1) == 3)
                %              errormatrix = [errormatrix, err(:,2)];
                %             end
                
            end
            
            
        end
    end
    iiter
    toc
    results = results + results_temp;
end
results = (niteration)^-1*results;

results = results/nvertices;

figure

hold on

plot(cmatrix, results(1,:),'-*b','LineWidth',2,'MarkerSize',10);
plot(cmatrix, results(4,:),'-or','LineWidth',2,'MarkerSize',10);
plot(cmatrix, results(6,:),'-^g','LineWidth',2,'MarkerSize',10);
plot(cmatrix, results(11,:),'->m','LineWidth',2,'MarkerSize',10);

%
plot(bound([1 4 6]),0.3*[1 1 1], 'o','MarkerSize',10)

plot([bound(1) bound(1)], [0.3, 1],'--b')
plot([bound(4) bound(4)], [0.3, 1],'--r')
plot([bound(6) bound(6)], [0.3, 1],'--g')
%
% %
% % line([17.885, 17.885], [0,0.7])
% %
legend('s = 0', 's = 0.3','s = 0.5','s = 1')
xlabel ( 'Cin - Cout', 'FontSize', 16)
ylabel('Accuracy','FontSize', 16)
