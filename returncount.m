function cf = returncount(c)             

% in the case of two groups (3 nonoverlapping groups) check the group
% assignments to characterize the incorrectly assigned nodes

ctemp = [c ; 33];
                cf1 = tabulate(ctemp);
                %hist(y*10+groups',33)
                
            
                
                
                
                
                cset = [11 length(find(c == 11))
                    12 length(find(c == 12))
                    13 length(find(c == 13))
                    21 length(find(c == 21))
                    22 length(find(c == 22))
                    23 length(find(c == 23))
                    31 length(find(c == 31))
                    32 length(find(c == 32))
                    33 length(find(c == 33))];
                
                cf = [cf1(cset(:,1),1) cf1(cset(:,1),2)];
                cf(9,2) = cf(9,2)-1;
end