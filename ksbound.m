function bound = ksbound(smatrix,xmax) 
%compute the separaion threshold (c_in - c_out) for the overlap given by
%smatrix and the average degree given by xmax
bound = zeros(length(smatrix),1);

tic
for is = 1:length(smatrix)
    s = smatrix(is);
    
    bound(is) = (2/(1-s))*sqrt(xmax*(3*s + 1)/2);
        
end

