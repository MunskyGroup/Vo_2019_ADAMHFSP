function y = loglikelihood_rom(data_cell, P_cell, nmax, unobserved)

threshold = 1.0e-14;

nt = length(P_cell);
y = 0;
for t = 1:nt
    y = y + loglikelihood_1point(data_cell{t}, P_cell{t}, nmax, unobserved);
end
function y = loglikelihood_1point(data, P, nmax, unobserved)
% compute the log-likelihood of the distribution w on the hyper-rectangular
% state space specified in FSP, given the discrete data.
N = length(nmax);
if (~isempty(unobserved))
    % sum along the unobservable dimensions
    pdim = unique([unobserved (1:N)], 'stable');
    sizes = num2cell(nmax+1);
    tensor = reshape(P, sizes{1:N});
    P1 = permute(tensor, pdim);
    nrow = prod(nmax(unobserved)+1);
    P1 = reshape(P1, nrow, numel(P1)/nrow);
    P1 = sum(P1,1)';
else
    P1 = P;
end

nmax1 = nmax(~ismember((1:N), unobserved));
xdat = data(:,1:end-1);
fdat = data(:,end);
Nobs = size(xdat,2);
arg = cell(1,Nobs);

if (Nobs>1)
    for i = 1:Nobs
        arg{i} = xdat(:,i) + 1;
    end
    ind = sub2ind(nmax1+1, arg{1:end});
else
    ind = xdat+1;
end
w1 = P1;
w1(w1<threshold) = threshold;
y = fdat'*log(w1(ind));
end
end

