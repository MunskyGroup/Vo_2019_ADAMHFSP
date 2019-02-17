function C = orth_reduce( A, B )
%ORTH_REDUCE C = orth_reduce( A, B ) reduces the matrix B into an
%orthogonal matrix C with columns orthogonal to the subspace spanned by
%columns of A.
%
%Contact:
%=======
%Huy Vo. huydvo@colostate.edu

C = zeros( size(B) );
jc = 0;
for jb = 1:size(B,2)
    u = B(:, jb);
    for it = 1:2 % re-orthogonalization
        for i = 1:size(A,2)
            u = u - (u'*A(:,i))*A(:,i);
        end
        for i = 1:jc
            u = u - (u'*C(:,i))*C(:,i);
        end
    end
    unorm = norm(u);
    if ( unorm > 1.0e-14 ) % only keep non-redundant vectors
        C(:,jc+1) = u/unorm;
        jc = jc+1;
    end
end
end

