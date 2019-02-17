function [P_out] = rom_solve(parameters, t_array, PHIcell, t_out, reduced_matrices)
% rom_solve solve the reduced order model and at the same time
% calculate the estimated error.

M = size(reduced_matrices.B, 1);
nt = length(t_array);

tau = t_array(1);
qstart = reduced_matrices.q0;
iexport= ismember( t_array, t_out );
itex= 1;
P_out= cell(1, length(t_out));

for jt = 1:nt
    B_now = parameters(1)*reduced_matrices.B{1, jt};
    for k = 2:M
        B_now = B_now + parameters(k)*reduced_matrices.B{k, jt};
    end           
    
    if (size(B_now,1) >= 300)
        qend = expv( tau, B_now, qstart, 1e-10, 30);    
    else 
        qend = expm(tau*B_now)*qstart;
    end
    
    if ( iexport(jt) )
        P_out{itex}= PHIcell{jt}*qend;
        itex= itex+1;
    end
    
    % Project the reduced solution of current time interval to the subspace of the
    % next time interval
    if (jt<nt)
        tau = t_array(jt+1) - t_array(jt);
        qstart = reduced_matrices.T{jt}*qend;
    end    
end

end

