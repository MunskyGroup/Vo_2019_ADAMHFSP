function [wt, err, hump] = expvt( T, A, v, tol, m )
%  [wt, err, hump] = expvt( T, A, v, tol, m )
%  EXPVT computes an approximation of expv(t*A)*v where t ranges over the
%  values in T. This is a modification of the function expv in
%  R.B.Sidje's EXPOKIT to output intermediate solution vectors. See the
%  documentation of EXPOKIT for more information about the algorithm.
%
%  Arguments:
%  ==========
%
%  T : one-dimensional vector of time points to output the solutions.
%
%  A : the matrix.
%
%  v : initial vector.
%
%  tol : tolerance.
%
%  m : dimension of the Krylov basis.
%
%  wt (out) : n x nt array of intermediate solutions, n == length(v) and nt
%  == length(T).
%
%  err : sum of local errors.
%
%  hump : an estimation of the matrix exponential 'hump' (see the reference
%  for details).
%  
%  Reference:
%  ==========
%  Roger B. Sidje 
%  EXPOKIT: Software Package for Computing Matrix Exponentials. 
%  ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
%
%  Contact (about the modification):
%  ================================
%  Huy Vo. huydvo@colostate.edu

[n,n] = size(A);
if nargin == 3,
  tol = 1.0e-7;
  m = min(n,30);
end;
if nargin == 4,
  m = min(n,30);
end;

anorm = norm(A,'inf'); 
mxrej = 10000;  btol  = 1.0e-7; 
gamma = 0.9; delta = 1.2; 
mb    = m; t_out   = T(end);
nstep = 0; t_new   = 0;
t_now = 0; s_error = 0;
rndoff= anorm*eps;

k1 = 2; xm = 1/m; normv = norm(v); beta = normv;
fact = (((m+1)/exp(1))^(m+1))*sqrt(2*pi*(m+1));
t_new = (1/anorm)*((fact*tol)/(4*beta*anorm))^xm;
s = 10^(floor(log10(t_new))-1); t_new = ceil(t_new/s)*s; 
sgn = sign(t_out); nstep = 0;

nt = length(T);
wt = zeros(length(v),nt);

w = v;
hump = normv;
while t_now < t_out
  nstep = nstep + 1;
  t_step = min( t_out - t_now,t_new );
  V = zeros(n,m+1); 
  H = zeros(m+2,m+2);

  V(:,1) = (1/beta)*w;
  for j = 1:m
     p = A*V(:,j);
     for i = 1:j
        H(i,j) = V(:,i)'*p;
        p = p-H(i,j)*V(:,i);
     end;
     s = norm(p); 
     if s < btol,
        k1 = 0;
        mb = j;
        t_step = t_out - t_now;
%         disp('Happy breakdown.');
        break;
     end;
     H(j+1,j) = s;
     V(:,j+1) = (1/s)*p;
  end; 
  if k1 ~= 0, 
     H(m+2,m+1) = 1;
     avnorm = norm(A*V(:,m+1)); 
  end;
  ireject = 0;
  while ireject <= mxrej,
     mx = mb + k1;
     F = expm(sgn*t_step*H(1:mx,1:mx));
     if k1 == 0,
        err_loc = btol; 
        break;
     else
        phi1 = abs( beta*F(m+1,1) );
        phi2 = abs( beta*F(m+2,1) * avnorm );
        if phi1 > 10*phi2,
           err_loc = phi2;
           xm = 1/m;
        elseif phi1 > phi2,
           err_loc = (phi1*phi2)/(phi1-phi2);
           xm = 1/m;
        else
           err_loc = phi1;
           xm = 1/(m-1);
        end;
     end;
     if err_loc <= delta * t_step*tol,        
        break;
     else
        t_step = gamma * t_step * (t_step*tol/err_loc)^xm;
        s = 10^(floor(log10(t_step))-1);
        t_step = ceil(t_step/s) * s;
        if ireject == mxrej,
           error('The requested tolerance is too high.');
        end;
        ireject = ireject + 1;
     end;
  end;
  mx = mb + max( 0,k1-1 );
  w = V(:,1:mx)*(beta*F(1:mx,1));


  % Huy Vo's modification
  % find the output time points that we have passed with the current Krylov basis
  Ipassed = find( t_now <= T & T <= t_now+t_step ); 
  
  if ( ~isempty(Ipassed) )
      for ip = 1:length(Ipassed)   
        it = Ipassed(ip);
        texport = T(it) - t_now;
        F = expm(texport*H(1:mx,1:mx));
        wt(:,it) = V(:,1:mx)*(beta*F(1:mx,1));        
      end
  end
  % End of Huy Vo's modification
  
  beta = norm( w );
  hump = max(hump,beta);  
  
  t_now = t_now + t_step;
  t_new = gamma * t_step * (t_step*tol/err_loc)^xm;
  s = 10^(floor(log10(t_new))-1); 
  t_new = ceil(t_new/s) * s;  

  err_loc = max(err_loc,rndoff);
  s_error = s_error + err_loc;
end;
err = s_error;
hump = hump / normv;




