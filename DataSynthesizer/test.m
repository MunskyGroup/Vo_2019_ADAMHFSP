% Test the data generation on a toggle swtich model
%% Define the model
nmax = [50 50];

x0 = [31 5];
Pfull0 = zeros(prod(nmax+1),1); i0 = sub2ind(nmax+1, x0(1)+1, x0(2)+1); Pfull0(i0) = 1;

S = [1 0;1 0;-1 0; 0 1; 0 1; 0 -1]'; % stoichiometry matrix

ayx = 6.1e-3; axy = 2.6e-3; nyx = 2.1; nxy = 3; % we assume these parameters are known

ind_prop = @(x,y) [ones(length(x),1), 1./(1+ayx*(y.^nyx)), x, ones(length(y),1), 1./(1+axy*(x.^nxy)), y]; % the parameter-independent part of the propensities

params0 = [5e-5, 1e-2, 5e-4, 2e-3, 1e-2, 3e-4];

t_array = linspace( 1, 10000, 4 );

prop = @(x) params0.*ind_prop(x(1), x(2));
datcell = data_gen( S', prop, x0, t_array, 1000, [] );
%% Solve the FSP to compare with the generated data
Acell = fsp_get_matrices( S, ind_prop, nmax );
Pcell = fsp_solve( params0, t_array, Acell, Pfull0 );

%% plot and compare results

figure();
for it = 1:length(t_array)
    subplot( length(t_array), 2, 2*(it-1)+1);
    contourf( reshape(Pcell{it}, nmax(1)+1, nmax(2)+1 ) );
    subplot( length(t_array), 2, 2*(it-1)+2);
    xdat = datcell{it}(:, 1:end-1);
    wdat = datcell{it}(:,end);
    contourf( full( sparse( xdat(:,1)+1, xdat(:,2)+1, wdat ) ) );
    xlim ([0 nmax(2)+1] );
    ylim ( [0 nmax(1)+1] );
end