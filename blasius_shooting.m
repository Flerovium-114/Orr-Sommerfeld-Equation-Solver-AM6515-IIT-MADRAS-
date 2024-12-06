clear all;
close all;
clc;

%% INPUTS
h = 0.1;   % eta step size - use this to control accuracy of calculations

% Velocity Profile Inputs
xi = linspace(0,1,50);
use_Blasius = 1; % use Blasius velocity profile below
transform_u = 1; % transform velocity to xi coordinates

% Temporal Stability Inputs
run_temporal_analysis = 1 ; % 1 = run temporal analysis, 0 = skip
t_Re_min = 1e3; % minimum Re
t_Re_max = 1e5; % maximum Re
t_N_Re = 40;    % number of Re
t_Res = logspace(log10(t_Re_min), log10(t_Re_max), t_N_Re);

t_alpha_min = 0.001;    % minimum alpha for temporal analysis
t_alpha_max = 1.3;      % maximum alpha
t_N_alpha = 40;         % number of alpha to use
t_alphas = linspace(t_alpha_min,t_alpha_max,t_N_alpha);

% for m1 = 1:length(s_alphais)
%     for m2 = 1:length(s_alphars)
%         alpha = s_alphars(m2)+1i*s_alphais(m1);
%         s_alphas(m1,m2) = alpha;
%     end
% end


%% DEFINED FIGURE NUMBERS (do not correspond to report)
FIG_VEL_PROFILE = 1;
FIG_TEM_EIG = 2;
FIG_TEM_CONT = 3;
FIG_SPA_EIG = 4;
FIG_SPA_MAP = 5;
FIG_SPA_CONT = 6;

%% COMPUTE BLASIUS VELOCITY PROFILE
if ( use_Blasius )
    fprintf('Solving for Blasius Velocity Profile ... ');
    tic; % start counting elapsed time

    % Blasius ODE Definition:  f''' + ff'' = 0  ->  f''' = -ff''
    fun = @(eta,f) [ ...
        f(2); ...
        f(3); ...
        -f(1)*f(3); ...
        ];

    eta0 = 0;   % Initial eta station (@ wall )
    etaf = 10;  % Final eta station (@ "freestream")
    target = 1; % Freestream BC: f'(inf) = U(inf) = 1
    tol = 1e-6; % Shooting method tolerance

    f0 = [ 0; 0; .45 ]; % Initial condition for f, f', f'' at wall

    go = 1;
    while ( go == 1 ) % loop until converged
        % Standard 4th-order Runge-Kutta method
        [ eta, f ] = rk4( fun, eta0, etaf, f0, h );

        % compute error, adjust f''(0)
        err = f(2,end) - target;
        if ( abs(err) < tol )
            go = 0;
        else
            f0(3,1) = f0(3,1) * ( target-err );
        end
    end
    U = f(2,:);             % U   = f'(eta)
    Up = f(3,:);            % U'  = f''(eta)
    Upp = -f(1,:).*f(3,:);  % U'' = f'''(eta)
    elapsed_time = toc;     % compute time elapsed
    fprintf('Done (%f sec)\n',elapsed_time);
end

%% TRANSFORM ETA TO XI
if ( transform_u )
    fprintf('Transforming Coordinates ... ');
    tic; % start counting elapsed time

    % find delta in terms of eta
    minval = 1; % initialize minimization variable
    index = 0;  % index of the edge of the Boundary Layer
    delta_vel = 0.999;  % percent of freestream velocity

    % find edge of Boundary layer, where U == delta_vel
    for n = 1:length(U)
        if ( abs(U(n)-delta_vel) < minval )
            index = n;
            minval = abs(U(n)-delta_vel);
        end
    end

    delta = eta(index); % Boundary layer thickness
    xi = eta/delta;     % define xi s.t. xi = 1 where eta = delta

    % transform to xi coordinate system
    % U = U;            % U(xi)   = U(eta)
    Up = Up*delta;      % U'(xi)  = U'(eta)*delta
    Upp = Upp*delta^2;  % U''(xi) = U''(eta)*delta^2

    elapsed_time = toc; % compute elapsed time
    fprintf('Done (%f sec)\n',elapsed_time);
end
NMAX = length(xi);  % look at all eigenvalues
% NMAX = index;       % only look at eigenvalues in B.L.

%% TEMPORAL STABILITY ANALYSIS
if ( run_temporal_analysis )
    fprintf('Temporal Stability Analysis ...\n');
    tic; % start counting elapsed time

    h = xi(2)-xi(1);    % xi step size
    for l = 1:length(t_Res)
        Re = t_Res(l);
        for m = 1:length(t_alphas)
            alpha = t_alphas(m);
            fprintf('\ta = %.2f, Re = %d ... ',alpha,Re);

            % generate the A matrix
            A = zeros(NMAX); % initialize A

            n = 2; % start at n = 2, go to n = nmax-1

            % coefficients as defined in my paper
            a1 = -U(n)*alpha^3 - Upp(n)*alpha + 1i*alpha^4/Re;
            a2 = U(n)*alpha - 2*1i*alpha^2/Re;
            a3 = 1i/Re;

            % 2nd order Central Difference Method
            A(n,n-1) = a2 - 4*a3/h^2;
            A(n,n)   = a1*h^2 - 2*a2 + 6*a3/h^2;
            A(n,n+1) = a2 - 4*a3/h^2;
            A(n,n+2) = a3/h^2;

            for n = 3:NMAX-2
                % coefficients as defined in my paper
                a1 = -U(n)*alpha^3 - Upp(n)*alpha + 1i*alpha^4/Re;
                a2 = U(n)*alpha - 2*1i*alpha^2/Re;
                a3 = 1i/Re;
                % 2nd order Central Difference Method
                A(n,n-2) = a3/h^2;
                A(n,n-1) = a2 - 4*a3/h^2;
                A(n,n)   = a1*h^2 - 2*a2 + 6*a3/h^2;
                A(n,n+1) = a2 - 4*a3/h^2;
                A(n,n+2) = a3/h^2;
            end

            n = NMAX-1;
            a1 = -U(n)*alpha^3 - Upp(n)*alpha + 1i*alpha^4/Re;
            a2 = U(n)*alpha - 2*1i*alpha^2/Re;
            a3 = 1i/Re;
            % 2nd order Central Difference Method
            A(n,n-2) = a3/h^2;
            A(n,n-1) = a2 - 4*a3/h^2;
            A(n,n)   = a1*h^2 - 2*a2 + 6*a3/h^2;
            A(n,n+1) = a2 - 4*a3/h^2;

            A = A(2:end-1,2:end-1); % remove first/last rows/cols of A (0 BCs)

            % generate B
            B = zeros(NMAX);    % initialize B
            for n = 2:NMAX-1
                % coefficients as defined in my paper
                b1 = -alpha^3;
                b2 = alpha;
                % 2nd order Central Difference Method
                B(n,n-1) = b2;
                B(n,n)   = b1*h^2 - 2*b2;
                B(n,n+1) = b2;
            end
            B = B(2:end-1,2:end-1); % remove first/last rows/cols of B (0 BCs)


            [V,e] = eig(B\A);   % invert B, compute eigenvalues of LHS
            e = diag(e);    % turn diagonal eigenvalue matrix into vector

            % store most unstable eigenvalue
            [maxe,cindex] = max(imag(e));
            fprintf('c = %f + %fi\n',real(e(cindex)),imag(e(cindex)));
            t_c(m,l) = e(cindex);

            % Plot eigenvalues of the discrete system
            figure(FIG_TEM_EIG)
            plot(real(e),imag(e),'.','MarkerSize',6)
            xlabel('c_r')
            ylabel('c_i')
            title(sprintf('alpha = %.2f, Re = %d',alpha,Re));
            axis tight

        end
    end
    elapsed_time = toc; % compute elapsed time
    %fprintf('Done (%f sec)\n',elapsed_time);
end

%% OUTPUT RESULTS
if ( run_temporal_analysis )
    cilevels = [0,.005,.01,.015,.0019];
    figure(FIG_TEM_CONT);
    contour(t_Res,t_alphas,imag(t_c),cilevels);
    set(gca,'XScale','log');
    xlabel('$Re_\delta$','Interpreter','latex','Fontsize',14);
    ylabel('$\bar{\alpha}$','Interpreter', 'latex','Fontsize',14);
    title('Neutral Stability Curves')
    colorbar
    grid on;
end

% uprop = ['k- ';'k: ';'k--'];
uprop = ['b';'r';'g'];
figure(FIG_VEL_PROFILE);
plot(U,xi,uprop(1,:),Up,xi,uprop(2,:),Upp,xi,uprop(3,:))
title('Blasius Velocity Profile')
xlabel('U, U'', U''''')
ylabel('\xi = y/\delta')%,'Interpreter','latex','FontSize',16)
legend('U','U''','U''''','location','best')
grid on;
