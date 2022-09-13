function memo = brg( test )
% BURGERS's equation:
%
% Here we solve the 2-dimensional Burger's equation
%
% u_t + u_x^2 / 2 + u_y^2 / 2 = epsilon \Delta*u,
%
% for Omega = [-5/4,5/4]^2, t in [0,2], and epsilon = 4/100. The initial data
% is u(0) = exp(-10\sin^2(x/2)-10\sin^2(y/2)), while we set homogeneous
% Dirichlet boundary conditions over {(x,y) in \partial Omega: y <= -x} and
% homogeneous Neumann conditions elsewhere.
% For the space discretization of this problem we set 500 discretization points
% for each dimension while for time marching we perform N_t = 8, 16, 24, 32, 56,
% and 64 time steps with the fourth order exponential integrator EPIRK4s3A.
%
% This integrator consists in computing 2 different linear combinations of
% phi-functions: we can compute each one of them with a single call to bamphi.
% The total cost amounts then to 2*N_t calls to the bamphi/kiops routine.
%
% The exponential integrator EPIRK4s3A is such that linearity changes
%   x (obviously) never inside any call of bamphi
%   x (obviously) never inside any exponential integration step
%   x at each exponential integration step
% therefore with bamphi we could fully integrate this equation in time with N_t
% calls to the Arnoldi process.
% Nevertheless, we found advantageous to compute a dedicated set of Ritz's
% values for each of the two calls to the bamphi routine amounting to a grand
% total of 2*N_t calls to the Arnoldi process.
%
% The number r of (IOM)-Arnoldi steps, which is initially set to its maximum,
% 64, adaptively changes after each call of bamphi for efficiency.
% To keep the tests simple, we did not explore mixed strategies such as
% launching the (IOM)-Arnoldi procedure at each stage but not at every
% exponential integration step.
% Such strategies might be rewarding if the solution is in a slowly varying
% region or in those situations where practitioners use to refresh the
% Jacobian matrix not at each time step but once in a while.

  % Equation parameters
  t0 = 0.00e+000;
  tf = 2.00e+000;
  omega.x.l = - 1.25; omega.x.r = 1.25;
  omega.y.l = - 1.25; omega.y.r = 1.25;
  epsilon.x = 4.00e-002;
  epsilon.y = 4.00e-002;
  alpha.x   = 1.00e+000;
  alpha.y   = 1.00e+000;

  % Space discretization
  Ns = test.Ns; N = prod( Ns );
  b_c.x.l = 'hom_dirichlet'; b_c.x.r = 'hom_neumann';
  b_c.y.l = 'hom_dirichlet'; b_c.y.r = 'hom_neumann';
  [  D1,  D2, h, X, Y ] = space_discr( omega, Ns, b_c );
  [ DD1, DD2          ] = kronsum_struct( Ns, D1, D2, alpha, epsilon );
  L = DD2;
  Lfun = @( x ) L * x;

  g  = @( u ) - DD1 * ( u.^2 / 2 );
  dg = @( u ) - DD1 * spdiags( u, 0, N, N );

  % Initial datum
  u0 = exp( - 10 * sin( X / 2 ).^2 - 10 * sin( Y / 2 ).^2 );
  u0 = u0(:);

  % Time discretization
   Nt = test.Nt;
  tol = min( h.x^2, h.y^2 ) * 1.00e-003;
  jacobian = strncmp( test.integrator, 'epi', 3 ) || strncmp( test.integrator, 'exprb', 5 );
  if not( jacobian )
       linearity = @(   x ) matfun( @( z ) Lfun( z ), x );
    nonlinearity = @( t,x ) g( x );
  else
    nonlinearity = @( x ) matfun( @( z ) Lfun( z ), x ) + g( x );
  end

  % Compute reference
  if test.compute_error
    Ntref = 2 * Nt( end );
    opts = []; info = [];
    if     strcmp( test.routines{ 1 }, 'bamphi' )
      opts.tol = tol;
    elseif strcmp( test.routines{ 1 }, 'kiops'  )
      opts = tol;
    end
    k = ( tf - t0 ) / Ntref;
    t = 0;
    uref = u0;
    fprintf(['- Compute reference with ', num2str( Ntref ),' timesteps of ', test.integrator, ' in combination with ', test.routines{ 1 },': %3.0f%%'], 100 * 1 / Ntref );
    for iter = 1 : Ntref
      fprintf('\b\b\b\b%3.0f%%', 100 * iter / Ntref );
      if jacobian
        J = dg( uref );
        linearity = @( x ) matfun( @( z ) Lfun( z ), x ) + J * x;
      end
      [ uref, info ] = feval( [ test.integrator, '_', test.routines{ 1 } ], uref, k, t, linearity, nonlinearity, opts, info );
      t = t + k;
    end
  end

  % LAUNCH TESTS
  fprintf('\n');
  for rout = 1 : length( test.routines )
    disp(['- Running tests with ', test.routines{ rout } ]);
    aux = NaN;
    for l = 1 : length( Nt )
      opts = []; info = []; clear_integrators()
      if     strcmp( test.routines{ rout }, 'bamphi' )
        opts.tol = tol; % opts.skew = 1;
      elseif strcmp( test.routines{ rout }, 'kiops' )
        opts = tol;
      end
      clear matfun
      k = ( tf - t0 ) / Nt( l );
      t = 0;
      u = u0;
      tic;
      for iter = 1 : Nt( l )
        if jacobian
          J = dg( u );
          linearity = @( x ) matfun( @( z ) Lfun( z ), x ) + J * x;
        end
        [ u, info ] = feval( [ test.integrator, '_', test.routines{ rout } ], u, k, t, linearity, nonlinearity, opts, info );
        t = t + k;
      end % marching
      [ ~, mv ] = matfun( NaN,NaN );
      memo.(test.integrator).(test.routines{ rout }).tstep( l ) = Nt( l );
      memo.(test.integrator).(test.routines{ rout }).clock( l ) = toc;
      memo.(test.integrator).(test.routines{ rout }).matve( l ) = mv;
      if test.compute_error
        err = norm( u - uref, inf ) / norm( uref, inf );
        memo.(test.integrator).(test.routines{ rout }).error( l ) = err;
      else
        err = norm( u - aux, inf ) / norm( aux, inf );
        aux = u;
      end
      if strcmp( test.routines{ rout },'bamphi' )
        memo.(test.integrator).(test.routines{ rout }).info{ l } = info;
      end
      disp(['     ', num2str( l ), ' on ', num2str( length( Nt ) ),' - ', num2str( Nt(l) ),' timesteps... ', num2str( toc ),' seconds.' ]);
    end % timesteps
  end % routines

end
