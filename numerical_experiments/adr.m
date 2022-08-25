function memo = adr( test )
% ADR equation:
%
% Here we solve the 2-dimensional Advection--Diffusion--Reaction (ADR) equation
%
% u_t = epsilons Delta*u - alpha Nabla*u + gamma u(u-1/2)(1-u)
%
% for Omega = [0,1]^2, t in [0,0.25], epsilon = 1/8, alpha=-3/2, gamma = 1000.
% The initial data is u(0) = 256x^2y^2(1-x)^2(1-y)^2, while we set homogeneous
% Dirichlet boundary conditions over { x,y in \partial Omega: xy = 0 } and
% homogeneous Neumann conditions elsewhere.
%
% For the space discretization of this problem we set N_x = N_y = 500
% discretization points for each dimension while for time marching we perform
% N_t = 256, 512, 768, 1024, 1280, 1536, 1792, and 2048 time steps with the
% fourth order Runge-Kutta exponential integrator expRK4s6.
%
% This integrator consists in computing 4 different linear combinations of
% phi-functions: we can compute each one of them with a single call to bamphi.
% The total cost amounts then to 4*N_t calls to the bamphi/kiops routine.
% The exponential integrator expRK4s6 is such that linearity changes
%   x (obviously) never inside any call of bamphi
%   x (obviously) never inside any exponential integration step
%   x never throughout the whole exponential integration process
% therefore with bamphi we could fully integrate this equation in time with
% just ONE call of the Arnoldi process.
% Nevertheless, we found advantageous to compute a dedicated set of Ritz's
% values for each of the four calls to the bamphi routine during the FIRST
% exponential integration step amounting to a grand  total of 4 calls to the
% Arnoldi process.

  % Equation parameters
  t0 = 0.00e+000;
  tf = 0.25e+000;
  omega.x.l = 0; omega.x.r = 1;
  omega.y.l = 0; omega.y.r = 1;
  if isfield( test, 'epsilon')
    epsilon.x = test.epsilon;
    epsilon.y = test.epsilon;
  else
    epsilon.x = 1 / 8;
    epsilon.y = 1 / 8;
  end
  alpha.x   = - 3 / 2;
  alpha.y   = - 3 / 2;
  gamma     = 1.00e+003;

  g  = @( u ) gamma * u .* ( u - 1 / 2 ) .* ( 1 - u );
  dg = @( u ) ( - 3 * gamma ) * ( u.^2 - u + 1 / 6 );

  % Space discretization
  Ns = test.Ns;
  b_c.x.l = 'hom_neumann'; b_c.x.r = 'hom_dirichlet';
  b_c.y.l = 'hom_neumann'; b_c.y.r = 'hom_dirichlet';
  [  D1,  D2, h, X, Y ] = space_discr( omega, Ns, b_c );
  [ DD1, DD2          ] = kronsum_struct( Ns, D1, D2, alpha, epsilon );
  L = DD2 - DD1;
  Lfun = @( x ) L * x;

  % Initial datum
  u0 = 256 * (    ( X(:) - omega.x.l ) .* ( omega.x.r - X(:) ) ...
               .* ( Y(:) - omega.y.l ) .* ( omega.y.r - Y(:) ) ).^2;

  % Time discretization
   Nt = test.Nt;
  tol = min( h.x^2, h.y^2 ) * 1.00e-003;
  jacobian = strncmp( test.integrator, 'epi', 3 ) || strncmp( test.integrator, 'exprb', 5 );
  if not( jacobian )
       linearity = @(   x ) matfun( @( z ) Lfun( z ), x );
    nonlinearity = @( t,x ) g( x );
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
           linearity = @( x ) matfun( @( z ) Lfun( z ), x ) + J .* x;
        nonlinearity = @( x ) matfun( @( z ) Lfun( z ), x ) +  g( x );
      end
      [ uref, info ] = feval( [ test.integrator, '_', test.routines{ 1 } ], uref, k, t, linearity, nonlinearity, opts, info );
      t = t + k;
    end
  end

  % LAUNCH TESTS
  fprintf('\n');
  for rout = 1 : length( test.routines )
    disp(['- Running tests with ', test.routines{ rout } ]);
    for l = 1 : length( Nt )
      opts = []; info = [];
      if     strcmp( test.routines{ rout }, 'bamphi' )
        opts.tol = tol;
        clear epirk4s3a_bamphi
        clear exprk4s6_bamphi
      elseif strcmp( test.routines{ rout }, 'kiops' )
        opts = tol;
        clear epirk4s3a_kiops
        clear exprk4s6_kiops
      end
      clear matfun
      k = ( tf - t0 ) / Nt( l );
      t = 0;
      u = u0;
      tic;
      for iter = 1 : Nt( l )
        if jacobian
          J = dg( u );
             linearity = @( x ) matfun( @( z ) Lfun( z ), x ) + J .* x;
          nonlinearity = @( x ) matfun( @( z ) Lfun( z ), x ) +  g( x );
        end
        [ u, info ] = feval( [ test.integrator, '_', test.routines{ rout } ], u, k, t, linearity, nonlinearity, opts, info );
        t = t + k;
      end % marching
      [ ~, mv ] = matfun( NaN,NaN );
      memo.(test.integrator).(test.routines{ rout }).tstep( l ) = Nt( l );
      memo.(test.integrator).(test.routines{ rout }).clock( l ) = toc;
      memo.(test.integrator).(test.routines{ rout }).matve( l ) = mv;
      if test.compute_error
        memo.(test.integrator).(test.routines{ rout }).error( l ) = norm( u - uref, inf ) / norm( uref, inf );
      end
      if strcmp( test.routines{ rout },'bamphi' )
        memo.(test.integrator).(test.routines{ rout }).info{ l } = info;
      end
      disp(['     ', num2str( l ), ' on ', num2str( length( Nt ) ),' - ', num2str( Nt(l) ),' timesteps... ', num2str( toc ),' seconds.' ]);
    end % timesteps
  end % routines

end
