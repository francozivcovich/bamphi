function memo = sch_S( test )
% Smooth cubic Schrodinger equation:
%
% Here we solve the 2-dimensional cubic Schrodinger equation
%
% i * u_t = - Delta * u + mu |u|^2 u
%
% with Omega = [0,1]^2, t \in [0,0.1], mu = 1. The initial data is
% u(0) = 256x^2y^2(1-x)^2(1-y)^2 and we set homogeneous Dirichlet boundary
% conditions. For the space discretization of this problem we set N_x = N_y = 200
% discretization points for each dimension while for time marching we perform
% N_t = 2, 4, 8, 16, 32, 64, 128, 256, and 512 time steps with the second order
% strang splitting integrator schss68.
%

  % Equation parameters
  t0 = 0.00e+000;
  tf = 0.10e+000;
  omega.x.l =  0; omega.x.r = 1;
  omega.y.l =  0; omega.y.r = 1;
  epsilon.x = 1;
  epsilon.y = 1;
  alpha.x   = 1;
  alpha.y   = 1;
  mu = 1;

  g  = @( u ) mu;

  % Space discretization
  Ns = test.Ns;
  b_c.x.l = 'hom_dirichlet'; b_c.x.r = 'hom_dirichlet';
  b_c.y.l = 'hom_dirichlet'; b_c.y.r = 'hom_dirichlet';
  [  D1,  D2, h, X, Y ] = space_discr( omega, Ns, b_c );
  [ DD1, DD2          ] = kronsum_struct( Ns, D1, D2, alpha, epsilon );
  L = DD2;
  Lfun = @( x ) L * x;

  % Initial datum
  u0 = ( ( X(:) - omega.x.l ) .* ( omega.x.r - X(:) ) .* ( Y(:) - omega.y.l ) .* ( omega.y.r - Y(:) ) ).^2;
  u0 = 256 * u0;

  % Time discretization
   Nt = test.Nt;
  tol = min( h.x^2, h.y^2 ) * 1.00e-004;
     linearity = @(   x ) matfun( @( z ) Lfun( z ), x );
  nonlinearity = @( t,x ) g( x );

  % Compute reference
  if test.compute_error
    Ntref = 8 * Nt( end );
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
      [ uref, info ] = feval( [ test.integrator, '_', test.routines{ 1 } ], uref, k, t, linearity, nonlinearity, opts, info );
      t = t + k;
    end
  end

  % LAUNCH TESTS
  fprintf('\n');
  for rout = 1 : length( test.routines )
    disp(['- Running tests with ', test.routines{ rout } ]);
    for l = 1 : length( Nt )
      opts = []; info = []; clear_integrators()
      if     strcmp( test.routines{ rout }, 'bamphi' )
        opts.tol = tol;
      elseif strcmp( test.routines{ rout }, 'kiops' )
        opts = tol;
      end
      clear matfun
      k = ( tf - t0 ) / Nt( l );
      t = 0;
      u = u0;
      tic;
      for iter = 1 : Nt( l )
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
      end
      if strcmp( test.routines{ rout },'bamphi' )
        memo.(test.integrator).(test.routines{ rout }).info{ l } = info;
      end
      disp(['     ', num2str( l ), ' on ', num2str( length( Nt ) ),' - ', num2str( Nt(l) ),' timesteps... ', num2str( toc ),' seconds.' ]);
    end % timesteps
  end % routines

end
