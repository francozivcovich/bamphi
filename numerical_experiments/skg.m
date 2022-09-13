function memo = skg( test )
% Sine Klein-Gordon equation:
%
% Here we solve the 2-dimensional sine Klein-Gordon equation
%
% u_tt = Delta * u + sin( u )
%
% for Omega = [0,1]^2, t \in [0,1] and initial data u_t(0) = 256x^2y^2(1-x)^2 (1-y)^2
% u(0) = 0. We set homogeneous Dirichlet boundary conditions. For the space
% discretization of this problem we set 500 discretization points for each
% dimension while for time marching we perform N_t = 16, 32, 64, 128, 256, 512,
% 1024, and 2048 time steps with the second order trigonometric integrator exptg2s1.
%

  mode = test.mode;

  % Equation parameters
  t0        = 0.00e+000;
  tf        = 1.00e+000;
  omega.x.l = 0; omega.x.r = 1;
  omega.y.l = 0; omega.y.r = 1;
  epsilon.x = 1;
  epsilon.y = 1;
    alpha.x = 1;
    alpha.y = 1;

  % Space discretization
  Ns = test.Ns;
  N = prod( Ns );
  g = @( u ) [ zeros( N, 1 ); sin( u( 1 : N ) ) ];
  b_c.x.l = 'hom_dirichlet'; b_c.x.r = 'hom_dirichlet';
  b_c.y.l = 'hom_dirichlet'; b_c.y.r = 'hom_dirichlet';
  [  D1,  D2, h, X, Y ] = space_discr( omega, Ns, b_c );
  [ DD1, DD2          ] = kronsum_struct( Ns, D1, D2, alpha, epsilon );
  L = DD2;

  [ R, ~, rs ] = chol( - L, 'vector' ); Rt = R';

  function x = Lfun  ( x )
    x       = L  * x;
  end

  function x = Rfun  ( x )
    x       = R  * x( rs );
  end

  function x = Rtfun ( x )
    x( rs ) = Rt * x;
  end

  function x = iRfun ( x )
    x( rs ) = R  \ x;
  end

  function x = iRtfun( x )
    x       = Rt \ x( rs );
  end

  function x = RRS( x )
    x = [ x( N + 1 : 2 * N ); Lfun( x( 1 : N ) ) ];
  end

  function x = LRLR( x )
    x = [ Rtfun( x( N + 1 : 2 * N ) ); - Rfun( x( 1 : N ) ) ];
  end


  % Initial datum
  u0 = [ 0 * X(:);
         256 * (    ( X(:) - omega.x.l ) .* ( omega.x.r - X(:) ) ...
                 .* ( Y(:) - omega.y.l ) .* ( omega.y.r - Y(:) ) ).^2 ];
  switch mode
    case {'AA'}
         linearity = @(   x ) matfun( @( z ) RRS( z ), x );
      nonlinearity = @( t,x ) g( x );
    case {'RR'}
         u0 = [ u0( 1 : N ); iRtfun( u0( N + 1 : 2 * N ) ) ];
         linearity = @(   x ) matfun( @( z ) LRLR( z ), x );
      nonlinearity = @( t,x ) [ zeros( N, 1 ); iRtfun( sin( x( 1 : N ) ) ) ];
    case {'RA'}
         linearity = @(   x ) matfun( @( z ) RRS( z ), x );
         aux_linrt = @(   x ) matfun( @( z ) LRLR( z ), x );
      nonlinearity = @( t,x ) g( x );
           aux_fun = @(   x ) matfun( @( z ) iRtfun( z ), x );
  end

  % Time discretization
  Nt = test.Nt;
  tol = min( h.x^2, h.y^2 ) * 1.00e-003;
  jacobian = strncmp( test.integrator, 'epi', 3 ) || strncmp( test.integrator, 'exprb', 5 );

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
      if strcmp( test.mode,'RA' ) && strcmp( test.routines{ 1 },'bamphi' )
        [ uref, info ] = feval( [ test.integrator, '_', test.routines{ 1 } ], uref, k, t, linearity, nonlinearity, opts, info, aux_linrt, aux_fun );
      else
        [ uref, info ] = feval( [ test.integrator, '_', test.routines{ 1 } ], uref, k, t, linearity, nonlinearity, opts, info );
      end
      t = t + k;
    end
    if strcmp( mode,'RR' )
      uref = [ uref( 1 : N ); Rtfun( uref( N + 1 : 2 * N ) ) ];
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
        if strcmp( test.mode,'RA' ) && strcmp( test.routines{ rout },'bamphi' )
          [ u, info ] = feval( [ test.integrator, '_', test.routines{ rout } ], u, k, t, linearity, nonlinearity, opts, info, aux_linrt, aux_fun );
        else
          [ u, info ] = feval( [ test.integrator, '_', test.routines{ rout } ], u, k, t, linearity, nonlinearity, opts, info );
        end
        t = t + k;
      end % marching
      if strcmp( mode,'RR' )
        u = [ u( 1 : N ); Rtfun( u( N + 1 : 2 * N ) ) ];
      end
      [ ~, mv ] = matfun( NaN,NaN );
      memo.(test.integrator).(test.routines{ rout }).tstep( l ) = Nt( l );
      memo.(test.integrator).(test.routines{ rout }).clock( l ) = toc;
      memo.(test.integrator).(test.routines{ rout }).matve( l ) = mv;
      if test.compute_error
        err = norm( u - uref, inf ) / norm( uref, inf );
        memo.(test.integrator).(test.routines{ rout }).error( l ) = err;
      else
        err = NaN;
      end
      if strcmp( test.routines{ rout },'bamphi' )
        memo.(test.integrator).(test.routines{ rout }).info{ l } = info;
      end
      disp(['     ', num2str( l ), ' on ', num2str( length( Nt ) ),' - ', num2str( Nt(l) ),' timesteps... ', num2str( toc ),' seconds.' ]);
    end % timesteps
  end % routines

end
