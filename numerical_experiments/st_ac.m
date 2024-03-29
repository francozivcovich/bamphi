function memo = st_ac( test )
% Stochastic Allen-Cahn equation:
%
% Here we solve the 2-dimensional stochastic extension of the Allen--Cahn
% equation
%
% du = ( epsilon Delta * u + u - u^3 ) * dt + sin( u ) * dW
%
% for Omega = [-1,1]^2, t \in [0,0.025], and epsilon = 0.1. The initial data is
% u(0) = sin( 2 pi x )sin( 2 pi y ), the noise term W(x,t) is a Q-Wiener process
% defined on a filtered probability space that is white in time. The noise can
% be represented as a series in the eigenfunctions of the covariance operator
% Q given by
%
% W(x,t) = sum_{i \in \mathbb{N}^2} \sqrt{q_i} e_i(x) beta_i(t)
%
% where (q_i,e_i), i \in \mathbb{N}^2 are the eigenvalues and eigenfunctions of
% the covariance operator $Q$ and $\beta_i$ are independent and identically
% distributed standard Brownian motions. We set homogeneous Dirichlet boundary
% conditions. For the space discretization of this problem we set N_x = N_y = 200
% discretization points for each dimension while for time marching we perform
% N_t = 2, 4, 8, 16, 32, 64, 128, and 256 time steps with the first order
% stochastic exponential integrator SETDM1, in combination with a Monte Carlo
% approach consisting of L = 100 launches.
%

  % Equation parameters
  t0 = 0.00e+000;
  tf = 0.125;%0.25e-001 * 10;
  omega.x.l =  - 1; omega.x.r = 1;
  omega.y.l =  - 1; omega.y.r = 1;
  epsilon.x = 0.25; %0.1;
  epsilon.y = 0.25; %0.1;
  alpha.x   = 1;
  alpha.y   = 1;

  g = @( u ) u - u.^3;
  w = @( u ) sin( u );

  % Space discretization
  Ns = test.Ns;
  b_c.x.l = 'hom_dirichlet'; b_c.x.r = 'hom_dirichlet';
  b_c.y.l = 'hom_dirichlet'; b_c.y.r = 'hom_dirichlet';
  [  D1,  D2, h, X, Y ] = space_discr( omega, Ns, b_c );
  [ DD1, DD2          ] = kronsum_struct( Ns, D1, D2, alpha, epsilon );
  L = DD2;

  size_u = [ prod( Ns ), test.MoCa ];
  Lfun = @( x ) reshape( L * reshape( x, size_u ), [], 1 );

  % Initial datum
  u0 = sin( 2 * pi * X(:) ) .* sin( 2 * pi * Y(:) );
  u0 = u0 * ones( 1, test.MoCa );

  % Time discretization
   Nt = test.Nt;
  tol = min( h.x^2, h.y^2 ) * 1.00e-003;
  jacobian = strncmp( test.integrator, 'epi', 3 ) || strncmp( test.integrator, 'exprb', 5 );
  if not( jacobian )
       linearity   = @(   x ) matfun( @( z ) Lfun( z ), x );
       aux_linrt   = @(   x ) matfun( @( z ) L * z, x );
    nonlinearity.g = @( t,x ) g( x );
    nonlinearity.w = @( t,x ) w( x );
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
    dW_stash = stoch_init( k, Nt, Ntref, Ns, omega, test.MoCa ); % initialize stochasticity
    memo_stoch( { dW_stash{ length( Nt ) + 1,: } } );
    t = 0;
    uref = u0;
    fprintf(['- Compute reference with ', num2str( Ntref ),' timesteps of ', test.integrator, ' in combination with ', test.routines{ 1 },': %3.0f%%'], 100 * 1 / Ntref );
    for iter = 1 : Ntref
      fprintf('\b\b\b\b%3.0f%%', 100 * iter / Ntref );
      if strcmp( test.routines{ 1 },'bamphi' )
        [ uref, info ] = feval( [ test.integrator, '_', test.routines{ 1 } ], uref, k, t, linearity, nonlinearity, opts, info, aux_linrt );
      else
        [ uref, info ] = feval( [ test.integrator, '_', test.routines{ 1 } ], uref, k, t, linearity, nonlinearity, opts, info );
      end
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
      memo_stoch( { dW_stash{ l,: } } ); % initialize stochasticity
      t = 0;
      u = u0;
      tic;
      for iter = 1 : Nt( l )
        if strcmp( test.routines{ rout },'bamphi' )
          [ u, info ] = feval( [ test.integrator, '_', test.routines{ rout } ], u, k, t, linearity, nonlinearity, opts, info, aux_linrt );
        else
          [ u, info ] = feval( [ test.integrator, '_', test.routines{ rout } ], u, k, t, linearity, nonlinearity, opts, info );
        end
        t = t + k;
      end % marching
      [ ~, mv ] = matfun( NaN,NaN );
      memo.(test.integrator).(test.routines{ rout }).tstep( l ) = Nt( l );
      memo.(test.integrator).(test.routines{ rout }).clock( l ) = toc;
      memo.(test.integrator).(test.routines{ rout }).matve( l ) = mv;
      if test.compute_error
        memo.(test.integrator).(test.routines{ rout }).error( l ) = norm( mean( u, 2 ) - mean( uref, 2 ), inf ) / norm( mean( uref, 2 ), inf );
      end
      if strcmp( test.routines{ rout },'bamphi' )
        memo.(test.integrator).(test.routines{ rout }).info{ l } = info;
      end
      disp(['     ', num2str( l ), ' on ', num2str( length( Nt ) ),' - ', num2str( Nt(l) ),' timesteps... ', num2str( toc ),' seconds.' ]);
    end % timesteps
  end % routines

end

function dW_stash = stoch_init( dtref, Nt, Ntref, Ns, omega, MoCa )

  load('seed.mat'), rng( my_seed ),
  disp([ '      stoch_init: resetted seed: ', num2str( randn ) ]);
  tic

  bj = get_2D_bj( dtref, Ns( 1 ), Ns( 2 ), omega.x.l, omega.x.r, omega.y.l, omega.y.r, 0.5 ); % Get form of noise
  for i = 1 : Ntref,
    dW_stash{ length( Nt ) + 1, i } = get_2D_dW( bj, 1, MoCa ); % generate noise for each timestep, for each gridpoint and for each Monte Carlo run
  end

  for j = 1 : length( Nt )
    for i = 1 : Nt( j )
      k = Ntref / Nt(j);
      dW_stash{ j, i } = 0;
      for l = 1 : k
        dW_stash{ j, i } = dW_stash{ j, i } + dW_stash{ length( Nt ) + 1, l + ( i - 1 ) * k };
      end
    end
  end

  %
  disp(['      stoch_init: done initializing stochasticity in ', num2str(toc),' secs.']);

end

function bj = get_2D_bj( dtref, Nx, Ny, xl, xr, yl, yr, alpha )

  lambdax = 2 * pi * [ 0 : Nx / 2, - Nx / 2 + 1 : -1 ]' / ( xr - xl );
  lambday = 2 * pi * [ 0 : Ny / 2, - Ny / 2 + 1 : -1 ]' / ( yr - yl );
  % corrected TS Dec 2015
  [ lambdaxx, lambdayy ] = ndgrid( lambdax, lambday );
  root_qj = exp( - alpha * ( lambdaxx.^2 + lambdayy.^2 ) / 2 );   % smooth noise
  bj = root_qj * sqrt( dtref ) * Nx * Ny / sqrt( ( xr - xl ) * ( yr - yl ) );

end

function [ dW1, dW2 ] = get_2D_dW( bj, kappa, M )

  J = size( bj );
  if ( kappa == 1 )
    nnr = randn( J(1), J(2), M );
    nnc = randn( J(1), J(2), M );
  else
    nnr = squeeze( sum( randn( J(1), J(2), M, kappa ), 4 ) );
    nnc = squeeze( sum( randn( J(1), J(2), M, kappa ), 4 ) );
  end
  nn2 = nnr + 1i * nnc;
  TMPHAT = bsxfun( @times, bj, nn2 );
  tmp = ifft2( TMPHAT );
  dW1 = real( tmp );
  dW2 = imag( tmp );

end
