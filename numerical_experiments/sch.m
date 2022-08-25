function memo = sch( test )
% NONSMOOTH CUBIC SCHRODINGER equation:
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

  % g  = @( u ) mu * ( u .* abs( u ).^2 );
  g  = @( u ) mu; % c'e un motivo ma non me lo ricordo (guardare integratore LR)

  % Space discretization
  Ns = test.Ns;
  b_c.x.l = 'hom_dirichlet'; b_c.x.r = 'hom_dirichlet';
  b_c.y.l = 'hom_dirichlet'; b_c.y.r = 'hom_dirichlet';
  [  D1,  D2, h, X, Y ] = space_discr( omega, Ns, b_c );
  [ DD1, DD2          ] = kronsum_struct( Ns, D1, D2, alpha, epsilon );
  L = DD2;
  Lfun = @( x ) L * x;

  % Initial datum
  u0 = generate_random_nonsmooth_u0( Ns, omega, test.theta ); u0 = u0(:);
  u0 = u0 .* ( ( X(:) - omega.x.l ) .* ( omega.x.r - X(:) ) .* ( Y(:) - omega.y.l ) .* ( omega.y.r - Y(:) ) ).^2;
  u0 = u0 / norm( u0, inf );

  % Time discretization
   Nt = test.Nt;
  tol = min( h.x^2, h.y^2 ) * 1.00e-003;
     linearity = @(   x ) matfun( @( z ) Lfun( z ), x );
  nonlinearity = @( t,x ) g( x );

  % Compute reference
  if test.compute_error
    Ntref = 3 * Nt( end );
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
      opts = []; info = [];
      if     strcmp( test.routines{ rout }, 'bamphi' )
        opts.tol = tol;
        clear cuschRS21_bamphi
        clear cuschss68_bamphi
      elseif strcmp( test.routines{ rout }, 'kiops' )
        opts = tol;
        clear cuschRS21_kiops
        clear cuschss68_kiops
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
        memo.(test.integrator).(test.routines{ rout }).error( l ) = norm( u - uref, inf ) / norm( uref, inf );
      end
      if strcmp( test.routines{ rout },'bamphi' )
        memo.(test.integrator).(test.routines{ rout }).info{ l } = info;
      end
      disp(['     ', num2str( l ), ' on ', num2str( length( Nt ) ),' - ', num2str( Nt(l) ),' timesteps... ', num2str( toc ),' seconds.' ]);
    end % timesteps
  end % routines

end

function u0 = generate_random_nonsmooth_u0( Ns, omega, theta )

  lambda_x = 2i * pi * ( - Ns( 1 ) / 2 : Ns( 1 ) / 2 - 1 )' / ( omega.x.r - omega.x.l );
  lambda_y = 2i * pi * ( - Ns( 2 ) / 2 : Ns( 2 ) / 2 - 1 )' / ( omega.y.r - omega.y.l );
  [ Lambda_x, Lambda_y ] = ndgrid( lambda_x, lambda_y );

  cof = ( sqrt( omega.x.r - omega.x.l ) / Ns( 1 ) ) * ( sqrt( omega.y.r - omega.y.l ) / Ns( 2 ) );
  my_fft  = @( u    ) fftshift(      fftn( u    ) ) * cof;
  my_ifft = @( uhat )    ifftn( ifftshift( uhat ) ) / cof;

  u0 = randn( Ns( 1 ), Ns( 2 ) ); % u0 = sin( pi * X ) .* sin( pi * Y );
  u0hat = my_fft( u0 );

  u0hat = u0hat .* abs( Lambda_x .* Lambda_y ) .^ ( - theta );
  u0hat( isinf( u0hat ) ) = 0;
  u0 = my_ifft( u0hat );

  u0_L2 = sqrt( sum( abs( u0(:) ).^2 ) * ( ( omega.x.r - omega.x.l ) / Ns( 1 ) * ( omega.y.r - omega.y.l ) / Ns( 2 ) ) );
  u0_L2 = u0_L2 + ( u0_L2 == 0 );
     u0 = u0 / u0_L2; % L2 normal

end
