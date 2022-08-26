function [ u, info ] = epirk4s3a_bamphi( u, k, t, Jfun, f, opts, info )

  persistent suggested_arnoldi_size

  fn = f( u );
  rfun = @( v ) f( v ) - fn - Jfun( v - u );

  if ( length( suggested_arnoldi_size ) > 0 )
    opts.r_arn = suggested_arnoldi_size( 1 );
  end
  [ U, info ] = bamphi( k * [ 1/2, 2/3 ], Jfun, [], [zeros(length(fn),1),fn], opts, info );
  suggested_arnoldi_size( 1 ) = bamphi_suggest_arnoldi_size( opts, info );
  info = rmfield( info, 'A' );

  r2 = rfun( u + U( :,1 ) );
  r3 = rfun( u + U( :,2 ) );

  if ( length( suggested_arnoldi_size ) > 1 )
    opts.r_arn = suggested_arnoldi_size( 2 );
  end
  [ fn, info ] = bamphi( k, Jfun, [], [zeros(length(fn),1),fn,zeros(length(fn),1),1/k^2*(32*r2-27/2*r3),1/k^3*(-144*r2+81*r3)], opts, info );
  suggested_arnoldi_size( 2 ) = bamphi_suggest_arnoldi_size( opts, info );
  info = rmfield( info, 'A' );

  u = u + fn;

end
