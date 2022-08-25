function [ u, info ] = cuschRS21_bamphi( u, k, t, A, g, opts, info )

  persistent store_info

  mu = g( NaN, NaN );
  ut = u;

  if length( store_info ) > 0,
    info = store_info{ 1 };
  end
  [ U, info ] = bamphi( [ -1i * k, 1i * k ], A, [], ut, opts, info );
  store_info{ 1 } = info;

  if length( store_info ) > 1,
    info = store_info{ 2 };
  end
  [ u, info ] = bamphi( -2i * k, A, [], [ zeros( size( ut ) ), conj( U(:,1) ) / ( -2i * k ), conj( U( :,1 ) ) / ( 4 * k^2 ) ], opts, info );
  store_info{ 2 } = info;

  U( :,2 ) = ( -1i * mu * k ) * ( U( :,2 ).^2 .* u );

  if length( store_info ) > 2,
    info = store_info{ 3 };
  end
  [ u, info ] = bamphi( -2i * k, A, [], [ zeros( length( ut ), 2 ), -conj( ut ) / ( 4 * k^2 ) ], opts, info );
  store_info{ 3 } = info;

  u = ut - ( 1i * mu * k ) * ut.^2 .* u - ( mu^2 * k^2 * (1/2) ) * ( ut .* conj( ut ) ).^2 .* ut;

  if length( store_info ) > 3,
    info = store_info{ 4 };
  end
  [ u, info ] = bamphi( 1i * k, A, [], u, opts, info );
  store_info{ 4 } = info;
  
  u = u + U( :,2 );

end
