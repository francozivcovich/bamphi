function [ u, info ] = stochexprk1_bamphi( u, k, t, A, nonlin, opts, info, B )

  persistent store_info

  if length( store_info ) > 0
    info = store_info{ 1 };
  else
    if ( nargin == 8 )
      aux = bamphi_fov( size( u,1 ), B ); % B is a matrix we rather use for computing A's FoV
      info = aux;
    end
  end


  step = round( t / k ) + 1;
  dW = memo_stoch( step );

  [ du, info ]  = bamphi( k, A, [], [ zeros( numel( u ),1 ),  A( u(:) + nonlin.w( t,u(:) ) .* dW(:) ) + nonlin.g( t, u(:) ) ], opts, info );
  store_info{ 1 } = info;
  u = u + reshape( du, size( u ) );

end
