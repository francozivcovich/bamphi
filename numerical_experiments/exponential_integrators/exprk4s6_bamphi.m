function [ u, info ] = exprk4s6_bamphi( u, k, t, A, g, opts, info, B )

  persistent store_info

  gn = g( t, u );
  Fn = A( u ) + gn;
  z = zeros( length( u ), 1 );
  U = zeros( length( u ), 2 );

  if length( store_info ) > 0
    info = store_info{ 1 };
  else
    if ( nargin == 8 )
      aux = bamphi_fov( length( u ), B ); % B is a matrix we rather use for computing A's FoV
      info = aux;
    end
  end
  [ U(:,1), info ] = bamphi( 1/2 * k, A, [], [ z, Fn ], opts, info );
  store_info{ 1 } = info;
  U(:,1) = u + U(:,1);

  if length( store_info ) > 1
    info = store_info{ 2 };
  else
    if ( nargin == 8 )
      info = aux;
    end
  end
  [ U, info ] = bamphi( k * [ 1/3, 1/2 ], A, [], [z,Fn,2/k*(g(t+1/2*k,U(:,1)) - gn)], opts, info );
  store_info{ 2 } = info;
  info = [];
  U(:,2) = g( t + 1/2 * k, u + U(:,2) ) - gn;
  U(:,1) = g( t + 1/3 * k, u + U(:,1) ) - gn;

  if length( store_info ) > 2
    info = store_info{ 3 };
  else
    if ( nargin == 8 )
      info = aux;
    end
  end
  [ U,info ] = bamphi( k * [ 1/3, 5/6 ], A, [], [z,Fn,6/k*(3/2*U(:,1)-2/3*U(:,2)),12/(k^2)*(2*U(:,2)-3*U(:,1))], opts, info );
  store_info{ 3 } = info;
  info = [];

  U(:,2) = g( t + 5/6 * k, u + U(:,2) ) - gn;
  U(:,1) = g( t + 1/3 * k, u + U(:,1) ) - gn;

  if length( store_info ) > 3
    info = store_info{ 4 };
  else
    if ( nargin == 8 )
      info = aux;
    end
  end
  [ z,info ] = bamphi( k, A, [], [z,Fn,2/k*(5/2*U(:,1)-2/5*U(:,2)),4/(k^2)*(6/5*U(:,2)-3*U(:,1))], opts, info );
  store_info{ 4 } = info;
  info = [];
  u = u + z;


end
