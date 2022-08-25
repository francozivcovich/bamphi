function [ u, info ] = exptg2s1_bamphi( u, k, t, A, g, opts, info, B )

  persistent store_info

  if length( store_info ) > 0
    info = store_info{ 1 };
  else
    if ( nargin == 8 )
      aux = bamphi_fov( length( u ), B ); % B is a matrix we rather use for computing A's FoV
      info = aux;
      info.A.pts.rho(:) = 1i * imag( info.A.pts.rho(:) );
    end
  end

  [ u, info ] = bamphi( k, A, [], u + ( k / 2 ) * g( t, u ), opts, info );
  store_info{ 1 } = info;

  u = u + ( k / 2 ) * g( t, u );


end
