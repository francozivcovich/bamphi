function [ u, m_opt ] = cuschRS21_kiops( u, k, t, A, g, tol, m_opt )

  mu = g( NaN, NaN );
  if isempty( m_opt )
    m_opt = 10 * ones( 1,5 ); % set to default
  end
  ut = u;
  [ u,   m_opt( 1 ) ] = kiops(  1i * k, A, conj( ut ), tol, m_opt( 1 ), [], [], false );
  [ u,   m_opt( 2 ) ] = kiops( -2i * k, A, [ zeros( size( ut ) ), u / ( -2i * k ), u / ( 4 * k^2 ) ], tol, m_opt( 2 ),[],[], false );
  [ ut2, m_opt( 3 ) ] = kiops(  1i * k, A, ut, tol, m_opt( 3 ), [], [], false );
  ut2 = ( -1i * mu * k ) * ( ut2.^2 .* u );
  [ u,   m_opt( 4 ) ] = kiops( -2i * k, A, [ zeros( length( ut ), 2 ), -conj( ut ) / ( 4 * k^2 ) ], tol, m_opt( 4 ), [], [], false );
  u = ut - ( 1i * mu * k ) * ut.^2 .* u - ( mu^2 * k^2 * (1/2) ) * ( ut .* conj( ut ) ).^2 .* ut;
  [ u,   m_opt( 5 ) ] = kiops(  1i * k, A, u, tol, m_opt( 5 ), [], [], false );
  u = u + ut2;

end
