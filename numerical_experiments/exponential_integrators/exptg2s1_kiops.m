function [ u, m_opt ] = exptg2s1_kiops( u, k, t, A, g, tol, m_opt )

  if isempty( m_opt )
    m_opt = 10 * ones( 1,1 ); % set to default
  end

  [ u, m_opt ] = kiops( k, A, u + ( k / 2 ) * g( t, u ), tol, m_opt, [], [], false );

  u = u + ( k / 2 ) * g( t, u );


end
