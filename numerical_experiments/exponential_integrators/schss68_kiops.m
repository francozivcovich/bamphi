function [ u, m_opt ] = cuschss68_kiops( u, k, t, A, g, tol, m_opt )

  mu = g( NaN, NaN );
  if isempty( m_opt )
    m_opt = 10 * ones( 1,1 ); % set to default
  end

  u = exp( ( - 1i * k / 2 ) * ( u .* conj( u ) ) ) .* u;
  [ u, m_opt ] = kiops(  1i * k, A, u, tol, m_opt, [], [], false );
  u = exp( ( - 1i * k / 2 ) * ( u .* conj( u ) ) ) .* u;


end
