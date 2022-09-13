function [ u, m_opt ] = stochexprk1_kiops( u, k, t, A, nonlin, tol, m_opt )

  if isempty( m_opt )
    m_opt = 10 * ones( 1,1 ); % set to default
  end
  step = round( t / k ) + 1;
  dW = memo_stoch( step );

  [ du, m_opt ] = kiops( k, A, [ zeros( numel( u ),1 ),  A( u(:) + nonlin.w( t,u(:) ) .* dW(:) ) + nonlin.g( t, u(:) ) ], tol, m_opt, [], [], false );
  u = u + reshape( du, size( u ) );

end
