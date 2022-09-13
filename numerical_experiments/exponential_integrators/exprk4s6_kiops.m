function [ u, m_opt ] = exprk4s6_kiops( u, k, t, A, g, tol, m_opt )
% exprk4s6_kiops:
%
% There is only one strategy to run this exponential integrator in combination
% with kiops: the standard one.
% It is extremely advantageous to remember the suggested Arnoldi size m_opt,
% that is why we do output it.
%

  if isempty( m_opt )
    m_opt = 10 * ones( 1,4 ); % set to default
  end

  gn = g( t, u );
  Fn = A( u ) + gn;
  z = zeros( length( u ), 1 );

  U = zeros( length( u ), 2 );
  [ U( :,1 ), m_opt( 1 ), stats ] = kiops( 1/2 * k, A, [ z, Fn ], tol, m_opt( 1 ), [], [], false );
  U( :,1 ) = u + U( :, 1 );

  [ U, m_opt( 2 ), stats ] = kiops( k * [ 1/3, 1/2 ], A,[ z, Fn, 2 / k * ( g( t + 1/2 * k, U( :,1 ) ) - gn ) ], tol, m_opt( 2 ), [], [], false );
  U( :,2 ) = g( t + 1 / 2 * k, u + U( :,2 ) ) - gn;
  U( :,1 ) = g( t + 1 / 3 * k, u + U( :,1 ) ) - gn;

  [ U, m_opt( 3 ), stats ] = kiops( k * [ 1/3, 5/6 ], A, [ z, Fn, 6 / k * ( - 2/3 * U( :,2 ) + 3/2 * U( :,1 ) ), 12 / ( k^2 ) * ( 2 * U( :,2 ) - 3 * U( :,1 ) ) ], tol, m_opt( 3 ), [], [], false );

  U( :,2 ) = g( t + 5 / 6 * k, u + U( :,2 ) ) - gn;
  U( :,1 ) = g( t + 1 / 3 * k, u + U( :,1 ) ) - gn;

  [ z, m_opt( 4 ), stats ] = kiops( k, A, [ z, Fn, 2 / k * ( -2/5 * U( :,2 ) + 5/2 * U( :,1 ) ), 4 / ( k^2 ) * ( 6/5 * U( :,2 ) - 3 * U( :,1 ) ) ], tol, m_opt( 4 ), [], [], false );

  u = u + z;

end
