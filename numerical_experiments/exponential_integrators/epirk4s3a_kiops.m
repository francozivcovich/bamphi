function [ u, m_opt ] = epirk4s3a_kiops( u, k, t, Jfun, f, tol, m_opt );
% epirk4s3a_kiops
%
% There is only one strategy to run this exponential integrator in combination
% with kiops: the standard one.
% It is extremely advantageous to remember the suggested Arnoldi size m_opt,
% that is why we do output it.
%
  if isempty( m_opt )
    m_opt = 10 * ones( 1,2 ); % set to default
  end

  fn = f( u );
  rfun = @(v) f( v ) - fn - Jfun( v - u );

  [ U, m_opt( 1 ) ] = kiops( k * [ 1 / 2, 2 / 3 ], Jfun, [ zeros(length(fn),1), k * fn ], tol, m_opt( 1 ) ); %[Un2,Un3]

  r2 = rfun( u + 1 / 2 * U( :, 1 ) );
  r3 = rfun( u + 2 / 3 * U( :, 2 ) );

  [ fn, m_opt( 2 ) ] =  kiops( k, Jfun, [ zeros(length(fn),1), fn, zeros(length(fn),1), 1/k^2*(32*r2-27/2*r3),1/k^3*(-144*r2+81*r3) ], tol, m_opt( 2 ),[],[],false);
  u = u + fn;


end
