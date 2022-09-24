function [ u, m_opt ] = setdm1_kiops( u, k, t, A, nonlin, tol, m_opt )
% setdm1_kiops
%
% There is only one strategy to run this exponential integrator in combination
% with kiops: the standard one.
% It is extremely advantageous to remember the suggested Arnoldi size m_opt,
% that is why we do output it.
%
% The real trick is applied externally: one could loop over the L = 100 Monte
% Carlo run or rather form a huge matrix as described in the paper. The second
% approach is rewarding even though for growing L it might finally cause memory
% limitations to manifest.
%
  if isempty( m_opt )
    m_opt = 10 * ones( 1,1 ); % set to default
  end
  step = round( t / k ) + 1;
  dW = memo_stoch( step );

  [ du, m_opt ] = kiops( k, A, [ zeros( numel( u ),1 ),  A( u(:) + nonlin.w( t,u(:) ) .* dW(:) ) + nonlin.g( t, u(:) ) ], tol, m_opt, [], [], false );
  u = u + reshape( du + nonlin.w( t,u(:) ) .* dW(:), size( u ) );

end
