function [ u, m_opt ] = exptg2s1_kiops( u, k, t, A, g, tol, m_opt )
% exptg2s1_kiops
%
% There is only one strategy to run this exponential integrator in combination
% with kiops: the standard one.
% It is extremely advantageous to remember the suggested Arnoldi size m_opt,
% that is why we do output it.
%
% Instead of A we could have give Lhat (see paper), on the other hand the much
% worsen sparsity pattern makes this choice not clever (try to change test.mode
% in main.m inside the "if strcmp( which_tests{ test_idx }, 'skg' )"" statement).
%

  if isempty( m_opt )
    m_opt = 10 * ones( 1,1 ); % set to default
  end

  [ u, m_opt ] = kiops( k, A, u + ( k / 2 ) * g( t, u ), tol, m_opt, [], [], false );

  u = u + ( k / 2 ) * g( t, u );


end
