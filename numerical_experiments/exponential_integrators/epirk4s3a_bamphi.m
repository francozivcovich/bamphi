function [ u, info ] = epirk4s3a_bamphi( u, k, t, Jfun, f, opts, info )
% epirk4s3a_bamphi
%
% There are two strategies to carry on this exponential integration with bamphi:
%
% 1. (default) compute (IOM)-Arnoldi at each call of bamphi (two calls per step)
% 2. compute (IOM)-Arnoldi once per step at the first call of bamphi
%
% If the user is really adverse to (IOM)-Arnoldi then 2. is motivated, otherwise
% the strategy 1. might be a much better alternative. Only the 1st one is left
% uncommented as it is the suggested one for this kind of task
%
% Another strategy is to recycle info as we do in strategy 2. of exprk4s6_bamphi.
% Surprisingly such a strategy works even better for this example but it's too
% unjustified (J changes at each timestep) to be used in a piece of scientific
% research.
% Practitioners, on the other hand, often refresh the Jacobian matrix ever once
% in a while, in that case it would be justified to use strategy 2. as in
% exprk4s6_bamphi.
%


  % 1st STRATEGY:
  persistent suggested_arnoldi_size

  fn = f( u );
  rfun = @( v ) f( v ) - fn - Jfun( v - u );

  if ( length( suggested_arnoldi_size ) > 0 )
    opts.r_arn = suggested_arnoldi_size( 1 );
  end
  [ U, info ] = bamphi( k * [ 1/2, 2/3 ], Jfun, [], [zeros(length(fn),1),fn], opts, info );
  suggested_arnoldi_size( 1 ) = bamphi_suggest_arnoldi_size( opts, info );
  info = rmfield( info, 'A' );

  r2 = rfun( u + U( :,1 ) );
  r3 = rfun( u + U( :,2 ) );

  if ( length( suggested_arnoldi_size ) > 1 )
    opts.r_arn = suggested_arnoldi_size( 2 );
  end
  [ fn, info ] = bamphi( k, Jfun, [], [zeros(length(fn),1),fn,zeros(length(fn),1),1/k^2*(32*r2-27/2*r3),1/k^3*(-144*r2+81*r3)], opts, info );
  suggested_arnoldi_size( 2 ) = bamphi_suggest_arnoldi_size( opts, info );
  info = rmfield( info, 'A' );

  % % 2nd STRATEGY:
  % persistent suggested_arnoldi_size
  %
  % fn = f( u );
  % rfun = @( v ) f( v ) - fn - Jfun( v - u );
  %
  % if ( length( suggested_arnoldi_size ) > 0 )
  %   opts.r_arn = suggested_arnoldi_size( 1 );
  % end
  % [ U, info ] = bamphi( k * [ 1/2, 2/3 ], Jfun, [], [zeros(length(fn),1),fn], opts, info );
  % suggested_arnoldi_size( 1 ) = bamphi_suggest_arnoldi_size( opts, info );
  %
  % r2 = rfun( u + U( :,1 ) );
  % r3 = rfun( u + U( :,2 ) );
  %
  % [ fn, info ] = bamphi( k, Jfun, [], [zeros(length(fn),1),fn,zeros(length(fn),1),1/k^2*(32*r2-27/2*r3),1/k^3*(-144*r2+81*r3)], opts, info );
  % suggested_arnoldi_size( 2 ) = bamphi_suggest_arnoldi_size( opts, info );
  % info = rmfield( info, 'A' );
  %
  % u = u + fn;

end
