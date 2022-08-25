function [ x, varargout ] = matfun( L, x )
% Function merely needed to precisely count the number of matrix-vector products
% in a given run.
  persistent i
  if isempty( i )
    i = 0;
  end
  varargout{ 1 } = i;
  if not( isa( L, 'function_handle' ) )
    x = [];
    return
  end
  i = i + 1;
  x = L( x );
end
