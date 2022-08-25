function [ DD1, DD2 ] = kronsum_struct( N, D1, D2, alpha, epsilon )
% Given the one-dimensional discretization of differential operators D1 and D2
% returns the required d-dimensional differential operators.
% The parameters alpha and epsilon stand, respectively, for the advection and
% diffusion coefficients.

  space_vars = { 'x', 'y', 'z' };
  DIM = length( N );
  if ( nargin < 4 ) || isempty( alpha )
    for i = 1 : DIM
      alpha.(space_vars{ i }) = 1;
    end
  end
  if ( nargin < 5 ) || isempty( epsilon )
    for i = 1 : DIM
      epsilon.(space_vars{ i }) = 1;
    end
  end

  DD1 = sparse( prod( N ), prod( N ) );
  DD2 = sparse( prod( N ), prod( N ) );
  for i = DIM : - 1 : 1
    %
    temp = 1;
    for j = DIM : -1 : i + 1
      temp = kron( temp, speye( N( j ) ) );
    end
    temp = kron( temp, alpha.(space_vars{ i }) * D1.(space_vars{ i }) );
    for j = i - 1 : -1 : 1
      temp = kron( temp, speye( N( j ) ) );
    end
    DD1 = DD1 + temp;
    %
    temp = 1;
    for j = DIM : -1 : i + 1
      temp = kron( temp, speye( N( j ) ) );
    end
    temp = kron( temp, epsilon.(space_vars{ i }) * D2.(space_vars{ i }) );
    for j = i - 1 : -1 : 1
      temp = kron( temp, speye( N( j ) ) );
    end
    DD2 = DD2 + temp;
  end

end
