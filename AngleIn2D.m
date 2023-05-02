function Angle = AngleIn2D(V1, V2)
% Source:
% Jan on 17 Jul 2018
% via
% https://www.mathworks.com/matlabcentral/answers/410724-measure-angles-
% between-two-vectors-solely-counter-clockwise
% Input: V1, V2: [N x 3] vectors, each one can be [1 x 3] also.
%        The angle between the vectors V1 and V2 is calculated
% 
% Method: atan2(norm(N1 x N2), DOT(N1, N2))
% W. Kahan suggested in "Mindeless.pdf":
% 2 * atan(norm(x*norm(y) - norm(x)*y) / norm(x * norm(y) + norm(x) * y))
N1 = V1 ./ sqrt(sum(V2 .* V2, 2));   % >= R2016b: arithmetic expanding!
N2 = V2 ./ sqrt(sum(V2 .* V2, 2));   % >= R2016b: arithmetic expanding!
% Calculate dot and cross product of vectors:
N1dotN2 = N1(:, 1) .* N2(:, 1) + N1(:, 2) .* N2(:, 2);
N1xN2   = (N1(:, 1) .* N2(:, 2) - N1(:, 2) .* N2(:, 1));
% Angle between N1xN2 and view vector:
signLXo = sign(N1xN2);
% Care about anti-parallel N1 and N2:
antiN1N2 = (N1dotN2 < -0.999999999999993);       % -1 + 3 * EPS
if any(antiN1N2)
   antiN1N2 = and(antiN1N2, isfinite(N1dotN2));  % Catch N1dotN2=-Inf
   signLXo(antiN1N2) = 1.0;
end
normN1xN2 = abs(N1xN2);
Angle     = signLXo .* atan2(normN1xN2, N1dotN2);
end