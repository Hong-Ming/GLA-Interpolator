function Observation = MakeObservation(GroundTrue,sigma_n)

% Make observation
Noise = sigma_n*randn(size(GroundTrue));
Observation = GroundTrue + Noise;

end
