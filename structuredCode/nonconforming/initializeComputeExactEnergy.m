function initializeComputeExactEnergy(fct, parAlpha, parBeta)

  if nargin < 3, parBeta = NaN; end
  geometry = 'BigSquare';
  minNrInnerEdges = 1e9;
  minPrecision = 8;
  degree4Integrate = 20;

  switch fct
    case 'f01'
      if isnan(parBeta), error('missing parameter parBeta'); end
      fStr = 'f01';
      fStrParams = [parAlpha, parBeta];
      uStr = 'f01ExactSolution';
      uStrParams = parBeta;
      gradUStr = 'f01ExactSolutionGradient';
      gradUStrParams = parBeta;
    case 'f02'
      if isnan(parBeta), error('missing parameter parBeta'); end
      fStr = 'f02';
      fStrParams = [parAlpha, parBeta];
      uStr = 'f02ExactSolution';
      uStrParams = parBeta;
      gradUStr = 'f02ExactSolutionGradient';
      gradUStrParams = parBeta;
    case 'f03'
      fStr = 'f03';
      fStrParams = parAlpha;
      uStr = 'f03ExactSolution';
      uStrParams = [];
      gradUStr = 'f03ExactSolutionGradient';
      gradUStrParams = [];
    case 'f04'
      fStr = 'f04';
      fStrParams = parAlpha;
      uStr = 'f04ExactSolution';
      uStrParams = [];
      gradUStr = 'f04ExactSolutionGradient';
      gradUStrParams = [];
    otherwise
      error('unknown function');
  end

  computeExactEnergyBV(geometry, fStr, fStrParams, uStr, uStrParams, ...
    gradUStr, gradUStrParams, parAlpha, ...,
    minNrInnerEdges, minPrecision, degree4Integrate)
end
