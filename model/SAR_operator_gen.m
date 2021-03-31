function y =SAR_operator_gen(x,mode,pulseSelectIDx,numRangeBins,numFreqBins,azimuthVals,BasisFunc,Amod)

velLight = 3e8;

az =azimuthVals;
numTotalPulse = length(azimuthVals);
pulseSel = zeros(numTotalPulse,1);
pulseSel(pulseSelectIDx) = 1;

pulseSelRep = repmat(pulseSel.',numRangeBins^2,1);
pulseSelRep = pulseSelRep(:);
numPulsesSel = nnz(pulseSel);
numBasisFunctions = size(BasisFunc,2);

pulseCounts=1;
if mode == 1
    xx = reshape(x,numRangeBins^2,numBasisFunctions);
    xx_interp = (BasisFunc*xx.').';
    xx_interp = reshape(xx_interp,numRangeBins^2,1,numTotalPulse);
    xx_interp1 = xx_interp(:,:,pulseSel>0); 

    y=mtimesx(Amod, xx_interp1,'MATLAB');
    y=y(:);
%     Amod=permute(Amod,[3,1,2]);
%     y=double(py.numpy.matmul(Amod,py.numpy.transpose(xx_interp1,[2,0,1])));
%     y=y.';
%     y=y(:);    
else
    xx = reshape(x,numFreqBins,1,numPulsesSel);
    yy = zeros(numRangeBins^2*numTotalPulse,1);
    yTemp = mtimesx(Amod,'C', xx,'MATLAB');

    yy(pulseSelRep > 0) = yTemp(:);
    yy_mod = reshape(yy,numRangeBins^2,numTotalPulse);
    y = (BasisFunc'*yy_mod.').';
    y = y(:);
end
end