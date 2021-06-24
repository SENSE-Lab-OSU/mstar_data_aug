function ee = reconError(y1,gaussWidth,pulseSelectIDx,dist,azimuthVals,numFreqBins,numRangeBins,A_mod,groups_l12,sigma_n)
    BasisFunc = exp(-0.5*dist.^2/gaussWidth^2);
    normBasisFunc = diag((sum(BasisFunc.^2,1)).^0.5);
    BasisFunc = BasisFunc/normBasisFunc;
    A = @(x,mode) SAR_operator_gen(x,mode,pulseSelectIDx,numRangeBins,numFreqBins,azimuthVals,BasisFunc,A_mod(:,:,pulseSelectIDx));
    opts = spgSetParms('iscomplex',1,'verbosity',0);
    %Find Optimum C
    C = spg_group(A, y1, groups_l12, sigma_n, opts );
    ee = 0.5*norm(A(C,1)-y1)^2;
end