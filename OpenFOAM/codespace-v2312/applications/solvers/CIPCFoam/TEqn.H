{
    // 잠열 (Latent Heat, 예: 물의 경우 약 2.26e6 J/kg)
    // 실제로는 mixture나 thermo에서 적절한 값을 가져와야 합니다.
    dimensionedScalar L("L", dimEnergy/dimMass, 2.26e6); 

    fvScalarMatrix TEqn
    (
        fvm::ddt(rho, T) + fvm::div(rhoPhi, T) - fvm::Sp(contErr, T)
      - fvm::laplacian(turbulence.alphaEff(), T)
      + (
             divUp()
          +  (fvc::ddt(rho, K) + fvc::div(rhoPhi, K))() - contErr*K
        )
       *(
           alpha1()/mixture.thermo1().Cv()()
         + alpha2()/mixture.thermo2().Cv()()
        )
     ==
        fvOptions(rho, T)
      - (mDot * L) / mixture.Cv() // <--- 잠열 소스항 추가!
    );

    TEqn.relax();
    fvOptions.constrain(TEqn);
    TEqn.solve();
    fvOptions.correct(T);

    mixture.correctThermo();
    mixture.correct();
}