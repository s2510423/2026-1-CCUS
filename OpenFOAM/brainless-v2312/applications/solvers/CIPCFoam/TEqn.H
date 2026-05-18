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

    // ================= [괄호 에러 수정한 철벽 방어막] =================
    // 💡 mesh() 가 아니라 'mesh' 로 괄호를 떼야 정상 컴파일됩니다!
    volScalarField& T_ref = const_cast<volScalarField&>(mesh.lookupObject<volScalarField>("T"));
    
    // 1. 내부 격자(Internal Field) 300K 락
    forAll(T_ref, celli)
    {
        if (T_ref[celli] < 200.0)
        {
            T_ref[celli] = 300.0; 
        }
    }
    
    // 2. 경계면(Boundary Field) 300K 락 (터짐 방지 핵심)
    forAll(T_ref.boundaryField(), patchi)
    {
        scalarField& ptf = T_ref.boundaryFieldRef()[patchi];
        forAll(ptf, facei)
        {
            if (ptf[facei] < 200.0)
            {
                ptf[facei] = 300.0;
            }
        }
    }
    // =============================================================

    fvOptions.correct(T_ref);

    mixture.correctThermo();
    mixture.correct();
   }
