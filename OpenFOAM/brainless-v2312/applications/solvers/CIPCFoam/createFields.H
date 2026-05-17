#include "createRDeltaT.H"

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

Info<< "Constructing twoPhaseMixtureThermo\n" << endl;
twoPhaseMixtureThermo mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

Info<< "Forcing registry of transportProperties for turbulence model\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",           // 메모리 창고(objectRegistry)에 등록될 이름
        runTime.constant(),              // 파일 위치 (constant 폴더)
        mesh,                            // 메쉬 참조
        IOobject::MUST_READ_IF_MODIFIED, // 파일이 있으면 무조건 읽어라
        IOobject::NO_WRITE               // 연산 중에 파일에 다시 쓰지는 않음
    )
);

Info<< "Reading thermophysical properties\n" << endl;

const volScalarField& rho1 = mixture.thermo1().rho();
const volScalarField& rho2 = mixture.thermo2().rho();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2
);


dimensionedScalar pMin
(
    "pMin",
    dimPressure,
    mixture
);

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


// Mass flux
// Initialisation does not matter because rhoPhi is reset after the
// alpha1 solution before it is used in the U equation.
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

volScalarField dgdt(alpha1*fvc::div(phi));

#include "createAlphaFluxes.H"

// Construct compressible turbulence model
compressibleInterPhaseTransportModel turbulence
(
    rho,
    U,
    phi,
    rhoPhi,
    alphaPhi10,
    mixture
);

#include "createK.H"

#include "createMRF.H"
#include "createFvOptions.H"

// Lee 모델 상수들
dimensionedScalar Tsat("Tsat", dimTemperature, 373.15); // 포화온도
dimensionedScalar L("L", dimEnergy/dimMass, 2257e3);    // 잠열
dimensionedScalar coeff("coeff", dimless/dimTime, 0.1); // 증발 계수 (작게 시작!)

// mDot(증발량)을 눈으로 확인하기 위한 필드 변수
volScalarField mDot
(
    IOobject("mDot", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("mDot", dimMass/dimVolume/dimTime, 0)
);