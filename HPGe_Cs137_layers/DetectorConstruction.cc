#include "DetectorConstruction.hh"
#include "DetectorParameters.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Polycone.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Exception.hh"

#include <cmath>

namespace
{
    using namespace DetectorParameters;

    constexpr G4double kGeomTolerance = 1.0e-6 * mm;

    constexpr G4double kJoinOverlap = 0.01 * mm;

    G4Material* MakeDeadLayerMaterial()
    {
        auto* existing = G4Material::GetMaterial("DeadLayerGe", false);
        if (existing != nullptr)
        {
            return existing;
        }

        auto* nist = G4NistManager::Instance();
        auto* material = new G4Material(
            "DeadLayerGe",
            5.23 * g / cm3,
            1
        );
        material->AddElement(nist->FindOrBuildElement("Ge"), 1);
        return material;
    }

    // The crystal
    G4VSolid* MakeCrystalSolid(
        const G4String& name,
        G4double diameter,
        G4double height,
        G4double topDiameter,
        G4double chamferHeight,
        G4double holeDiameter,
        G4double holeDepth,
        G4double grooveInnerDiameter,
        G4double grooveOuterDiameter,
        G4double grooveDepth)
    {
        const G4double zBottom = -height / 2.0;
        const G4double zChamfer = height / 2.0 - chamferHeight;
        const G4double zTop = height / 2.0;

        G4double zPlanes[3] = {zBottom, zChamfer, zTop};
        G4double rInner[3] = {0.0, 0.0, 0.0};
        G4double rOuter[3] = {
            diameter / 2.0,
            diameter / 2.0,
            topDiameter / 2.0
        };

        G4VSolid* solid = new G4Polycone(
            name + "_Outer",
            0.0,
            twopi,
            3,
            zPlanes,
            rInner,
            rOuter
        );

        auto* bore = new G4Tubs(
            name + "_Bore",
            0.0,
            holeDiameter / 2.0,
            holeDepth / 2.0 + kGeomTolerance,
            0.0,
            twopi
        );
        const G4double boreCenterZ = zBottom + holeDepth / 2.0;

        solid = new G4SubtractionSolid(
            name + "_WithBore",
            solid,
            bore,
            nullptr,
            G4ThreeVector(0.0, 0.0, boreCenterZ)
        );

        auto* groove = new G4Tubs(
            name + "_Groove",
            grooveInnerDiameter / 2.0,
            grooveOuterDiameter / 2.0,
            grooveDepth / 2.0 + kGeomTolerance,
            0.0,
            twopi
        );
        const G4double grooveCenterZ = zBottom + grooveDepth / 2.0;

        solid = new G4SubtractionSolid(
            name + "_WithGroove",
            solid,
            groove,
            nullptr,
            G4ThreeVector(0.0, 0.0, grooveCenterZ)
        );

        return solid;
    }

    G4VSolid* MakeCylindricalShell(
        const G4String& name,
        G4double innerRadius,
        G4double outerRadius,
        G4double height)
    {
        return new G4Tubs(
            name,
            innerRadius,
            outerRadius,
            height / 2.0,
            0.0,
            twopi
        );
    }
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    auto* nist = G4NistManager::Instance();

    auto* worldMaterial = nist->FindOrBuildMaterial("G4_AIR");
    auto* germanium = nist->FindOrBuildMaterial("G4_Ge");
    auto* aluminium = nist->FindOrBuildMaterial("G4_Al");
    auto* carbon = nist->FindOrBuildMaterial("G4_C");
    auto* ptfe = nist->FindOrBuildMaterial("G4_TEFLON");
    auto* deadGermanium = MakeDeadLayerMaterial();

    auto* worldSolid = new G4Box(
        "WorldSolid",
        worldHalfSize,
        worldHalfSize,
        worldHalfSize
    );
    auto* worldLogical = new G4LogicalVolume(
        worldSolid,
        worldMaterial,
        "WorldLogical"
    );
    auto* worldPhysical = new G4PVPlacement(
        nullptr,
        G4ThreeVector(),
        worldLogical,
        "World",
        nullptr,
        false,
        0,
        true
    );
    worldLogical->SetVisAttributes(G4VisAttributes::GetInvisible());

    auto* fullCrystalSolid = MakeCrystalSolid(
        "FullCrystal",
        D,
        H,
        Q,
        T,
        P,
        centralHoleDepth,
        S1,
        S2,
        M
    );

    const G4double activeD = D - 2.0 * d2;
    const G4double activeH = H - d1 - d3;
    const G4double activeQ = Q - 2.0 * d2;
    const G4double activeT = T;
    const G4double activeP = P + 2.0 * d4;
    const G4double activeL = centralHoleDepth + d5 - d3;
    const G4double activeS1 = S1 - 2.0 * d6;
    const G4double activeS2 = S2 + 2.0 * d7;
    const G4double activeM = M + d8 - d3;
    const G4double activeShiftZ = (d3 - d1) / 2.0;

    if (activeD <= 0.0 || activeH <= 0.0 || activeQ <= 0.0 ||
        activeT <= 0.0 || activeP <= 0.0 || activeL <= 0.0 ||
        activeS1 <= activeP || activeS2 <= activeS1 || activeM <= 0.0)
    {
        G4Exception(
            "DetectorConstruction::Construct",
            "InvalidActiveCrystal",
            FatalException,
            "The calculated active crystal dimensions are invalid."
        );
    }

    auto* activeCrystalSolid = MakeCrystalSolid(
        "ActiveCrystal",
        activeD,
        activeH,
        activeQ,
        activeT,
        activeP,
        activeL,
        activeS1,
        activeS2,
        activeM
    );
    auto* activeCrystalLogical = new G4LogicalVolume(
        activeCrystalSolid,
        germanium,
        "ActiveCrystalLogical"
    );
    new G4PVPlacement(
        nullptr,
        G4ThreeVector(0.0, 0.0, activeShiftZ),
        activeCrystalLogical,
        "ActiveCrystal",
        worldLogical,
        false,
        0,
        true
    );

    auto* deadCrystalSolid = new G4SubtractionSolid(
        "DeadCrystalSolid",
        fullCrystalSolid,
        activeCrystalSolid,
        nullptr,
        G4ThreeVector(0.0, 0.0, activeShiftZ)
    );
    auto* deadCrystalLogical = new G4LogicalVolume(
        deadCrystalSolid,
        deadGermanium,
        "DeadCrystalLogical"
    );
    new G4PVPlacement(
        nullptr,
        G4ThreeVector(),
        deadCrystalLogical,
        "DeadCrystal",
        worldLogical,
        false,
        0,
        true
    );

    const G4double crystalTopZ = H / 2.0;
    const G4double crystalBottomZ = -H / 2.0;
    
    constexpr G4double capWindowOpeningDiameter = 69.7 * mm;

    constexpr G4double capWindowAluminiumThickness = 1.0 * mm;
    constexpr G4double capWindowRecessDepth = capWindowAluminiumThickness;
    constexpr G4double capFlangeHeight = 1.5 * mm;

    if (capMaxDiameter <= capStepDiameter ||
        capStepDiameter <= capInnerDiameter ||
        capInnerDiameter <= capWindowOpeningDiameter ||
        capWindowRecessDepth >= capFlangeHeight ||
        housingSideThickness <= 0.0)
    {
        G4Exception(
            "DetectorConstruction::Construct",
            "InvalidCapDimensions",
            FatalException,
            "The stepped cap dimensions are inconsistent."
        );
    }

    
    const G4double capTopClearance = windowGapAboveCrystal;
    const G4double capFlangeBottomTargetZ =
        crystalTopZ + ptfeSleeveThickness + kJoinOverlap
        + capTopClearance;

    const G4double windowRadius = capInnerDiameter / 2.0;
    
    const G4double windowBottomZ =
        capFlangeBottomTargetZ + capFlangeHeight
        - capWindowRecessDepth - thinWindowThickness;
    const G4double windowTopZ =
        windowBottomZ + thinWindowThickness;
    const G4double windowCenterZ =
        (windowBottomZ + windowTopZ) / 2.0;

    auto* windowSolid = new G4Tubs(
        "CarbonWindowSolid",
        0.0,
        windowRadius,
        thinWindowThickness / 2.0,
        0.0,
        twopi
    );
    auto* windowLogical = new G4LogicalVolume(
        windowSolid,
        carbon,
        "CarbonWindowLogical"
    );
    new G4PVPlacement(
        nullptr,
        G4ThreeVector(0.0, 0.0, windowCenterZ),
        windowLogical,
        "CarbonWindow",
        worldLogical,
        false,
        0,
        true
    );

    auto* carbonWindowAluminiumCoverSolid = new G4Tubs(
        "CarbonWindowAluminiumCoverSolid",
        0.0,
        capWindowOpeningDiameter / 2.0,
        capWindowAluminiumThickness / 2.0,
        0.0,
        twopi
    );
    auto* carbonWindowAluminiumCoverLogical = new G4LogicalVolume(
        carbonWindowAluminiumCoverSolid,
        aluminium,
        "CarbonWindowAluminiumCoverLogical"
    );
    new G4PVPlacement(
        nullptr,
        G4ThreeVector(
            0.0,
            0.0,
            windowTopZ + capWindowAluminiumThickness / 2.0
        ),
        carbonWindowAluminiumCoverLogical,
        "CarbonWindowAluminiumCover",
        worldLogical,
        false,
        0,
        true
    );

    const G4double capOuterRadius = capMaxDiameter / 2.0;

    constexpr G4double capSideClearance = 5.0 * mm;
    constexpr G4double capWallThickness = 1.0 * mm;
    const G4double capSideInnerRadius =
        plateDiameter / 2.0 + capSideClearance;
    const G4double capSideOuterRadius =
        capSideInnerRadius + capWallThickness;
    const G4double capClearanceToPlate =
        capSideInnerRadius - plateDiameter / 2.0;

    if (capWallThickness <= 0.0 ||
        capClearanceToPlate < capSideClearance - kGeomTolerance ||
        capSideOuterRadius >= capOuterRadius)
    {
        G4Exception(
            "DetectorConstruction::Construct",
            "InvalidCapClearance",
            FatalException,
            "The cap must preserve the requested 5 mm radial clearance and fit within the top flange."
        );
    }

   
    const G4double capTopZ = windowTopZ + capWindowRecessDepth;
    const G4double capFlangeBottomZ = capTopZ - capFlangeHeight;
    if (std::abs(capFlangeBottomZ - capFlangeBottomTargetZ) > kGeomTolerance)
    {
        G4Exception(
            "DetectorConstruction::Construct",
            "InvalidCapTopClearance",
            FatalException,
            "The cap flange is not at the requested top-clearance position."
        );
    }
    const G4double shoulderHeight =
        capFlangeHeight - capWindowRecessDepth;
    const G4double shoulderCenterZ =
        (capFlangeBottomZ + windowTopZ) / 2.0;
    const G4double lipHeight = capWindowRecessDepth;
    const G4double lipCenterZ =
        (windowTopZ + capTopZ) / 2.0;

    const G4double capCoverBottomZ =
        crystalTopZ - plateTopOffsetFromCrystalTop - plateThickness;


    const G4double capWallBottomZ = capCoverBottomZ;
    const G4double capWallTopZ = capFlangeBottomZ + kJoinOverlap;
    const G4double capWallHeight = capWallTopZ - capWallBottomZ;
    const G4double capWallCenterZ =
        (capWallBottomZ + capWallTopZ) / 2.0;

    auto* capWallSolid = MakeCylindricalShell(
        "CapWallSolid",
        capSideInnerRadius,
        capSideOuterRadius,
        capWallHeight
    );


    auto* capShoulderSolid = MakeCylindricalShell(
        "CapShoulderSolid",
        capInnerDiameter / 2.0,
        capOuterRadius,
        shoulderHeight + kJoinOverlap
    );
    const G4double capShoulderCenterZ =
        capFlangeBottomZ + (shoulderHeight + kJoinOverlap) / 2.0;

    auto* capUpperLipSolid = MakeCylindricalShell(
        "CapUpperLipSolid",
        capWindowOpeningDiameter / 2.0,
        capOuterRadius,
        lipHeight
    );

    auto* capSteppedRimSolid = new G4UnionSolid(
        "CapSteppedRimSolid",
        capShoulderSolid,
        capUpperLipSolid,
        nullptr,
        G4ThreeVector(
            0.0,
            0.0,
            lipCenterZ - capShoulderCenterZ
        )
    );

    auto* capAssemblySolid = new G4UnionSolid(
        "CapAssemblySolid",
        capWallSolid,
        capSteppedRimSolid,
        nullptr,
        G4ThreeVector(
            0.0,
            0.0,
            capShoulderCenterZ - capWallCenterZ
        )
    );

    auto* capAssemblyLogical = new G4LogicalVolume(
        capAssemblySolid,
        aluminium,
        "CapAssemblyLogical"
    );
    new G4PVPlacement(
        nullptr,
        G4ThreeVector(0.0, 0.0, capWallCenterZ),
        capAssemblyLogical,
        "CapAssembly",
        worldLogical,
        false,
        0,
        true
    );

    const G4double ptfeInnerRadius = D / 2.0;
    const G4double ptfeOuterRadius = ptfeInnerRadius + ptfeSleeveThickness;

    auto* ptfeSleeveSolid = MakeCylindricalShell(
        "PTFESleeveSolid",
        ptfeInnerRadius,
        ptfeOuterRadius,
        H
    );
    auto* ptfeSleeveLogical = new G4LogicalVolume(
        ptfeSleeveSolid,
        ptfe,
        "PTFESleeveLogical"
    );
    new G4PVPlacement(
        nullptr,
        G4ThreeVector(),
        ptfeSleeveLogical,
        "PTFESleeve",
        worldLogical,
        false,
        0,
        true
    );


    const G4double plateRadius = plateDiameter / 2.0;
    const G4double plateTopZ = crystalTopZ - plateTopOffsetFromCrystalTop;
    const G4double plateBottomZ = plateTopZ - plateThickness;
    const G4double plateCenterZ = (plateTopZ + plateBottomZ) / 2.0;

    auto* plateSolid = new G4Tubs(
        "MountingPlateSolid",
        0.0,
        plateRadius,
        plateThickness / 2.0,
        0.0,
        twopi
    );

    const G4double coldFingerBlockRadius = coldFingerBlockDiameter / 2.0;
    const G4double blockTopZ = plateBottomZ;
    const G4double blockBottomZ = blockTopZ - coldFingerBlockHeight;
    const G4double blockCenterZ = (blockTopZ + blockBottomZ) / 2.0;

    auto* coldFingerBlockSolid = new G4Tubs(
        "ColdFingerBlockSolid",
        0.0,
        coldFingerBlockRadius,
        coldFingerBlockHeight / 2.0,
        0.0,
        twopi
    );

    const G4double coldFingerRodRadius = coldFingerRodDiameter / 2.0;
    const G4double rodTopZ = blockBottomZ;
    const G4double rodBottomZ = rodTopZ - coldFingerRodHeight;
    const G4double rodCenterZ = (rodTopZ + rodBottomZ) / 2.0;

    auto* coldFingerRodSolid = new G4Tubs(
        "ColdFingerRodSolid",
        0.0,
        coldFingerRodRadius,
        coldFingerRodHeight / 2.0,
        0.0,
        twopi
    );

    auto* plateAndBlockSolid = new G4UnionSolid(
        "PlateAndBlockSolid",
        plateSolid,
        coldFingerBlockSolid,
        nullptr,
        G4ThreeVector(
            0.0,
            0.0,
            blockCenterZ - plateCenterZ + kJoinOverlap
        )
    );


    auto* plateBlockRodSolid = new G4UnionSolid(
        "PlateBlockRodSolid",
        plateAndBlockSolid,
        coldFingerRodSolid,
        nullptr,
        G4ThreeVector(
            0.0,
            0.0,
            rodCenterZ - plateCenterZ + 2.0 * kJoinOverlap
        )
    );

    const G4double ptfeTopRingHeight = ptfeSleeveThickness;
    const G4double holderInnerRadius = ptfeOuterRadius;
    const G4double holderOuterRadius =
        holderInnerRadius + housingSideThickness;
    const G4double holderTopZ =
        crystalTopZ + ptfeTopRingHeight + kJoinOverlap;
    const G4double holderBottomZ = plateTopZ - kJoinOverlap;
    const G4double holderHeight = holderTopZ - holderBottomZ;
    const G4double holderCenterZ = (holderTopZ + holderBottomZ) / 2.0;

    auto* holderSideSolid = MakeCylindricalShell(
        "HolderSideSolid",
        holderInnerRadius,
        holderOuterRadius,
        holderHeight
    );

    auto* completeAluminiumSolid = new G4UnionSolid(
        "CompleteAluminiumSolid",
        plateBlockRodSolid,
        holderSideSolid,
        nullptr,
        G4ThreeVector(0.0, 0.0, holderCenterZ - plateCenterZ)
    );

    auto* completeAluminiumLogical = new G4LogicalVolume(
        completeAluminiumSolid,
        aluminium,
        "CompleteAluminiumLogical"
    );
    new G4PVPlacement(
        nullptr,
        G4ThreeVector(0.0, 0.0, plateCenterZ),
        completeAluminiumLogical,
        "CompleteAluminium",
        worldLogical,
        false,
        0,
        true
    );

    const G4double ptfeTopRingInnerRadius = Q / 2.0;
    const G4double ptfeTopRingOuterRadius = ptfeOuterRadius;

    auto* ptfeTopRingSolid = MakeCylindricalShell(
        "PTFETopRingSolid",
        ptfeTopRingInnerRadius,
        ptfeTopRingOuterRadius,
        ptfeTopRingHeight
    );
    auto* ptfeTopRingLogical = new G4LogicalVolume(
        ptfeTopRingSolid,
        ptfe,
        "PTFETopRingLogical"
    );
    new G4PVPlacement(
        nullptr,
        G4ThreeVector(
            0.0,
            0.0,
            crystalTopZ + ptfeTopRingHeight / 2.0
        ),
        ptfeTopRingLogical,
        "PTFETopRing",
        worldLogical,
        false,
        0,
        true
    );

    auto* activeGeVis = new G4VisAttributes(G4Colour(0.2, 0.75, 0.9));
    activeGeVis->SetVisibility(true);
    activeGeVis->SetForceSolid(true);
    activeCrystalLogical->SetVisAttributes(activeGeVis);

    auto* deadGeVis = new G4VisAttributes(G4Colour(0.85, 0.25, 0.25));
    deadGeVis->SetVisibility(true);
    deadGeVis->SetForceSolid(true);
    deadCrystalLogical->SetVisAttributes(deadGeVis);

    auto* carbonVis = new G4VisAttributes(G4Colour(0.25, 0.25, 0.25));
    carbonVis->SetVisibility(true);
    carbonVis->SetForceSolid(true);
    windowLogical->SetVisAttributes(carbonVis);

    auto* aluminiumVis = new G4VisAttributes(G4Colour(0.75, 0.75, 0.78));
    aluminiumVis->SetVisibility(true);
    aluminiumVis->SetForceSolid(true);
    capAssemblyLogical->SetVisAttributes(aluminiumVis);
    carbonWindowAluminiumCoverLogical->SetVisAttributes(aluminiumVis);
    completeAluminiumLogical->SetVisAttributes(aluminiumVis);

    auto* ptfeVis = new G4VisAttributes(G4Colour(0.9, 0.8, 0.35));
    ptfeVis->SetVisibility(true);
    ptfeVis->SetForceSolid(true);
    ptfeSleeveLogical->SetVisAttributes(ptfeVis);
    ptfeTopRingLogical->SetVisAttributes(ptfeVis);

    return worldPhysical;
}

