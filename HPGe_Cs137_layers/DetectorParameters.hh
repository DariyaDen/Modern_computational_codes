#ifndef DetectorParameters_h
#define DetectorParameters_h 1

#include "G4SystemOfUnits.hh"
#include "globals.hh"

namespace DetectorParameters
{

    inline constexpr G4double D = 71.0 * mm;
    inline constexpr G4double H = 54.9 * mm;
    inline constexpr G4double Q = 66.0 * mm;
    inline constexpr G4double T = 2.9 * mm;

    inline constexpr G4double P = 9.2 * mm;
    inline constexpr G4double centralHoleDepth = 40.4 * mm;

    inline constexpr G4double M = 5.9 * mm;
    inline constexpr G4double S1 = 16.5 * mm;
    inline constexpr G4double S2 = 27.0 * mm;

    inline constexpr G4double d1 = 0.15 * mm;
    inline constexpr G4double d2 = 1.0 * mm;
    inline constexpr G4double d3 = 1.0 * mm;

    inline constexpr G4double d4 = 0.0003 * mm;

    inline constexpr G4double d5 = 0.0003 * mm;

    inline constexpr G4double d6 = 3.6496 * mm;

    inline constexpr G4double d7 = 1.0 * mm;

    inline constexpr G4double d8 = 0.0003 * mm;


    inline constexpr G4double capMaxDiameter = 95.0 * mm;
    inline constexpr G4double capStepDiameter = 85.0 * mm;
    inline constexpr G4double capInnerDiameter = 74.0 * mm;

    inline constexpr G4double carbonDiameter = 85.0 * mm;
    inline constexpr G4double thinWindowThickness = 0.15 * mm;

    inline constexpr G4double windowGapAboveCrystal = 7.0 * mm;

    inline constexpr G4double topAssemblyHeight = 1.0 * mm;

    inline constexpr G4double deadLayerThickness = 0.7 * mm;
    inline constexpr G4double housingSideThickness = 1.5 * mm;
    inline constexpr G4double housingBottomThickness = 1.5 * mm;

    inline constexpr G4double ptfeSleeveThickness = 1.0 * mm;

    inline constexpr G4double plateDiameter = 80.0 * mm;
    inline constexpr G4double plateThickness = 8.0 * mm;
    inline constexpr G4double plateTopOffsetFromCrystalTop = 61.5 * mm;

    inline constexpr G4double coldFingerBlockDiameter = 27.0 * mm;
    inline constexpr G4double coldFingerBlockHeight = 31.0 * mm;
    inline constexpr G4double coldFingerRodDiameter = 11.0 * mm;
    inline constexpr G4double coldFingerRodHeight = 28.0 * mm;


    inline constexpr G4double worldHalfSize = 0.5 * m;
    inline constexpr G4double sourceDistance = 250.0 * mm;
    inline constexpr G4double sourceEnergy = 661.657 * keV;
}

#endif
