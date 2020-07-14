import pytest
from adduct_calculator.adduct_rules import AdductTransformer
from adduct_calculator.adduct_rules import PROTON_MASS

def test_mass2ion_H():
    at = AdductTransformer()
    mz = 100.
    adduct = '[M+H]+'
    ion = at.mass2ion(mz,adduct)
    assert ion == pytest.approx(mz + PROTON_MASS,0.0001)

def test_mass2ion_K():
    at = AdductTransformer()
    mz = 100.
    adduct = '[M+K]+'
    ion = at.mass2ion(mz,adduct)
    assert ion == pytest.approx(mz + 38.9631581000907,0.0001)

def test_mass2ion_M_2H():
    at = AdductTransformer()
    mz = 100.
    adduct = '[M+2H]2+'
    ion = at.mass2ion(mz,adduct)
    assert ion == pytest.approx((mz + 2*PROTON_MASS)/2,0.0001)