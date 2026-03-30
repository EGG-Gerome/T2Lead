"""Regression tests for MD energy/positions type compatibility.

Covers bugs where OpenMM return values (energy, positions, element.mass)
are sometimes ``Quantity`` and sometimes plain Python types, causing
``'float' object has no attribute 'value_in_unit'``.
"""

from __future__ import annotations

import math
import numpy as np
import pytest

from drugpipe.lead_optimization.md_simulation import (
    _energy_as_kcal,
    _positions_as_angstrom,
    _is_heavy_atom,
    _KJ_TO_KCAL,
    _OPENMM_OK,
)

try:
    import openmm.unit as unit
except ImportError:
    unit = None


class TestEnergyAsKcal:
    """Ensure ``_energy_as_kcal`` handles both Quantity and float."""

    def test_plain_float_converted_from_kj(self):
        raw_kj = -1000.0
        result = _energy_as_kcal(raw_kj)
        expected = raw_kj * _KJ_TO_KCAL
        assert math.isclose(result, expected, rel_tol=1e-6)

    def test_plain_int(self):
        result = _energy_as_kcal(-500)
        expected = -500 * _KJ_TO_KCAL
        assert math.isclose(result, expected, rel_tol=1e-6)

    def test_zero(self):
        assert _energy_as_kcal(0.0) == 0.0

    def test_positive_energy(self):
        result = _energy_as_kcal(100.0)
        assert result > 0

    def test_no_attribute_error_on_float(self):
        """The original bug: calling .value_in_unit() on a float."""
        try:
            _energy_as_kcal(42.0)
        except AttributeError:
            pytest.fail("_energy_as_kcal raised AttributeError on float input")

    @pytest.mark.skipif(not _OPENMM_OK, reason="OpenMM not installed")
    def test_real_openmm_quantity_kj(self):
        q = -1234.5 * unit.kilojoules_per_mole
        result = _energy_as_kcal(q)
        expected = float(q.value_in_unit(unit.kilocalories_per_mole))
        assert math.isclose(result, expected, rel_tol=1e-6)

    @pytest.mark.skipif(not _OPENMM_OK, reason="OpenMM not installed")
    def test_real_openmm_quantity_kcal(self):
        q = -295.0 * unit.kilocalories_per_mole
        result = _energy_as_kcal(q)
        assert math.isclose(result, -295.0, rel_tol=1e-6)

    @pytest.mark.skipif(not _OPENMM_OK, reason="OpenMM not installed")
    def test_consistency_quantity_vs_float(self):
        """Quantity and manually-converted float should agree."""
        kj_val = -2000.0
        q = kj_val * unit.kilojoules_per_mole
        from_quantity = _energy_as_kcal(q)
        from_float = _energy_as_kcal(kj_val)
        assert math.isclose(from_quantity, from_float, rel_tol=1e-4)


class TestPositionsAsAngstrom:
    """Ensure ``_positions_as_angstrom`` handles array-like and Quantity."""

    def test_plain_numpy_array(self):
        arr = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        result = _positions_as_angstrom(arr)
        assert result.shape == (2, 3)
        assert result.dtype == np.float64

    def test_plain_list(self):
        lst = [[0.1, 0.2, 0.3]]
        result = _positions_as_angstrom(lst)
        assert result.shape == (1, 3)

    @pytest.mark.skipif(not _OPENMM_OK, reason="OpenMM not installed")
    def test_openmm_quantity(self):
        from openmm import Vec3
        pos = [Vec3(1, 2, 3)] * unit.angstrom
        result = _positions_as_angstrom(pos)
        assert math.isclose(result[0, 0], 1.0, rel_tol=1e-6)


class TestIsHeavyAtom:
    """Ensure ``_is_heavy_atom`` avoids Quantity comparison pitfalls."""

    def test_fake_carbon(self):
        class FakeElem:
            symbol = "C"
            mass = 12.0
        class FakeAtom:
            element = FakeElem()
        assert _is_heavy_atom(FakeAtom()) is True

    def test_fake_hydrogen(self):
        class FakeElem:
            symbol = "H"
            mass = 1.008
        class FakeAtom:
            element = FakeElem()
        assert _is_heavy_atom(FakeAtom()) is False

    def test_no_element(self):
        class FakeAtom:
            element = None
        assert _is_heavy_atom(FakeAtom()) is False

    def test_no_attribute_error_on_float_mass(self):
        """Regression: Quantity(mass).__gt__(1.5) calls float.value_in_unit()."""
        class FakeElem:
            symbol = "N"
            mass = 14.0
        class FakeAtom:
            element = FakeElem()
        try:
            _is_heavy_atom(FakeAtom())
        except AttributeError:
            pytest.fail("_is_heavy_atom raised AttributeError")
