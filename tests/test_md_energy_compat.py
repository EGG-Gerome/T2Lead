"""Regression tests for MD energy return-value type compatibility.

Covers the bug where ``getPotentialEnergy()`` sometimes returns a plain
``float`` instead of an OpenMM ``Quantity``, causing
``'float' object has no attribute 'value_in_unit'``.
"""

from __future__ import annotations

import math
import pytest

from drugpipe.lead_optimization.md_simulation import (
    _energy_as_kcal,
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
