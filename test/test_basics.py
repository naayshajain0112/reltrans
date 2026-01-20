import os

import pytest
import numpy as np
from wrapper import DCP_Parameters


def _debug_plot(energy, output, title=""):
    """Used when creating new tests for quickly looking at the data to make
    sure it is sensible."""
    import matplotlib.pyplot as plt

    plt.clf()
    plt.plot(energy[0:-1], output)
    plt.xscale("log")
    plt.show()


def test_basic_invocation(reltrans, assert_snapshot):
    """A smoke test to check whether the default values are working."""
    energy = np.logspace(np.log10(0.1), np.log10(100), 501)
    output = reltrans.dcp(energy, DCP_Parameters())
    assert_snapshot(output)


def test_re_im_parameter(reltrans, assert_snapshot, telescope):
    """Test the re_im parameter to assert that all the different outputs of the
    model are working. This test requires an RMF and ARF, which is provided by
    the `telescope` fixture."""
    energy = np.logspace(np.log10(0.1), np.log10(100), 101)

    os.environ["RMF_SET"] = telescope.rmf_path
    os.environ["ARF_SET"] = telescope.arf_path
    os.environ["EMIN_REF"] = "0.3"
    os.environ["EMAX_REF"] = "10.0"

    xrb1 = DCP_Parameters(mass=10.0, flo_hz=0.122, fhi_hz=0.224, re_im=4.0)
    output = reltrans.dcp(energy, xrb1)
    assert_snapshot(output, name="time_lag")

    xrb1.re_im = 1
    output = reltrans.dcp(energy, xrb1)
    assert_snapshot(output, name="real_part")

    xrb1.re_im = 2
    output = reltrans.dcp(energy, xrb1)
    assert_snapshot(output, name="imaginary_part")

    xrb1.re_im = 3
    output = reltrans.dcp(energy, xrb1)
    assert_snapshot(output, name="magnitude")
