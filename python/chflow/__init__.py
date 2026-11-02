"""ChuanFlow Python package (README M1/M2 binding layer).

The heavy lifting lives in the C++ kernel; this package re-exports the
pybind11 module `chflow_core` (built by the ChuanFlow CMake tree) with a
Python API. Put both this directory's parent and the build output
directory on PYTHONPATH (CTest does this for the regression test).
"""

from chflow_core import (  # noqa: F401
    FluidProperties,
    NetworkPipe,
    OneDNetwork,
    PipeFlowResult,
    SteadyControls,
    SteadyNetworkResult,
    TransientCase,
    TransientNetwork,
    TransientNodeConfig,
    TransientPipeConfig,
    TransientSample,
    TransientSeries,
    TransientSource,
    TransientTimeControls,
    churchill_friction,
    hydrate_kinetics_step,
    load_network,
    load_steady_controls,
    load_transient_case,
    load_transient_controls,
    run_transient,
    solve_steady_network,
)

__all__ = [
    "FluidProperties",
    "NetworkPipe",
    "OneDNetwork",
    "SteadyControls",
    "PipeFlowResult",
    "SteadyNetworkResult",
    "TransientCase",
    "TransientNetwork",
    "TransientNodeConfig",
    "TransientPipeConfig",
    "TransientSample",
    "TransientSeries",
    "TransientSource",
    "TransientTimeControls",
    "churchill_friction",
    "hydrate_kinetics_step",
    "load_network",
    "load_steady_controls",
    "load_transient_case",
    "load_transient_controls",
    "run_transient",
    "solve_steady_network",
]
