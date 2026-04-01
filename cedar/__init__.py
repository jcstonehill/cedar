from cedar.log import Log
from cedar.mesh import *
from cedar.field import *
from cedar.transfer import Transfer
from cedar.term import Term
from cedar.fluid import Fluid
from cedar.material import Material
from cedar.model import Model
from cedar.problem import Problem
from cedar.benchmark import Benchmark

# Fluids
import cedar.fluids as fluids

# Functions
import cedar.functions as functions

# Materials
import cedar.materials as materials

# Models
from cedar.models.flow.bc import *
from cedar.models.flow.source import *
from cedar.models.flow.flow import *

from cedar.models.heat_transfer.bc import *
from cedar.models.heat_transfer.source import *
from cedar.models.heat_transfer.heat_transfer import *

# Benchmarks
import cedar.models.heat_transfer.verification as _heat_transfer_verification

# Benchmarks
_benchmarks: dict[str, list[cedar.Benchmark]] = {
    "Heat Transfer": [
        _heat_transfer_verification.AdiabaticTransient0D,
        _heat_transfer_verification.DirichletSteady1D,
        _heat_transfer_verification.NAFEMSTransient1D,
    ]
}

import time


def run_benchmarks():
    cedar.Log.message("Running Benchmarks!")
    cedar.Log.line_break()

    for section in _benchmarks:
        cedar.Log.message(section)
        cedar.Log.add_level()

        for benchmark in _benchmarks[section]:
            # Computation takes place here.
            cedar.Log.is_active = False
            start = time.time()
            error = benchmark.run()
            end = time.time()
            cedar.Log.is_active = True

            status = "Passed" if error <= benchmark.tol else "FAILED!!!"

            cedar.Log.message(
                f"{benchmark.__name__:<20}     {status:^9}     {f"Error: {error:.1f} [%]"}     {cedar.functions.format_computation_time(end-start):30}"
            )

        cedar.Log.remove_level()
