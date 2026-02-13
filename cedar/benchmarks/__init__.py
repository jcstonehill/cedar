import time

from cedar.benchmarks.heat_transfer.adiabatic_transient_0d import AdiabaticTransient0D
from cedar.benchmarks.heat_transfer.dirichlet_steady_1d import DirichletSteady1D
from cedar.benchmarks.heat_transfer.nafems_transient_1d import NAFEMSTransient1D

import cedar

_benchmarks: dict[str, list[cedar.Benchmark]] = {
    "Heat Transfer" : [
        AdiabaticTransient0D,
        DirichletSteady1D,
        NAFEMSTransient1D
    ]
}

def run():
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
            
            cedar.Log.message(f"{benchmark.__name__:<20}     {status:^9}     {f"Error: {error:.1f} [%]"}     {cedar.functions.format_computation_time(end-start):30}")

        cedar.Log.remove_level()