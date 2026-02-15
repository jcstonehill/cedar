import time

import cedar.benchmarks.heat_transfer as heat_transfer

import cedar

_benchmarks: dict[str, list[cedar.Benchmark]] = {
    "Heat Transfer" : [
        heat_transfer.AdiabaticTransient0D,
        heat_transfer.DirichletSteady1D,
        heat_transfer.NAFEMSTransient1D
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