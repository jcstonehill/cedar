import time
import os, shutil

import cedar.benchmarks.flow as flow
import cedar.benchmarks.thermal as thermal
import cedar

_benchmarks: dict[list[cedar.base.Benchmark]] = {
    "Flow" : [
        flow.THMCase1
    ],
    
    "Thermal" : [
        thermal.AdiabaticTransient0D,
        thermal.DirichletSteady1D,
        thermal.NAFEMSTransient1D,
    ]
}

def plot():
    """
    Generate and save plots for all registered benchmarks.

    This function iterates over all benchmark sections and benchmarks,
    calls each benchmark's ``plot`` method, and saves the resulting
    figures to the ``outputs/benchmarks/`` directory.

    Notes
    -----
    - All plots are saved to ``outputs/benchmarks/<BenchmarkName>/``.
    - Existing files and directories under ``outputs/benchmarks`` are
      deleted before plotting begins.
    """
    cedar.Log.print_messages = True
    cedar.Log.message("Plotting Benchmarks!")

    output_path = "outputs/benchmarks"

    if(os.path.isdir(output_path)):
        for item in os.listdir(output_path):
            item_path = os.path.join(output_path, item)
            if os.path.isfile(item_path):
                os.remove(item_path)

            elif os.path.isdir(item_path):
                shutil.rmtree(item_path)

    else:
        os.mkdir("outputs/benchmarks")

    for section in _benchmarks.keys():
        cedar.Log.message(section)
        cedar.Log.add_level()

        benchmark: cedar.base.Benchmark
        for benchmark in _benchmarks[section]:
            os.mkdir("outputs/benchmarks/" + benchmark.__name__)
            cedar.Log.print_messages = False
            benchmark.plot(show = False, save = True)
            cedar.Log.print_messages = True

            cedar.Log.message(benchmark.__name__)

        cedar.Log.remove_level()

    cedar.Log.message("All plots saved to /outputs/benchmarks/")


def run():
    """
    Run all registered benchmarks and report their results.

    This function executes each benchmark's ``run`` method, compares the
    resulting error against the benchmark tolerance, and reports pass or
    fail status along with execution time.

    Benchmarks are grouped by section, and a formatted summary is printed
    to the console. Failed benchmarks are tracked internally for reporting.

    Notes
    -----
    - A benchmark is considered **passed** if ``error <= benchmark.tol``.
    - Execution time is measured for each benchmark individually.
    - Logging output is temporarily suppressed during benchmark execution
      to reduce console noise.
    """
    cedar.Log.print_messages = True
    cedar.Log.message("Running Benchmarks!")
    cedar.Log.line_break()

    N_failed = 0
    N_passed = 0

    failed_tests = ""
    cedar.Log.add_level()
    for section in _benchmarks.keys():
        cedar.Log.message(section)
        cedar.Log.add_level()

        benchmark: cedar.base.Benchmark
        for benchmark in _benchmarks[section]:
            name = benchmark.__name__

            cedar.Log.print_messages = False
            start = time.time()
            error = benchmark.run()
            end = time.time()
            cedar.Log.print_messages = True

            if error <= benchmark.tol:
                status = "Passed"
                N_passed += 1

            else:
                status = "FAILED!!!"
                N_failed += 1

                if not failed_tests == "":
                    failed_tests += ", "
                failed_tests += section + "." + name

            error_str = f"Error: {error:.1f} [%]"
            cedar.Log.message(f"{name:>30}     {status:^9}     {error_str}     {cedar.helper.format_computation_time(end-start):30}")

        cedar.Log.line_break()
        cedar.Log.remove_level()
    cedar.Log.remove_level()