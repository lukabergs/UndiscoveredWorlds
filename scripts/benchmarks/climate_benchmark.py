"""Run one registered Earth climate benchmark and write its gallery report.

The native flags are those of `world_pipeline climate-benchmark`. This module
adds the benchmark defaults, captures the log under runs/logs/climate/<id>.log
and calls write_benchmark_report.py for the new run.
"""
import hashlib
import os
from pathlib import Path
import re
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from uwtools import paths
import run_log

# `uv run` exports its own path; subprocesses may not find uv on PATH.
UV = os.environ.get('UV', 'uv')
ROOT = Path(__file__).resolve().parents[2]
EXECUTABLE = paths.executable('world_pipeline')
DEFAULT_MAPS = ('koppen,temp,precip,precip-tif,jan-s-wind-speed,jul-s-wind-speed,jan-u-wind-speed,jul-u-wind-speed,'
    'jan-s-wind-lic,jan-s-wind-part,jul-s-wind-lic,jul-s-wind-part,jan-u-wind-lic,jan-u-wind-part,jul-u-wind-lic,'
    'jul-u-wind-part,jan-sst,jul-sst,jan-ocean-current,jul-ocean-current')
# Run 373 grids: 1024x512 world, 128-wide atmosphere/ocean, 256-wide moisture.
DEFAULTS = {'--seed': '1', '--resolution': '1024', '--climate-resolution': '128',
            '--climate-hydrology-resolution': '256', '--climate-coupling-iterations': '8', '--map': DEFAULT_MAPS}
FIXED = ['--no-rivers', '--no-lakes', '--no-deltas']
# Earth benchmark input assist used by runs up to 378: ERA5 vegetation bulk
# drag contrast (see its metadata.json). Not used for procedural worlds.
DRAG_DIRECTORY = ROOT / 'refs/processed/climate/vegetation-drag-20260914'
DRAG_SOURCE = DRAG_DIRECTORY / 'land-drag-contrast-1024.txt'
DRAG_VARIABLE = 'UW_CLIMATE_SURFACE_DRAG_INPUT'
# The console shows stage progress, warnings and the result; solver diagnostics
# (and everything else) stay in runs/logs/climate/<id>.log, which reports parse.
PROGRESS = re.compile(r'(\.\.\.$|\[[0-9.]+ms\]$|^Generating |^Climate benchmark run ID|did not converge'
                      r'|^World generation completed)')
PROBLEM = re.compile(r'\b(error|warning|failed)\b', re.IGNORECASE)
DIAGNOSTIC = re.compile(r'\w=\S')


def console_line(line):
    line = line.strip()
    if DIAGNOSTIC.search(line):
        return False
    return bool(PROGRESS.search(line) or PROBLEM.search(line))

# Parameter names of the retired run-climate-benchmark.ps1, kept so archived
# command.json receipts can be replayed.
LEGACY_PARAMETERS = {
    'Resolution': '--resolution', 'WindMapResolution': '--wind-map-resolution',
    'ClimateResolution': '--climate-resolution', 'HydrologyResolution': '--climate-hydrology-resolution',
    'CouplingIterations': '--climate-coupling-iterations', 'Seed': '--seed', 'Maps': '--map',
    'TemperatureMode': '--climate-temperature-mode', 'ClimateMode': '--climate-mode',
    'Transport': '--climate-transport', 'TropicalClosure': '--climate-tropical-closure',
    'OceanStepHours': '--climate-ocean-step-hours', 'OceanCirculationCache': '--climate-ocean-circulation-cache',
    'OceanCirculationSolver': '--climate-ocean-circulation-solver', 'OceanHeatKernel': '--climate-ocean-heat-kernel',
    'OceanHeatWorkers': '--climate-ocean-heat-workers', 'OceanCoastalScheme': '--climate-ocean-coastal-scheme',
    'OceanCapture': '--climate-ocean-capture', 'AtmosphereSolver': '--climate-atmosphere-solver',
    'AtmosphereCapture': '--climate-atmosphere-capture', 'Executable': '--executable'}


def bench_command(arguments):
    """Return the uw.py command list for native climate arguments."""
    return [sys.executable, str(ROOT / 'scripts/uw.py'), 'climate', 'run', *arguments]


def from_legacy_command(command):
    """Convert a run-climate-benchmark.ps1 receipt to uw.py arguments; pass others through."""
    if not any(str(part).endswith('run-climate-benchmark.ps1') for part in command):
        return command
    start = next(i for i, part in enumerate(command) if str(part).endswith('run-climate-benchmark.ps1')) + 1
    arguments = []
    for name, value in zip(command[start::2], command[start + 1::2]):
        flag = LEGACY_PARAMETERS[name.lstrip('-')]
        if flag == '--climate-hydrology-resolution' and str(value) == '0':
            continue
        if flag == '--wind-map-resolution':
            continue  # retired: maps now follow the world grid
        arguments += [flag, str(value)]
    return bench_command(arguments)


def set_option(command, flag, value):
    """Set or append one flag value in a uw.py command list."""
    if flag in command:
        command[command.index(flag) + 1] = str(value)
    else:
        command += [flag, str(value)]
    return command


def surface_drag_input(width):
    """Return the vegetation drag contrast for a `width` world, area-averaging the 1024 source once if needed."""
    path = DRAG_DIRECTORY / f'land-drag-contrast-{width}.txt'
    if path.exists():
        return path
    import numpy as np
    sys.path.insert(0, str(ROOT / 'scripts/refs'))
    from cell_grid import conservative_remap
    values = DRAG_SOURCE.read_text(encoding='utf-8').split()
    if values[0] != 'UWROUGH1':
        raise SystemExit(f'Unexpected drag input header: {DRAG_SOURCE}')
    columns, rows = int(values[1]), int(values[2])
    source = np.asarray(values[3:], dtype=np.float64).reshape(rows, columns)
    # Stress factors are remapped by area, as in the 1024 preparation.
    field = np.clip(conservative_remap(source, width), 0.25, 4.0)
    temporary = path.with_suffix('.tmp')
    with temporary.open('w', encoding='utf-8') as stream:
        stream.write(f'UWROUGH1 {width} {width // 2}\n')
        for row in field:
            stream.write(' '.join(f'{value:.6g}' for value in row) + '\n')
    temporary.replace(path)
    print(f'Derived {path.name} from {DRAG_SOURCE.name} by conservative area remap.', flush=True)
    return path


def native_arguments(arguments, label=''):
    arguments = list(arguments)
    for flag, value in DEFAULTS.items():
        if flag not in arguments and not (flag == '--map' and '--maps' in arguments):
            arguments += [flag, value]
    arguments += [flag for flag in FIXED if flag not in arguments]
    pairs = ' '.join(arguments)
    information = f'{label}; {pairs}' if label else pairs
    return arguments + ['--benchmark-info', information], information


def run(arguments, executable=EXECUTABLE, label='', surface_drag='earth'):
    """Run the benchmark and its report; return the registered run ID.

    surface_drag: 'earth' (vegetation drag contrast for the world width),
    'none' (uniform land drag) or a UWROUGH1 file path. A value already set in
    the environment (as climate_experiment.py does) takes precedence.
    """
    executable = Path(executable)
    if not executable.is_file():
        raise SystemExit(f'Build world_pipeline first: {executable}')
    # Archived pre-split executables take the old flag instead of a command.
    command = ['--earth-climate-benchmark'] if executable.name == 'UndiscoveredWorlds.exe' else ['climate-benchmark']
    width = int(DEFAULTS['--resolution'])
    if '--resolution' in arguments:
        width = int(arguments[arguments.index('--resolution') + 1])
    environment = os.environ.copy()
    if not environment.get(DRAG_VARIABLE):
        if surface_drag == 'none':
            environment.pop(DRAG_VARIABLE, None)
        else:
            drag = surface_drag_input(width) if surface_drag == 'earth' else Path(surface_drag).resolve()
            environment[DRAG_VARIABLE] = str(drag)
    drag_label = Path(environment[DRAG_VARIABLE]).name if environment.get(DRAG_VARIABLE) else 'uniform'
    native, information = native_arguments(arguments, f'{label}; surface drag {drag_label}' if label else f'surface drag {drag_label}')
    log_root = ROOT / 'runs/logs/climate'
    log_root.mkdir(parents=True, exist_ok=True)
    pending = log_root / f'pending-{time.strftime("%Y%m%d-%H%M%S")}-{time.perf_counter_ns() % 1000000}.log'
    print(f'Benchmark: {information}', flush=True)
    start = time.perf_counter()
    with pending.open('x', encoding='utf-8') as log:
        log.write(run_log.header('climate', executable, [executable, *command, *native]))
        log.write(f'{DRAG_VARIABLE}: {environment.get(DRAG_VARIABLE, "(unset)")}\n')
        process = subprocess.Popen([str(executable), *command, *native], cwd=ROOT, env=environment, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT, text=True, encoding='utf-8', errors='replace')
        for line in process.stdout:
            log.write(line)
            log.flush()
            if console_line(line):
                print(line, end='', flush=True)
        exit_code = process.wait()
    elapsed = time.perf_counter() - start
    if exit_code:
        raise SystemExit(f'Benchmark exited {exit_code}; log: {pending}')
    ids = re.findall(r'Climate benchmark run ID: (\d+)', pending.read_text(encoding='utf-8'))
    if not ids:
        raise SystemExit(f'No completed run ID in {pending}')
    run_id = int(ids[-1])
    final = log_root / f'{run_id}.log'
    if final.exists():
        raise SystemExit(f'Refusing to replace existing log: {final}')
    pending.rename(final)
    subprocess.run([UV, 'run', '--offline', 'python', str(ROOT / 'scripts/benchmarks/write_benchmark_report.py'), '--run-id', str(run_id),
                    '--elapsed-seconds', f'{elapsed:.7f}'], cwd=ROOT, check=True)
    print(f'Run {run_id} completed in {elapsed:.1f} s. Log: {final}', flush=True)
    return run_id
