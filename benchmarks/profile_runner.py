#!/usr/bin/env python3
"""
profile_runner.py - Memory (Peak RSS) and Time Profiling Runner
Part of Pompeiu-Hausdorff distance research and benchmarking.

Captures:
- Peak RSS resident memory (monitoring /proc/[pid]/status @ 10ms + getrusage)
- Wall-clock time and CPU time (user + system) with millisecond precision
- Subdivided face count and upper bound calls (u0, u1, u2, u3, u4)
- Exports structured metrics to CSV (e.g. results_preliminar.csv)
"""

import argparse
import csv
import datetime
import multiprocessing
import os
import re
import resource
import subprocess
import sys
import threading
import time
from typing import Dict, Optional, Tuple


class ProcStatusSampler(threading.Thread):
    """
    Background sampler thread monitoring /proc/[pid]/status at high frequency (10ms).
    Tracks VmHWM (Peak Resident Set Size) and instantaneous VmRSS.
    """

    def __init__(self, pid: int, sample_interval_sec: float = 0.010):
        super().__init__(daemon=True)
        self.pid = pid
        self.interval = sample_interval_sec
        self.stop_event = threading.Event()
        self.peak_vmhwm_kb = 0.0
        self.peak_vmrss_kb = 0.0

    def run(self):
        status_path = f"/proc/{self.pid}/status"
        while not self.stop_event.is_set():
            try:
                if not os.path.exists(status_path):
                    break
                with open(status_path, "r", encoding="utf-8") as f:
                    for line in f:
                        if line.startswith("VmHWM:"):
                            parts = line.split()
                            if len(parts) >= 2:
                                val = float(parts[1])
                                if val > self.peak_vmhwm_kb:
                                    self.peak_vmhwm_kb = val
                        elif line.startswith("VmRSS:"):
                            parts = line.split()
                            if len(parts) >= 2:
                                val = float(parts[1])
                                if val > self.peak_vmrss_kb:
                                    self.peak_vmrss_kb = val
            except (ProcessLookupError, FileNotFoundError, PermissionError):
                break
            except Exception:
                pass
            time.sleep(self.interval)

    def stop(self):
        self.stop_event.set()


def parse_binary_stdout(stdout_text: str) -> Dict[str, float]:
    """
    Parses key numerical results and profiling metrics emitted by the binary.
    """
    metrics: Dict[str, float] = {}

    patterns = {
        "dA": r"dA\s*=\s*([0-9.eE+-]+)",
        "lower": r"lower\s*=\s*([0-9.eE+-]+)",
        "upper_max": r"upper_max\s*=\s*([0-9.eE+-]+)",
        "lower_over_dA": r"lower/dA\s*=\s*([0-9.eE+-]+)",
        "upper_over_dA": r"upper_max/dA\s*=\s*([0-9.eE+-]+)",
        "bvh_time_ms": r"bvh_time\(ms\)\s*=\s*([0-9.eE+-]+)",
        "bound_time_ms": r"bound_time\(ms\)\s*=\s*([0-9.eE+-]+)",
        "subdivided_faces": r"subdivided_faces\s*=\s*([0-9]+)",
        "u0_calls": r"u0_calls\s*=\s*([0-9]+)",
        "u1_calls": r"u1_calls\s*=\s*([0-9]+)",
        "u2_calls": r"u2_calls\s*=\s*([0-9]+)",
        "u3_calls": r"u3_calls\s*=\s*([0-9]+)",
        "u4_calls": r"u4_calls\s*=\s*([0-9]+)",
    }

    for key, pattern in patterns.items():
        match = re.search(pattern, stdout_text)
        if match:
            try:
                metrics[key] = float(match.group(1))
            except ValueError:
                metrics[key] = float("nan")

    return metrics


def run_single_profile(
    executable: str,
    mesh_a: str,
    mesh_b: str,
    tol: str,
    max_factor: str,
    normalize: str,
    engine: str,
    sample_interval_sec: float = 0.010,
    verbose: bool = True,
) -> Dict[str, any]:
    """
    Executes one profiled execution of the pompeiu_hausdorff binary.
    """
    cmd = [executable, mesh_a, mesh_b, tol, max_factor, normalize, engine]
    if verbose:
        print(f"\n[profile_runner] Launching: {' '.join(cmd)}")

    ru_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    start_wall = time.perf_counter()

    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    sampler = ProcStatusSampler(process.pid, sample_interval_sec=sample_interval_sec)
    sampler.start()

    stdout, stderr = process.communicate()
    end_wall = time.perf_counter()

    sampler.stop()
    sampler.join(timeout=0.5)

    ru_after = resource.getrusage(resource.RUSAGE_CHILDREN)

    # Wall-clock time in milliseconds
    wall_time_ms = (end_wall - start_wall) * 1000.0

    # CPU time in milliseconds (User + System)
    user_cpu_sec = max(0.0, ru_after.ru_utime - ru_before.ru_utime)
    sys_cpu_sec = max(0.0, ru_after.ru_stime - ru_before.ru_stime)
    cpu_time_ms = (user_cpu_sec + sys_cpu_sec) * 1000.0

    # Peak RSS from /proc (kB -> MB)
    proc_peak_kb = max(sampler.peak_vmhwm_kb, sampler.peak_vmrss_kb)
    peak_rss_proc_mb = proc_peak_kb / 1024.0

    # Peak RSS from getrusage (Linux reports ru_maxrss in kB -> MB)
    rusage_peak_kb = float(ru_after.ru_maxrss)
    peak_rss_rusage_mb = rusage_peak_kb / 1024.0

    # Global peak RSS takes the maximum observed
    peak_rss_mb = max(peak_rss_proc_mb, peak_rss_rusage_mb)

    parsed_metrics = parse_binary_stdout(stdout)

    record = {
        "timestamp": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "mesh_a": os.path.basename(mesh_a),
        "mesh_b": os.path.basename(mesh_b),
        "tol": tol,
        "max_factor": max_factor,
        "normalize": normalize,
        "engine": engine,
        "wall_time_ms": round(wall_time_ms, 3),
        "cpu_time_ms": round(cpu_time_ms, 3),
        "peak_rss_mb": round(peak_rss_mb, 3),
        "peak_rss_proc_mb": round(peak_rss_proc_mb, 3),
        "peak_rss_rusage_mb": round(peak_rss_rusage_mb, 3),
        "subdivided_faces": int(parsed_metrics.get("subdivided_faces", 0)),
        "u0_calls": int(parsed_metrics.get("u0_calls", 0)),
        "u1_calls": int(parsed_metrics.get("u1_calls", 0)),
        "u2_calls": int(parsed_metrics.get("u2_calls", 0)),
        "u3_calls": int(parsed_metrics.get("u3_calls", 0)),
        "u4_calls": int(parsed_metrics.get("u4_calls", 0)),
        "lower": parsed_metrics.get("lower", float("nan")),
        "upper_max": parsed_metrics.get("upper_max", float("nan")),
        "dA": parsed_metrics.get("dA", float("nan")),
        "lower_over_dA": parsed_metrics.get("lower_over_dA", float("nan")),
        "upper_over_dA": parsed_metrics.get("upper_over_dA", float("nan")),
        "bvh_time_ms": parsed_metrics.get("bvh_time_ms", float("nan")),
        "bound_time_ms": parsed_metrics.get("bound_time_ms", float("nan")),
        "return_code": process.returncode,
    }

    if verbose:
        print("-" * 50)
        print(f"Return Code:       {record['return_code']}")
        print(f"Wall-Clock Time:   {record['wall_time_ms']:.2f} ms")
        print(f"CPU Time:          {record['cpu_time_ms']:.2f} ms")
        print(f"Peak RSS (Total):  {record['peak_rss_mb']:.2f} MB")
        print(f"  /proc/status:    {record['peak_rss_proc_mb']:.2f} MB")
        print(f"  getrusage:       {record['peak_rss_rusage_mb']:.2f} MB")
        print(f"Subdivided Faces:  {record['subdivided_faces']}")
        print(f"Upper Bound Calls: u0={record['u0_calls']}, u1={record['u1_calls']}, "
              f"u2={record['u2_calls']}, u3={record['u3_calls']}, u4={record['u4_calls']}")
        print(f"Distance Bounds:   [{record['lower']:.8f}, {record['upper_max']:.8f}]")
        print("-" * 50)

    if process.returncode != 1:
        details = stderr.strip() or stdout.strip()
        raise RuntimeError(
            f"Profiling failed (exit code {process.returncode}): {details}"
        )

    return record


def _worker_entrypoint(params: Dict[str, any], queue: multiprocessing.Queue):
    """
    Entrypoint executed inside a fresh worker process.
    Guarantees that getrusage(RUSAGE_CHILDREN) starts at 0 without residual metrics.
    """
    try:
        record = run_single_profile(**params)
        queue.put((True, record))
    except Exception as exc:
        queue.put((False, str(exc)))


def run_isolated_profile(**params) -> Dict[str, any]:
    """
    Executes a single run inside an isolated subprocess worker so that
    RUSAGE_CHILDREN reflects only this run's memory without bleeding into subsequent runs.
    """
    ctx = multiprocessing.get_context("spawn")
    queue = ctx.Queue()
    worker = ctx.Process(target=_worker_entrypoint, args=(params, queue))
    worker.start()

    # Read result from queue before joining to avoid pipe buffer deadlock
    result = None
    try:
        while worker.is_alive():
            try:
                result = queue.get(timeout=0.1)
                break
            except Exception:
                continue
        if result is None and not queue.empty():
            result = queue.get(timeout=1.0)
    finally:
        worker.join()

    if result is None:
        raise RuntimeError("Profile worker process exited unexpectedly without returning metrics")

    success, payload = result
    if not success:
        raise RuntimeError(payload)

    return payload


def append_record_to_csv(csv_path: str, record: Dict[str, any]) -> None:
    """
    Appends a profiled record to the specified CSV file, creating headers if needed.
    """
    fieldnames = [
        "timestamp",
        "mesh_a",
        "mesh_b",
        "tol",
        "max_factor",
        "normalize",
        "engine",
        "wall_time_ms",
        "cpu_time_ms",
        "peak_rss_mb",
        "peak_rss_proc_mb",
        "peak_rss_rusage_mb",
        "subdivided_faces",
        "u0_calls",
        "u1_calls",
        "u2_calls",
        "u3_calls",
        "u4_calls",
        "lower",
        "upper_max",
        "dA",
        "lower_over_dA",
        "upper_over_dA",
        "bvh_time_ms",
        "bound_time_ms",
        "return_code",
    ]

    file_exists = os.path.isfile(csv_path)

    os.makedirs(os.path.dirname(os.path.abspath(csv_path)), exist_ok=True)

    with open(csv_path, mode="a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if not file_exists or os.path.getsize(csv_path) == 0:
            writer.writeheader()
        writer.writerow(record)


def main():
    parser = argparse.ArgumentParser(
        description="Profile execution memory (Peak RSS) and runtime for Pompeiu-Hausdorff distance."
    )
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    default_exe = os.path.join(repo_root, "build", "pompeiu_hausdorff")
    default_mesh_a = os.path.join(repo_root, "meshes", "107100.obj")
    default_mesh_b = os.path.join(repo_root, "meshes", "107100_sf.obj")
    default_csv = os.path.join(repo_root, "results_preliminar.csv")

    parser.add_argument("--executable", default=default_exe, help="Path to pompeiu_hausdorff executable")
    parser.add_argument("--mesh-a", default=default_mesh_a, help="Path to mesh A (.obj)")
    parser.add_argument("--mesh-b", default=default_mesh_b, help="Path to mesh B (.obj)")
    parser.add_argument("--tol", default="1e-4", help="Tolerance threshold (default: 1e-4)")
    parser.add_argument("--max-factor", default="1000000", help="Max factor for subdivision")
    parser.add_argument("--normalize", default="1", help="Normalize bounds by dA (0 or 1)")
    parser.add_argument(
        "--engine",
        default="best_first",
        choices=["best_first", "idastar"],
        help="Search engine to evaluate (best_first | idastar)",
    )
    parser.add_argument("--runs", type=int, default=1, help="Number of repetitions to run")
    parser.add_argument(
        "--output-csv",
        default=default_csv,
        help="Target CSV file for structured metric export (default: results_preliminar.csv)",
    )
    parser.add_argument(
        "--sample-interval",
        type=float,
        default=0.010,
        help="Sampling interval in seconds for /proc/status polling (default: 0.010 = 10ms)",
    )

    args = parser.parse_args()

    if not os.path.isfile(args.executable):
        print(f"[ERROR] Executable not found at '{args.executable}'. Please build the project first.")
        sys.exit(1)

    if not os.path.isfile(args.mesh_a):
        print(f"[ERROR] Mesh A not found at '{args.mesh_a}'.")
        sys.exit(1)

    if not os.path.isfile(args.mesh_b):
        print(f"[ERROR] Mesh B not found at '{args.mesh_b}'.")
        sys.exit(1)

    print("=" * 60)
    print("POMPEIU-HAUSDORFF MEMORY & TIME PROFILER")
    print(f"Target Binary: {args.executable}")
    print(f"Mesh A:        {args.mesh_a}")
    print(f"Mesh B:        {args.mesh_b}")
    print(f"Engine:        {args.engine}")
    print(f"Output CSV:    {args.output_csv}")
    print(f"Runs:          {args.runs}")
    print("=" * 60)

    for run_idx in range(1, args.runs + 1):
        if args.runs > 1:
            print(f"\n>>> Running iteration {run_idx}/{args.runs}...")
        record = run_isolated_profile(
            executable=args.executable,
            mesh_a=args.mesh_a,
            mesh_b=args.mesh_b,
            tol=args.tol,
            max_factor=args.max_factor,
            normalize=args.normalize,
            engine=args.engine,
            sample_interval_sec=args.sample_interval,
            verbose=True,
        )

        append_record_to_csv(args.output_csv, record)
        print(f"[OK] Record successfully appended to: {args.output_csv}")

    print("\n[SUCCESS] Profiling completed successfully.")


if __name__ == "__main__":
    main()
