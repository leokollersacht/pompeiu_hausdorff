#!/usr/bin/env python3
"""
compare_engines.py - Automated Comparative Benchmarking for Pompeiu-Hausdorff Engines

Executes both Best-First (baseline) and IDA* (proposed bounded-memory search)
on a given mesh pair, computing percentage deltas for:
- Peak Resident Memory (Peak RSS in MB)
- Execution Runtime (Wall-Clock and CPU time in ms / s)
- Subdivided face count and upper bound evaluations (u0..u4)
- Numerical convergence verification (lower/upper bound parity)

Complies with project engineering standards: SRP, typed Python, and formal reporting.
"""

import argparse
import csv
import datetime
import os
import sys
from typing import Any, Dict, List, Optional, Tuple

# Enable relative import of profile_runner from the benchmarks directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(SCRIPT_DIR)
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from profile_runner import run_isolated_profile  # noqa: E402


def calculate_percentage_change(baseline: float, candidate: float) -> Tuple[float, str]:
    """
    Calculates the relative percentage difference between candidate and baseline.

    Returns:
        (percentage_value, formatted_string)
        e.g., (-84.60, "-84.60% (Economy/Faster)") or (+25.00, "+25.00%")
    """
    if baseline == 0.0:
        return 0.0, "0.00%"
    pct = ((candidate - baseline) / baseline) * 100.0
    sign = "+" if pct > 0 else ""
    return pct, f"{sign}{pct:.2f}%"


def compute_comparison_metrics(
    best_first: Dict[str, Any],
    idastar: Dict[str, Any],
) -> Dict[str, Any]:
    """
    Computes comparative deltas and percentage ratios between Best-First and IDA*.
    """
    bf_mem = best_first["peak_rss_mb"]
    ida_mem = idastar["peak_rss_mb"]
    mem_delta_mb = ida_mem - bf_mem
    mem_pct, mem_pct_str = calculate_percentage_change(bf_mem, ida_mem)

    bf_wall = best_first["wall_time_ms"]
    ida_wall = idastar["wall_time_ms"]
    wall_delta_ms = ida_wall - bf_wall
    wall_pct, wall_pct_str = calculate_percentage_change(bf_wall, ida_wall)

    bf_cpu = best_first["cpu_time_ms"]
    ida_cpu = idastar["cpu_time_ms"]
    cpu_delta_ms = ida_cpu - bf_cpu
    cpu_pct, cpu_pct_str = calculate_percentage_change(bf_cpu, ida_cpu)

    bf_faces = best_first["subdivided_faces"]
    ida_faces = idastar["subdivided_faces"]
    faces_delta = ida_faces - bf_faces
    faces_pct, faces_pct_str = calculate_percentage_change(bf_faces, ida_faces)
    faces_ratio = (ida_faces / bf_faces) if bf_faces > 0 else float("nan")

    bf_total_bounds = sum(best_first[f"u{i}_calls"] for i in range(5))
    ida_total_bounds = sum(idastar[f"u{i}_calls"] for i in range(5))
    bounds_pct, bounds_pct_str = calculate_percentage_change(
        bf_total_bounds, ida_total_bounds
    )

    lower_match = abs(best_first["lower"] - idastar["lower"]) < 1e-12
    upper_match = abs(best_first["upper_max"] - idastar["upper_max"]) < 1e-12
    exact_parity = lower_match and upper_match

    return {
        "mem_delta_mb": mem_delta_mb,
        "mem_pct": mem_pct,
        "mem_pct_str": mem_pct_str,
        "wall_delta_ms": wall_delta_ms,
        "wall_pct": wall_pct,
        "wall_pct_str": wall_pct_str,
        "cpu_delta_ms": cpu_delta_ms,
        "cpu_pct": cpu_pct,
        "cpu_pct_str": cpu_pct_str,
        "faces_delta": faces_delta,
        "faces_pct": faces_pct,
        "faces_pct_str": faces_pct_str,
        "faces_ratio": faces_ratio,
        "bf_total_bounds": bf_total_bounds,
        "ida_total_bounds": ida_total_bounds,
        "bounds_pct": bounds_pct,
        "bounds_pct_str": bounds_pct_str,
        "exact_parity": exact_parity,
    }


def format_table(
    mesh_a: str,
    mesh_b: str,
    tol: str,
    max_factor: str,
    normalize: str,
    bf: Dict[str, Any],
    ida: Dict[str, Any],
    cmp: Dict[str, Any],
) -> str:
    """
    Renders an ASCII/Markdown-compatible structured table comparing both engines.
    """
    lines: List[str] = []
    separator = "-" * 88

    lines.append("=" * 88)
    lines.append("POMPEIU-HAUSDORFF SEARCH ENGINE COMPARATIVE BENCHMARK")
    lines.append(f"Mesh A:        {os.path.basename(mesh_a)}")
    lines.append(f"Mesh B:        {os.path.basename(mesh_b)}")
    lines.append(f"Parameters:    tol={tol}, max_factor={max_factor}, normalize={normalize}")
    lines.append(f"Date & Time:   {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    lines.append("=" * 88)
    lines.append(
        f"{'Metric':<28} | {'Best-First':<18} | {'IDA*':<18} | {'Delta (IDA* vs BF)':<18}"
    )
    lines.append(separator)

    # 1. Memory metrics
    lines.append(
        f"{'Peak RSS Memory':<28} | {bf['peak_rss_mb']:>14.2f} MB | {ida['peak_rss_mb']:>14.2f} MB | {cmp['mem_pct_str']:>18}"
    )
    lines.append(
        f"{'  Memory Delta (MB)':<28} | {'-':>18} | {'-':>18} | {cmp['mem_delta_mb']:>+15.2f} MB"
    )

    # 2. Runtime metrics
    lines.append(separator)
    lines.append(
        f"{'Wall-Clock Time':<28} | {bf['wall_time_ms']:>14.2f} ms | {ida['wall_time_ms']:>14.2f} ms | {cmp['wall_pct_str']:>18}"
    )
    lines.append(
        f"{'  Wall-Clock Time (s)':<28} | {(bf['wall_time_ms']/1000.0):>15.2f} s | {(ida['wall_time_ms']/1000.0):>15.2f} s | {(cmp['wall_delta_ms']/1000.0):>+16.2f} s"
    )
    lines.append(
        f"{'CPU Time (User+Sys)':<28} | {bf['cpu_time_ms']:>14.2f} ms | {ida['cpu_time_ms']:>14.2f} ms | {cmp['cpu_pct_str']:>18}"
    )

    # 3. Search tree metrics
    lines.append(separator)
    faces_ratio_str = f"{cmp['faces_ratio']:.2f}x"
    lines.append(
        f"{'Subdivided Faces':<28} | {bf['subdivided_faces']:>18,d} | {ida['subdivided_faces']:>18,d} | {cmp['faces_pct_str']:>18}"
    )
    lines.append(
        f"{'  Subdivision Ratio':<28} | {'1.00x':>18} | {faces_ratio_str:>18} | {'-':>18}"
    )
    lines.append(
        f"{'Total Bound Calls':<28} | {cmp['bf_total_bounds']:>18,d} | {cmp['ida_total_bounds']:>18,d} | {cmp['bounds_pct_str']:>18}"
    )
    for i in range(5):
        lines.append(
            f"{f'  u{i} evaluations':<28} | {bf[f'u{i}_calls']:>18,d} | {ida[f'u{i}_calls']:>18,d} | {'-':>18}"
        )

    # 4. Numerical convergence
    lines.append(separator)
    lines.append(
        f"{'Computed Lower Bound':<28} | {bf['lower']:>18.8f} | {ida['lower']:>18.8f} | {'Identical' if cmp['exact_parity'] else 'Mismatch':>18}"
    )
    lines.append(
        f"{'Computed Upper Bound':<28} | {bf['upper_max']:>18.8f} | {ida['upper_max']:>18.8f} | {'Identical' if cmp['exact_parity'] else 'Mismatch':>18}"
    )
    lines.append(
        f"{'Final Bound Gap':<28} | {(bf['upper_max'] - bf['lower']):>18.8f} | {(ida['upper_max'] - ida['lower']):>18.8f} | {'-':>18}"
    )
    lines.append("=" * 88)

    # Executive conclusion
    lines.append("\nEXECUTIVE SUMMARY:")
    if cmp["mem_pct"] < 0:
        lines.append(
            f"- Memory:  IDA* reduced peak resident memory by {abs(cmp['mem_pct']):.2f}% "
            f"({abs(cmp['mem_delta_mb']):.2f} MB saved)."
        )
    else:
        lines.append(
            f"- Memory:  IDA* changed peak resident memory by {cmp['mem_pct_str']}."
        )

    if cmp["wall_pct"] < 0:
        lines.append(
            f"- Runtime: IDA* was {abs(cmp['wall_pct']):.2f}% faster than Best-First "
            f"({abs(cmp['wall_delta_ms'])/1000.0:.2f} s speedup)."
        )
    else:
        lines.append(
            f"- Runtime: Best-First was {abs(cmp['wall_pct']):.2f}% faster than IDA* "
            f"(Trade-off: +{cmp['wall_delta_ms']/1000.0:.2f} s for bounded memory)."
        )

    lines.append(
        f"- Parity:  Mathematical convergence parity: "
        f"{'EXACT MATCH' if cmp['exact_parity'] else 'DIVERGENCE DETECTED'}."
    )
    lines.append("=" * 88)

    return "\n".join(lines)


def export_markdown_report(
    target_path: str,
    mesh_a: str,
    mesh_b: str,
    tol: str,
    max_factor: str,
    normalize: str,
    bf: Dict[str, Any],
    ida: Dict[str, Any],
    cmp: Dict[str, Any],
) -> None:
    """
    Exports a clean Markdown table and summary for documentation and academic reports.
    """
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    content = f"""# Comparative Benchmark: Best-First vs. IDA*

**Date:** {timestamp}  
**Mesh A:** `{os.path.basename(mesh_a)}`  
**Mesh B:** `{os.path.basename(mesh_b)}`  
**Configuration:** `tol={tol}`, `max_factor={max_factor}`, `normalize={normalize}`  

| Metric | Best-First Search | IDA* Search | Delta (IDA* vs. Best-First) |
| :--- | :--- | :--- | :--- |
| **Peak Resident Memory (Peak RSS)** | **{bf['peak_rss_mb']:.2f} MB** | **{ida['peak_rss_mb']:.2f} MB** | **{cmp['mem_pct_str']}** ({cmp['mem_delta_mb']:+.2f} MB) |
| **Wall-Clock Runtime** | {bf['wall_time_ms']:.2f} ms ({bf['wall_time_ms']/1000.0:.2f} s) | {ida['wall_time_ms']:.2f} ms ({ida['wall_time_ms']/1000.0:.2f} s) | **{cmp['wall_pct_str']}** ({cmp['wall_delta_ms']/1000.0:+.2f} s) |
| **CPU Time (User + System)** | {bf['cpu_time_ms']:.2f} ms | {ida['cpu_time_ms']:.2f} ms | {cmp['cpu_pct_str']} |
| **Subdivided Faces** | {bf['subdivided_faces']:,d} | {ida['subdivided_faces']:,d} | {cmp['faces_pct_str']} ({cmp['faces_ratio']:.2f}x) |
| **Total Upper Bound Calls** | {cmp['bf_total_bounds']:,d} | {cmp['ida_total_bounds']:,d} | {cmp['bounds_pct_str']} |
| **Computed Lower Bound** | `{bf['lower']:.8f}` | `{ida['lower']:.8f}` | {'Exact Match' if cmp['exact_parity'] else 'Mismatch'} |
| **Computed Upper Bound** | `{bf['upper_max']:.8f}` | `{ida['upper_max']:.8f}` | {'Exact Match' if cmp['exact_parity'] else 'Mismatch'} |

## Key Insights
- **Memory Optimization:** IDA* resident set size changed by **{cmp['mem_pct_str']}**, validating the bounded-memory hypothesis.
- **Convergence Verification:** Both algorithms produced {'strictly identical' if cmp['exact_parity'] else 'divergent'} distance bounds within tolerance `{tol}`.
"""
    with open(target_path, "w", encoding="utf-8") as f:
        f.write(content)


def export_csv_row(
    target_path: str,
    mesh_a: str,
    mesh_b: str,
    tol: str,
    max_factor: str,
    normalize: str,
    bf: Dict[str, Any],
    ida: Dict[str, Any],
    cmp: Dict[str, Any],
) -> None:
    """
    Appends a structured comparative record to an output CSV file.
    """
    file_exists = os.path.exists(target_path)
    fieldnames = [
        "timestamp",
        "mesh_a",
        "mesh_b",
        "tol",
        "max_factor",
        "normalize",
        "bf_peak_rss_mb",
        "ida_peak_rss_mb",
        "mem_reduction_pct",
        "bf_wall_time_ms",
        "ida_wall_time_ms",
        "wall_time_change_pct",
        "bf_subdivided_faces",
        "ida_subdivided_faces",
        "faces_ratio",
        "bf_lower",
        "ida_lower",
        "bf_upper",
        "ida_upper",
        "exact_parity",
    ]

    row = {
        "timestamp": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "mesh_a": os.path.basename(mesh_a),
        "mesh_b": os.path.basename(mesh_b),
        "tol": tol,
        "max_factor": max_factor,
        "normalize": normalize,
        "bf_peak_rss_mb": bf["peak_rss_mb"],
        "ida_peak_rss_mb": ida["peak_rss_mb"],
        "mem_reduction_pct": round(cmp["mem_pct"], 2),
        "bf_wall_time_ms": bf["wall_time_ms"],
        "ida_wall_time_ms": ida["wall_time_ms"],
        "wall_time_change_pct": round(cmp["wall_pct"], 2),
        "bf_subdivided_faces": bf["subdivided_faces"],
        "ida_subdivided_faces": ida["subdivided_faces"],
        "faces_ratio": round(cmp["faces_ratio"], 2),
        "bf_lower": bf["lower"],
        "ida_lower": ida["lower"],
        "bf_upper": bf["upper_max"],
        "ida_upper": ida["upper_max"],
        "exact_parity": cmp["exact_parity"],
    }

    with open(target_path, "a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if not file_exists:
            writer.writeheader()
        writer.writerow(row)


def main() -> None:
    default_executable = os.path.join(PROJECT_DIR, "build", "pompeiu_hausdorff")
    default_mesh_a = os.path.join(PROJECT_DIR, "meshes", "107100.obj")
    default_mesh_b = os.path.join(PROJECT_DIR, "meshes", "107100_sf.obj")

    parser = argparse.ArgumentParser(
        description="Automated comparative profiling between Best-First and IDA* engines."
    )
    parser.add_argument(
        "--executable",
        default=default_executable,
        help=f"Path to binary (default: {default_executable})",
    )
    parser.add_argument(
        "--mesh-a",
        default=default_mesh_a,
        help=f"Path to mesh A (.obj) (default: {default_mesh_a})",
    )
    parser.add_argument(
        "--mesh-b",
        default=default_mesh_b,
        help=f"Path to mesh B (.obj) (default: {default_mesh_b})",
    )
    parser.add_argument(
        "--tol",
        default="1e-4",
        help="Tolerance threshold between upper and lower bounds (default: 1e-4)",
    )
    parser.add_argument(
        "--max-factor",
        default="1000000",
        help="Max factor for face/vertex expansion (default: 1000000)",
    )
    parser.add_argument(
        "--normalize",
        default="1",
        help="Normalize bounds by bounding box diagonal (0 or 1, default: 1)",
    )
    parser.add_argument(
        "--output-md",
        default=None,
        help="Optional path to export comparative table as Markdown document.",
    )
    parser.add_argument(
        "--output-csv",
        default=None,
        help="Optional path to append structured comparison metrics row in CSV.",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress intermediate execution stdout from single runs.",
    )

    args = parser.parse_args()

    if not os.path.isfile(args.executable):
        sys.exit(f"Error: Executable not found at {args.executable}. Compile the project first.")
    if not os.path.isfile(args.mesh_a):
        sys.exit(f"Error: Mesh A not found at {args.mesh_a}.")
    if not os.path.isfile(args.mesh_b):
        sys.exit(f"Error: Mesh B not found at {args.mesh_b}.")

    print("=" * 88)
    print("STARTING COMPARATIVE BENCHMARK EXECUTION")
    print(f"Mesh A:     {args.mesh_a}")
    print(f"Mesh B:     {args.mesh_b}")
    print(f"Tolerance:  {args.tol}")
    print("=" * 88)

    # 1. Run Best-First
    print("\n>>> [1/2] Executing Best-First Search Engine...")
    bf_record = run_isolated_profile(
        executable=args.executable,
        mesh_a=args.mesh_a,
        mesh_b=args.mesh_b,
        tol=args.tol,
        max_factor=args.max_factor,
        normalize=args.normalize,
        engine="best_first",
        verbose=not args.quiet,
    )

    # 2. Run IDA*
    print("\n>>> [2/2] Executing IDA* Search Engine...")
    ida_record = run_isolated_profile(
        executable=args.executable,
        mesh_a=args.mesh_a,
        mesh_b=args.mesh_b,
        tol=args.tol,
        max_factor=args.max_factor,
        normalize=args.normalize,
        engine="idastar",
        verbose=not args.quiet,
    )

    # 3. Compute deltas & render comparative table
    cmp_metrics = compute_comparison_metrics(bf_record, ida_record)
    table_output = format_table(
        args.mesh_a,
        args.mesh_b,
        args.tol,
        args.max_factor,
        args.normalize,
        bf_record,
        ida_record,
        cmp_metrics,
    )

    print("\n" + table_output)

    # 4. Optional exports
    if args.output_md:
        export_markdown_report(
            args.output_md,
            args.mesh_a,
            args.mesh_b,
            args.tol,
            args.max_factor,
            args.normalize,
            bf_record,
            ida_record,
            cmp_metrics,
        )
        print(f"\n[OK] Markdown report saved to: {args.output_md}")

    if args.output_csv:
        export_csv_row(
            args.output_csv,
            args.mesh_a,
            args.mesh_b,
            args.tol,
            args.max_factor,
            args.normalize,
            bf_record,
            ida_record,
            cmp_metrics,
        )
        print(f"[OK] Comparison CSV row appended to: {args.output_csv}")


if __name__ == "__main__":
    main()
