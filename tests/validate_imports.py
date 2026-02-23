#!/usr/bin/env python3
"""Validate all Python script imports without executing pipeline logic."""
import importlib
import sys
import os
from pathlib import Path


def check_stdlib_or_installed(module_name):
    """Check if a module is importable."""
    try:
        importlib.import_module(module_name)
        return True, None
    except ImportError as e:
        return False, str(e)


def check_script_syntax(script_path):
    """Check script for syntax errors without executing."""
    try:
        with open(script_path) as f:
            source = f.read()
        compile(source, script_path, "exec")
        return True, None
    except SyntaxError as e:
        return False, f"Line {e.lineno}: {e.msg}"


def extract_imports(script_path):
    """Extract top-level import module names from a Python script."""
    imports = []
    with open(script_path) as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith("import ") or stripped.startswith("from "):
                # Skip conditional/indented imports
                if line[0] in (' ', '\t'):
                    continue
                # Extract module name
                if stripped.startswith("from "):
                    module = stripped.split()[1].split(".")[0]
                else:
                    module = stripped.split()[1].split(".")[0].rstrip(",")
                imports.append(module)
    return list(set(imports))


def main():
    base_dir = sys.argv[1] if len(sys.argv) > 1 else "."
    scripts_dir = os.path.join(base_dir, "scripts")

    # Find all Python scripts
    py_files = sorted(
        list(Path(base_dir).glob("*.py")) + list(Path(scripts_dir).glob("*.py"))
    )

    results = {"pass": [], "fail": [], "warn": []}

    # Required core libraries
    core_libs = ["pandas", "numpy", "scipy", "matplotlib", "seaborn", "Bio"]
    print("=== Core Library Check ===")
    for lib in core_libs:
        ok, err = check_stdlib_or_installed(lib)
        status = "OK" if ok else "MISSING"
        print(f"  {lib}: {status}")
        if not ok:
            results["fail"].append(f"Core lib {lib}: {err}")

    # ete3 special handling (Python 3.13+)
    print("\n=== ete3 Compatibility Check ===")
    try:
        import ete3
        print(f"  ete3: OK (version {ete3.__version__})")
    except ImportError:
        print("  ete3: Not importable (expected on Python 3.13+, pipeline has fallback)")
        results["warn"].append("ete3 not importable - verify conditional import fallbacks")

    # Check each Python script
    print(f"\n=== Script Syntax + Import Check ({len(py_files)} files) ===")
    for script in py_files:
        name = script.name

        # Syntax check
        ok, err = check_script_syntax(str(script))
        if not ok:
            print(f"  FAIL {name}: {err}")
            results["fail"].append(f"{name}: syntax error - {err}")
            continue

        # Import extraction
        imports = extract_imports(str(script))
        missing = []
        for imp in imports:
            if imp in ("__future__",):
                continue
            ok, _ = check_stdlib_or_installed(imp)
            if not ok:
                missing.append(imp)

        if missing:
            print(f"  WARN {name}: missing imports: {', '.join(missing)}")
            results["warn"].append(f"{name}: missing {', '.join(missing)}")
        else:
            print(f"  OK   {name} ({len(imports)} imports)")
            results["pass"].append(name)

    # Summary
    print(f"\n=== Summary ===")
    print(f"  PASS: {len(results['pass'])}")
    print(f"  WARN: {len(results['warn'])}")
    print(f"  FAIL: {len(results['fail'])}")

    if results["fail"]:
        print("\nFAILURES:")
        for f in results["fail"]:
            print(f"  - {f}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
