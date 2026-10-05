"""Check repository paths after manual cleanup; does not move or delete files."""
from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from uwtools import palettes

ROOT = Path(__file__).resolve().parents[1]
ROOT_DIRECTORIES = {
    ".git", ".githooks", ".codex", ".agents", "assets", "configs", "refs",
    "definitions", "docs", "meta", "out", "references", "runs", "scripts", "src",
    "tests", "vcpkg-ports", "saved_worlds", ".venv", ".claude",
}
# Code outside the tectonics engine includes only these engine headers.
TECTONICS_ENGINE = ROOT / "src/simulations/geology/plate_tectonics"
TECTONICS_PUBLIC_HEADERS = {
    "tectonic_api.hpp", "tectonic_contract.hpp", "tectonic_recipe.hpp", "tectonic_scenario.hpp",
    "tectonic_world_grid.hpp", "topography_codec.hpp", "material.hpp", "utils.hpp",
}
RUNS_RESERVED = {"archive", "ledger", "registry", "reports", "work"}
ROOT_FILES = {
    ".gitignore", ".gitattributes", "AGENTS.md", "AGENTS.override.md",
    "README.md", "CMakeLists.txt", "CMakePresets.json", "vcpkg.json",
    "LICENSE", "version.txt", "pyproject.toml", "uv.lock",
}


def main():
    issues = []
    for path in ROOT.iterdir():
        allowed = ROOT_DIRECTORIES if path.is_dir() else ROOT_FILES
        if path.name not in allowed:
            issues.append(f"Unclassified root item: {path.name}; move it to its scope folder")

    for path in ((ROOT / "refs").iterdir() if (ROOT / "refs").exists() else ()):
        if path.name not in {"source", "processed", "README.md"}:
            issues.append(f"Unclassified dataset root item: {path.name}; generated metadata belongs in processed/metadata/")

    cmake_files = [ROOT / "CMakeLists.txt"]
    for directory in ("src", "tests"):
        cmake_files.extend((ROOT / directory).rglob("CMakeLists.txt"))
        cmake_files.extend((ROOT / directory).rglob("*.cmake"))
    for cmake_file in cmake_files:
        cmake = cmake_file.read_text(encoding="utf-8")
        # CMake itself checks expanded variables. Check literal module sources
        # here, including the root-relative world_modules.cmake include.
        for path in sorted(set(re.findall(r"(?<![\w${}])([\w./-]+\.(?:cpp|cu))\b", cmake))):
            base = ROOT if path.startswith(("src/", "tests/")) else cmake_file.parent
            if not (base / path).is_file():
                issues.append(f"Missing CMake source in {cmake_file.relative_to(ROOT)}: {path}")

    internal = {path.name for path in TECTONICS_ENGINE.glob("*.hpp")} - TECTONICS_PUBLIC_HEADERS
    for directory in ("src",):
        for path in (ROOT / directory).rglob("*.[ch]pp"):
            relative = path.relative_to(ROOT).as_posix()
            if TECTONICS_ENGINE in path.parents:
                continue
            for header in re.findall(r'#include\s+"([^"]+)"', path.read_text(encoding="utf-8", errors="replace")):
                if header in internal:
                    issues.append(f"Internal tectonics header {header} used outside the engine: {relative}")

    for path in (ROOT / "refs/processed").rglob("*.xml"):
        if path.name.endswith(".png.aux.xml"):
            if not path.with_name(path.name.removesuffix(".aux.xml")).is_file():
                issues.append(f"Orphan PNG georeferencing: {path.relative_to(ROOT)}")
        elif "<PAMDataset>" in path.read_text(encoding="utf-8", errors="replace"):
            issues.append(f"GDAL sidecar must retain the image filename and .aux.xml: {path.relative_to(ROOT)}")

    # Run layout (docs/architecture/consolidation-plan.md, section 4): runs/<module>/<id>/run.json;
    # the sealed layout under runs/archive/ is not checked.
    for module in ((ROOT / "runs").iterdir() if (ROOT / "runs").exists() else ()):
        if not module.is_dir() or module.name in RUNS_RESERVED:
            continue
        for path in module.iterdir():
            if not path.name.isdigit() or not path.is_dir():
                issues.append(f"Run directories are numeric ids: {path.relative_to(ROOT).as_posix()}")
            elif not (path / "run.json").is_file():
                issues.append(f"Run directory without run.json: {path.relative_to(ROOT).as_posix()}")
            else:
                # Every maps folder of a run names the palettes it was rendered with.
                for maps in sorted(folder for folder in path.rglob("maps") if folder.is_dir()):
                    problem = palettes.sidecar_issue(maps)
                    if problem:
                        issues.append(f"{maps.relative_to(ROOT).as_posix()}: {problem}")
    # Reference maps rendered before the sidecar existed have none; one that exists must be sound.
    for sidecar in (ROOT / "refs/processed").rglob(f"maps/{palettes.SIDECAR}"):
        problem = palettes.sidecar_issue(sidecar.parent)
        if problem:
            issues.append(f"{sidecar.parent.relative_to(ROOT).as_posix()}: {problem}")

    reference_maps = ROOT / "refs/processed/climate/maps"
    for path in reference_maps.rglob("*.png"):
        if not path.stem.isdigit() or path.parent == reference_maps:
            issues.append(f"Reference maps require semantic folders and width filenames: {path.relative_to(ROOT)}")

    if issues:
        print("\n".join(issues))
        return 1
    print("Layout checks passed: root items, CMake sources, tectonics boundary, reference sidecars, run directories, palette sidecars.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
