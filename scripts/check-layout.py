"""Check repository paths after manual cleanup; does not move or delete files."""
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
ROOT_DIRECTORIES = {
    ".git", ".githooks", ".codex", ".agents", "0", "assets", "configs", "refs",
    "definitions", "docs", "meta", "out", "references", "runs", "scripts", "src",
    "tests", "tools", "vcpkg-ports", "saved_worlds",
}
ROOT_FILES = {
    ".gitignore", ".gitattributes", "AGENTS.md", "AGENTS.override.md",
    "README.md", "CMakeLists.txt", "CMakePresets.json", "vcpkg.json",
    "LICENSE", "version.txt",
}


def main():
    issues = []
    for path in ROOT.iterdir():
        allowed = ROOT_DIRECTORIES if path.is_dir() else ROOT_FILES
        if path.name not in allowed:
            issues.append(f"Unclassified root item: {path.name}; review for 0/ or its scope folder")

    for path in ((ROOT / "refs").iterdir() if (ROOT / "refs").exists() else ()):
        if path.name not in {"source", "processed", "prepare.py", "README.md", "TODO.txt"}:
            issues.append(f"Unclassified dataset root item: {path.name}; generated metadata belongs in processed/metadata/")

    cmake_files = [ROOT / "CMakeLists.txt"]
    for directory in ("src", "tests", "tools"):
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

    for path in (ROOT / "refs/processed").rglob("*.xml"):
        if path.name.endswith(".png.aux.xml"):
            if not path.with_name(path.name.removesuffix(".aux.xml")).is_file():
                issues.append(f"Orphan PNG georeferencing: {path.relative_to(ROOT)}")
        elif "<PAMDataset>" in path.read_text(encoding="utf-8", errors="replace"):
            issues.append(f"GDAL sidecar must retain the image filename and .aux.xml: {path.relative_to(ROOT)}")

    for category in ("maps", "fields"):
        directory = ROOT / "runs" / category
        for path in (directory.iterdir() if directory.exists() else ()):
            if path.is_dir() and path.name != "climate":
                issues.append(f"Unclassified run category: {path.relative_to(ROOT)}")
        for path in (ROOT / "runs" / category / "climate").rglob("*"):
            if path.is_file() and path.suffix in {".png", ".tif"}:
                if not path.stem.isdigit() or "0" in path.relative_to(ROOT / "runs" / category).parts[:-1]:
                    issues.append(f"Run maps require numeric IDs and 1-based seasons: {path.relative_to(ROOT)}")

    reference_maps = ROOT / "refs/processed/climate/maps"
    for path in reference_maps.rglob("*.png"):
        if not path.stem.isdigit() or path.parent == reference_maps:
            issues.append(f"Reference maps require semantic folders and width filenames: {path.relative_to(ROOT)}")

    if issues:
        print("\n".join(issues))
        return 1
    print("Layout checks passed: root items, CMake sources, reference sidecars, run categories.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
