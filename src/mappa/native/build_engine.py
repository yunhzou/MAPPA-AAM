"""Build mappa._engine from this checkout.

Run: python src/mappa/native/build_engine.py

Requires pybind11 and a C++17 compiler. Vendored nauty lives separately under
third_party/nauty. MAPPA_NATIVE=0 explicitly selects the Python engine; a missing
selected native engine is an error.
"""
import shutil
from pathlib import Path

import pybind11
from setuptools import Distribution, Extension
from setuptools.command.build_ext import build_ext

NATIVE = Path(__file__).resolve().parent
ROOT = NATIVE.parents[2]
NAUTY = NATIVE / "third_party" / "nauty"
BUILD = ROOT / "build" / "engine"

nauty_sources = [str(NAUTY / name) for name in (
    "nauty.c", "nautil.c", "naugraph.c", "schreier.c", "naurng.c")]

ext = Extension(
    "mappa._engine",
    sources=[str(NATIVE / "engine.cpp"),
             str(NATIVE / "freeze.cpp"),
             str(NATIVE / "autgrp.cpp"),
             *nauty_sources],
    include_dirs=[pybind11.get_include(), str(NAUTY)],
    language="c++",
    extra_compile_args=["-O3", "-std=c++17", "-ffp-contract=off", "-fvisibility=hidden",
                        "-DMAXN=0", "-Wno-unused-result"],
)


class Build(build_ext):
    def build_extension(self, extension):
        # nauty is C; compile it with the C compiler flags, the .cpp files as C++
        c_flags = ["-O3", "-fPIC", "-DMAXN=0", "-Wno-unused-result"]
        objects = []
        compiler = self.compiler
        for source in extension.sources:
            is_c = source.endswith(".c")
            objects += compiler.compile(
                [source], output_dir=str(BUILD / "obj"),
                include_dirs=extension.include_dirs,
                extra_postargs=c_flags if is_c else extension.extra_compile_args)
        output = self.get_ext_fullpath(extension.name)
        Path(output).parent.mkdir(parents=True, exist_ok=True)
        compiler.link_shared_object(
            objects, output, extra_postargs=["-std=c++17"],
            target_lang="c++")


def main():
    dist = Distribution({"name": "mappa_engine", "ext_modules": [ext]})
    cmd = Build(dist)
    cmd.build_lib = str(BUILD / "lib")
    cmd.build_temp = str(BUILD / "tmp")
    cmd.ensure_finalized()
    cmd.run()
    built = next((BUILD / "lib" / "mappa").glob("_engine*.so"))
    target = ROOT / "src" / "mappa" / built.name
    shutil.copy2(built, target)
    print("built", target)


if __name__ == "__main__":
    main()
