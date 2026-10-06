from setuptools import Extension, setup
import pybind11


setup(
    ext_modules=[
        Extension(
            "mappa._group_ops",
            sources=["src/mappa/native/group_ops.cpp"],
            include_dirs=[pybind11.get_include()],
            language="c++",
            extra_compile_args=["-O3", "-std=c++17"],
        ),
        Extension(
            "mappa._native",
            sources=["src/mappa/native/paired_mapping.cpp"],
            language="c++",
            extra_compile_args=["-O3", "-std=c++17"],
        )
    ]
)
