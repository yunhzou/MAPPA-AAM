from setuptools import Extension, setup
import pybind11


setup(
    ext_modules=[
        Extension(
            "mappa._group_ops",
            sources=["src/mappa/mappa_core/native/group_ops.cpp"],
            depends=[f"src/mappa/mappa_core/native/{name}.h" for name in (
                "set_orbit", "assignment_search", "permutation_group",
                "action_program", "event_response", "history_projection")],
            include_dirs=[pybind11.get_include()],
            language="c++",
            extra_compile_args=["-O3", "-std=c++17"],
        ),
    ]
)
