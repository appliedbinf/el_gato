#!/usr/bin/env python
from setuptools import setup

setup(
    name="elgato",
    version="1.10.3",
    python_requires='>=3.8',
    scripts = [
        'el_gato/el_gato.py',
        'run_el_gato.nf'
    ],
    install_requires=[
        "colorama; platform_system == 'Linux'",
        "importlib-metadata; python_version <= '3.8'",
    ],
    packages = ["el_gato"],
    package_dir={"": "./"}
)
