#!/usr/bin/env python
from setuptools import setup

setup(
    name="elgato",
    version="1.0.0",
    scripts = ['el_gato/el_gato.py'],
    install_requires=[
        "colorama; platform_system == 'Linux'",
        "importlib-metadata; python_version < '3.8'",
    ],
    packages = ["el_gato"],
    package_dir={"": "./"}
)
