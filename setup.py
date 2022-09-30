#!/usr/bin/env python
from setuptools import setup

setup(
    name="elgato",
    scripts = ['el_gato.py'],
    install_requires=[
        "colorama; platform_system == 'Linux'",
        "importlib-metadata; python_version < '3.8'",
        "argparse",
        "inspect",
        "logging",
        "multiprocessing",
        "os",
        "re",
        "subprocess",
        "sys",
        "shutil",
        "shlex",
        "time",
    ],
)
