from setuptools import setup

setup(
    name="elgato",
    install_requires=[
        "colorama; platform_system == 'Linux'",
        "importlib-metadata; python_version < '3.8'",
    ],
)
