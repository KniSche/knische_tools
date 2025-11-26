from setuptools import setup, find_packages

setup(
    name='knische_tools',              # Name of the package
    version='0.1.0',
    packages=find_packages(),          # Automatically finds the 'knische_tools' folder
    install_requires=[                 # Dependencies
        'scanpy',
        'matplotlib',
        'pandas',
        'numpy'
    ],
    author='Vincent Knight-Schrijver (KniSche)',
    description='Custom tools for bioinformatics',
)
