from setuptools import setup, find_packages
import nenums

setup(name       = 'nenums',
    packages     = find_packages(),
    include_package_data=True,
    install_requires=['numpy', 'astropy'],
    scripts      = ['bin/nenums'],
    version      = nenums.__version__,
    description  = 'NenuFAR XST to MS converter',
    url          = 'https://github.com/AlanLoh/nenums.git',
    author       = 'Alan Loh',
    author_email = 'alan.loh@obspm.fr',
    license      = 'MIT',
    zip_safe     = False)
