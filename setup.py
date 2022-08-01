from setuptools import setup, find_packages
import os

setup(
    name = 'SWAMPE',
    version = '0.0.23',
    description = '2D Shallow-Water General Circulation Model for Exoplanet Atmospheres',
    long_description = open(os.path.join(
                            os.path.dirname(__file__), 'README.md')).read(),
    long_description_content_type = 'text/x-rst',
    author = 'Ekaterina Landgren',
    author_email = 'ek672@cornell.edu',
    license = 'MIT License',
    packages = ['SWAMPE'],
    include_package_data = True,
    python_requires = '<3.10',
    install_requires = ['numpy',
                        'scipy',
                        'imageio',
                        'matplotlib',
                        'astropy',
                        'pandas',
                        'jupyter'],
    zip_safe = False,
)