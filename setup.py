from setuptools import setup, find_packages
import os

setup(
    name = 'SWAMP-E',
    version = '0.0.1',
    description = '2D Shallow-Water General Circulation Model for Exoplanet Atmospheres',
    long_description = open(os.path.join(
                            os.path.dirname(__file__), 'README.rst')).read(),
    long_description_content_type = 'text/x-rst',
    author = 'Ekaterina Landgren',
    author_email = 'ek672@cornell.edu',
    license = 'MIT License',
    packages = ['SWAMP-E'],
    include_package_data = True,
    python_requires = '<3.10',
    install_requires = ['numpy',
                        'scipy',
                        'pickle',
                        'imageio',
                        'matplotlib<=3.5.1',
                        'astropy',
                        'pandas',
                        'jupyter'],
    zip_safe = False,
)