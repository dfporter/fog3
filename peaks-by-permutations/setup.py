from setuptools import setup, find_packages
setup(
    # Application name:
    name="peaks-by-permutations",
    description="Find peaks by permutations of reads.",
    long_description=open("README.rst").read(),

    # Version number (initial):
    version="0.1.0",

    # Application author details:
    author="dfporter",
    author_email='dfporter@wisc.edu',

    # Packages
    packages=find_packages(),

    # Include additional files into the package
    include_package_data=True,
    url='http://wisc.edu',
    #
    # license="LICENSE.txt",

    # Dependent packages (distributions)
)
