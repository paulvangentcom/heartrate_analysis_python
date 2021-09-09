import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="heartpy",
    version="1.2.7",
    author="Paul van Gent",
    author_email="P.vanGent@tudelft.nl",
    description="Heart Rate Analysis Toolkit",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/paulvangentcom/heartrate_analysis_python",
    packages=["heartpy"],
    install_requires=[
        "cycler==0.10.0;python_version<='3.5'",
        "kiwisolver==1.1.0;python_version<='3.5'",
        "pyparsing==2.4.7;python_version=='2.7'",
        "numpy<=1.15;python_version=='2.7'",
        "numpy<=1.15;python_version=='3.4'",
        "numpy<=1.15;python_version=='3.5'",
        "numpy<=1.17;python_version=='3.6'",
        "numpy>=1.17;python_version>='3.7'",
        "scipy<=1.1.0;python_version=='2.7'",
        "scipy<=1.2.0;python_version=='3.4'",
        "scipy<=1.4.1;python_version=='3.5'",
        "scipy<=1.4.1;python_version=='3.6'",
        "scipy;python_version>='3.7'",
        "matplotlib==2.2.3;python_version=='2.7'",
        "matplotlib<=2.2.3;python_version=='3.4'",
        "matplotlib<=3.0.0;python_version=='3.5'",
        "matplotlib<=3.1.2;python_version=='3.6'",
        "matplotlib<=3.3.4;python_version=='3.7'",
        "matplotlib;python_version>='3.8'"],
    include_package_data=True,
    package_data={
        '': ['data/*.csv', 'data/*.mat', 'data/*.log']       
    },
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
