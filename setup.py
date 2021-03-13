import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="heartpy",
    version="1.2.6",
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
        "pyparsing==2.4.7;python_version=='2.7'"],
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
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: OS Independent",
    ],
)
