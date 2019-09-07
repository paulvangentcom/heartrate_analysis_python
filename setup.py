import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="heartpy",
    version="1.1.7",
    author="Paul van Gent",
    author_email="info@paulvangent.com",
    description="Heart Rate Analysis Toolkit",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/paulvangentcom/heartrate_analysis_python",
    packages=["heartpy", "heartpy.datasets"],
    include_package_data=True,
    package_data={'heartpy': ['datasets/data/*.csv']},
    install_requires=["numpy", "scipy", "matplotlib"],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: OS Independent",
    ],
)
