from setuptools import setup, find_packages

def parse_requirements(filename):
    """Load requirements from a pip requirements file"""
    with open(filename, "r") as f:
        lines = f.readlines()
    reqs = [line.strip() for line in lines if line.strip() and not line.startswith("#")]
    return reqs

setup(
    name="DashingTurtle",
    version="0.1.3",
    author="J. White Bear",
    author_email="jwbear15@gmail.com",
    description="An applicaton for building structural landscapes of RNA sequence modifications with Nanopore data.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/jwbear/Dashing_Turtle.git",
    packages=find_packages(),
    install_requires=parse_requirements("requirements.txt"),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
    entry_points={
        "console_scripts": [
            "dt-cli = DashML.UI.DT_CLI:main",
            "dt-gui = DashML.GUI.DT:main",
        ],
    },
)
