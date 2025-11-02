from setuptools import setup, find_packages
import re

# version
def get_version():
    with open("smiles2mol/__init__.py") as f:
        match = re.search(r'__version__ = "([^"]+)"', f.read())
        if match:
            return match.group(1)
        raise RuntimeError("Version not found in smiles2mol/__init__.py")

# long description
with open("README.md", "r") as f:
    long_description = f.read()
    
# install requires
with open("requirements.txt", "r") as f:
    install_requires = f.readlines()

setup(
    name="smiles2mol",
    version=get_version(),
    description="Chemical structure generator from SMILES strings.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Shota Inoue",
    entry_points={
        "console_scripts": [
            "smiles2mol = smiles2mol.main:main",
        ],
    },
    url="https://github.com/s-inoue0108/smiles2mol",
    license="MIT",
    include_package_data=True,
    packages=find_packages(),
    install_requires=install_requires,
    python_requires=">= 3.10, < 3.12",
)

