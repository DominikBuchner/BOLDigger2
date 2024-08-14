import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="boldigger2",
    version="2.1.1",
    author="Dominik Buchner",
    author_email="dominik.buchner524@googlemail.com",
    description="An even better python package to query different databases of boldsystems.org",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/DominikBuchner/BOLDigger",
    packages=setuptools.find_packages(),
    license="MIT",
    install_requires=[
        "beautifulsoup4>=4.12.3",
        "biopython>=1.79",
        "joblib>=1.1.0",
        "more_itertools>=10.2.0",
        "numpy>=1.24.0, <2.0.0",
        "pandas>=2.2.2",
        "Requests>=2.32.2",
        "requests_html>=0.10.0",
        "tqdm>=4.66.4",
        "tqdm_joblib>=0.0.4",
        "urllib3>=1.26.12",
        "tables>=3.9.2",
        "html5lib>=1.1",
        "lxml>=4.9.1",
        "soupsieve>=2.5",
        "openpyxl>=3.1.1",
        "pyarrow>=11.0.0",
        "lxml_html_clean>=0.1.1",
        "free-proxy >= 1.1.1",
    ],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.10",
    entry_points={
        "console_scripts": [
            "boldigger2 = boldigger2.__main__:main",
        ]
    },
)
