import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="boldigger2",
    version="1.0.0",
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
        "Bio>=1.7.0",
        "biopython>=1.79",
        "joblib>=1.1.0",
        "more_itertools>=10.2.0",
        "numpy>=1.26.4",
        "pandas>=2.2.2",
        "Requests>=2.32.2",
        "requests_html>=0.10.0",
        "tqdm>=4.66.4",
        "tqdm_joblib>=0.0.4",
        "urllib3>=1.26.12",
    ],
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
    entry_points={
        "console_scripts": [
            "boldigger2 = boldigger2.__main__:main",
        ]
    },
)
