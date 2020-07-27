import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="mass_spec_utils", # Replace with your own username
    version="0.0.5",
    author="Simon Rogers",
    author_email="simon.d.rogers@gmail.com",
    description="Some useful MS code",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/sdrogers/mass-spec-utils",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires = [
    	'pymzml',
    	'molmass',
    	'numpy',
    	'requests',
    ],
)
