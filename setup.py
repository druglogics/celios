from setuptools import setup, find_packages

setup(
	name="celios",
	version="0.1.0",
	description="CELIOS: a DrugLogics pipeline module that processes CELl LIne OmicS (CELIOS) for calibration of Boolean models in the DrugLogics and TRAFIKK pipeline.",
	author='Viviam Solangeli Bermudez',
	author_email='viviamsb@ntnu.no',
	packages=find_packages(where="src", exclude=("tests",)),
	package_dir={"": "src"},
	include_package_data=True,
	install_requires=[
		"pandas>=1.0.0",
	],
	extras_require={
		"dev": [
			"pytest>=6.0.0",
		],
	},
	entry_points={
		"console_scripts": [
			"celios=celios.cli:main",
		]
	},
)
